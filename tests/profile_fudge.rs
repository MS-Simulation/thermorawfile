//! The per-chunk profile `fudge` is an additive m/z correction applied after the
//! frequency→m/z conversion (OpenTFRaw 420dcef). Checked on the Orbitrap Velos sample,
//! where every MS1 profile chunk carries one (up to ~1.2e-3 m/z).

use thermorawfile::{Calibration, RawFile};

const DATA: &str = "tests/data/small2.RAW";

fn ms1_scans(rf: &RawFile) -> Vec<u32> {
    (rf.first_scan..=rf.last_scan)
        .filter(|&s| rf.scan_event(s).map(|e| e.ms_order) == Some(1) && rf.profile(s).is_some())
        .collect()
}

fn calibration(rf: &RawFile, scan: u32) -> Calibration {
    rf.calibration_at_event(rf.scan_event_byte_offset(scan).unwrap())
        .unwrap()
}

/// The centroid list the instrument stored after the profile, sorted by m/z.
fn stored_centroids(rf: &RawFile, scan: u32) -> Vec<f64> {
    let e = &rf.index[(scan - rf.first_scan) as usize];
    let pkt = (rf.data_addr + e.offset) as usize;
    let u = |o: usize| u32::from_le_bytes(rf.bytes[o..o + 4].try_into().unwrap());
    let (profile_words, peaklist_words) = (u(pkt + 4) as usize, u(pkt + 8));
    let pl = pkt + 40 + profile_words * 4;
    let n = u(pl) as usize;
    let wide = peaklist_words as usize == 1 + 3 * n;
    let mut out: Vec<f64> = (0..n)
        .map(|i| {
            if wide {
                f64::from_le_bytes(rf.bytes[pl + 4 + i * 12..pl + 12 + i * 12].try_into().unwrap())
            } else {
                f32::from_le_bytes(rf.bytes[pl + 4 + i * 8..pl + 8 + i * 8].try_into().unwrap()) as f64
            }
        })
        .collect();
    out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    out
}

fn nearest(sorted: &[f64], mz: f64) -> f64 {
    let i = sorted.partition_point(|&c| c < mz);
    [i.wrapping_sub(1), i]
        .iter()
        .filter_map(|&k| sorted.get(k))
        .map(|&c| mz - c)
        .min_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap())
        .unwrap()
}

#[test]
fn points_apply_fudge_after_conversion() {
    let rf = RawFile::open(DATA).unwrap();
    let (mut with, mut without, mut n) = (0.0f64, 0.0f64, 0u32);
    for scan in ms1_scans(&rf) {
        let cal = calibration(&rf, scan);
        let prof = rf.profile(scan).unwrap();
        assert!(prof.chunks.iter().all(|c| c.fudge != 0.0), "scan {scan}: expected fudge on every chunk");
        let cents = stored_centroids(&rf, scan);
        let pts = prof.points(&cal);
        // Parabolic apex of each clear local maximum, compared with the nearest stored centroid.
        let mut k = 0;
        for ch in &prof.chunks {
            let s = &ch.signal;
            for j in 1..s.len().saturating_sub(1) {
                if s[j] > s[j - 1] && s[j] >= s[j + 1] && s[j] > 1000.0 {
                    let (y0, y1, y2) = (s[j - 1] as f64, s[j] as f64, s[j + 1] as f64);
                    let d = 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2);
                    let bin = ch.first_bin as f64 + j as f64 + d;
                    let grid = cal.mz(prof.first_value + bin * prof.step);
                    let r0 = nearest(&cents, grid);
                    if r0.abs() < 0.01 {
                        without += r0.abs();
                        with += nearest(&cents, grid + ch.fudge as f64).abs();
                        n += 1;
                    }
                }
            }
            // points() is grid m/z + this chunk's fudge, in order.
            for (j, &v) in s.iter().enumerate() {
                let (mz, iv) = pts[k + j];
                assert_eq!(iv, v);
                assert_eq!(mz, prof.mz_of_bin(ch.first_bin + j as u32, &cal) + ch.fudge as f64);
            }
            k += s.len();
        }
    }
    let (with, without) = (with / n as f64, without / n as f64);
    assert!(n > 1000, "too few apexes: {n}");
    // Measured: 1.2e-5 with the fudge, 2.8e-4 without.
    assert!(with < 3e-5 && with * 10.0 < without, "mean |apex - centroid|: with {with:e}, without {without:e}");
}

#[test]
fn overlay_keeps_real_signal_in_place() {
    let mut rf = RawFile::open(DATA).unwrap();
    let scan = ms1_scans(&rf)[0];
    let cal = calibration(&rf, scan);
    let before = rf.profile(scan).unwrap().points(&cal);
    rf.overlay_profile(scan, &[], &cal).unwrap();
    let after = rf.profile(scan).unwrap().points(&cal);
    assert_eq!(before, after, "overlay with no sim peaks moved the real profile");
}

#[test]
fn overlay_sim_peak_reads_back_at_its_mz() {
    let mut rf = RawFile::open(DATA).unwrap();
    let scan = ms1_scans(&rf)[0];
    let cal = calibration(&rf, scan);
    let prof = rf.profile(scan).unwrap();
    let half_bin = |mz: f64| {
        let f = cal.freq(mz).unwrap();
        0.5 * (cal.mz(f) - cal.mz(f + prof.step)).abs()
    };

    // One sim peak right after the most fudged real chunk (it joins that chunk), one
    // between two real chunks far from both (it gets its own chunk with fudge 0).
    let ch = prof
        .chunks
        .iter()
        .max_by(|a, b| a.fudge.abs().partial_cmp(&b.fudge.abs()).unwrap())
        .unwrap();
    let edge = ch.first_bin + ch.signal.len() as u32;
    let joined = prof.mz_of_bin(edge, &cal) + ch.fudge as f64;
    let gap = prof
        .chunks
        .windows(2)
        .find(|w| w[1].first_bin > w[0].first_bin + w[0].signal.len() as u32 + 20)
        .unwrap();
    let alone_bin = (gap[0].first_bin + gap[0].signal.len() as u32 + gap[1].first_bin) / 2;
    let alone = prof.mz_of_bin(alone_bin, &cal);

    let sentinel = 7.5e6f32;
    rf.overlay_profile(scan, &[(joined, sentinel), (alone, sentinel)], &cal)
        .unwrap();
    let after = rf.profile(scan).unwrap();
    let pts = after.points(&cal);
    for target in [joined, alone] {
        let (mz, _) = pts
            .iter()
            .filter(|p| p.1 >= sentinel)
            .min_by(|a, b| (a.0 - target).abs().partial_cmp(&(b.0 - target).abs()).unwrap())
            .unwrap();
        assert!((mz - target).abs() <= half_bin(target) + 1e-9, "sim peak at {target} read back at {mz}");
    }
    // The joined peak sits in a chunk with the real chunk's fudge; the lone one in a fudge-0 chunk.
    let fudge_at = |b: u32| {
        after
            .chunks
            .iter()
            .find(|c| c.first_bin <= b && b < c.first_bin + c.signal.len() as u32)
            .map(|c| c.fudge)
    };
    assert_eq!(fudge_at(edge), Some(ch.fudge));
    assert_eq!(fudge_at(alone_bin), Some(0.0));
}
