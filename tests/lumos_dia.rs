//! Scan-event alignment on a Fusion Lumos DIA file: 232-byte MS1 and 288-byte MS2 events.
//!
//! Original RAW from the public PRIDE study PXD031322 (Fusion Lumos, gas-phase-fractionated
//! DIA). It is not redistributed here; download it and point `THERMORAWFILE_LUMOS_RAW` at it:
//!
//! ```text
//! curl -O https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/OFL001513-YLL-GPF-15K-1.raw
//! THERMORAWFILE_LUMOS_RAW=$PWD/OFL001513-YLL-GPF-15K-1.raw cargo test --release --test lumos_dia -- --ignored
//! ```
//!
//! OpenTFRaw 5e7d26d found that this 232/288 stream also fits the 232/344 size equation, which
//! put isolation windows on the wrong scans. We walk the event grammar instead of solving for
//! sizes, so this test pins the result rather than a fix: every event is checked against the
//! per-scan trailer, which the instrument writes independently. The trailer's conversion
//! coefficients differ from scan to scan, so a table shifted by even one event fails.

use thermorawfile::RawFile;

#[test]
#[ignore = "requires public PXD031322 OFL001513-YLL-GPF-15K-1.raw; set THERMORAWFILE_LUMOS_RAW"]
fn lumos_dia_events_align_with_scan_trailers() {
    let path = std::env::var("THERMORAWFILE_LUMOS_RAW").expect("set THERMORAWFILE_LUMOS_RAW");
    let rf = RawFile::open(&path).unwrap();
    assert!(rf.has_scan_events(), "scan-event walk failed on this layout");

    let (mut ms1, mut ms2) = (0u32, 0u32);
    let mut windows_in_cycle = 0u32;
    let mut cycle_lengths = std::collections::BTreeSet::new();
    let mut prev_center = 0.0f64;
    for scan in rf.first_scan..=rf.last_scan {
        let ev = rf.scan_event(scan).unwrap();
        let rec = rf.scan_params(scan).expect("scan trailer").record();

        let cal = rf
            .calibration_at_event(rf.scan_event_byte_offset(scan).unwrap())
            .unwrap_or_else(|| panic!("scan {scan}: no calibration in event"));
        assert_eq!(
            (Some(cal.a), Some(cal.b), Some(cal.c)),
            (
                rec.get_f64("Conversion Parameter A:"),
                rec.get_f64("Conversion Parameter B:"),
                rec.get_f64("Conversion Parameter C:"),
            ),
            "scan {scan}: event calibration disagrees with trailer"
        );

        if ev.ms_order == 1 {
            ms1 += 1;
            if windows_in_cycle > 0 {
                cycle_lengths.insert(windows_in_cycle);
            }
            windows_in_cycle = 0;
            prev_center = 0.0;
            continue;
        }
        ms2 += 1;
        windows_in_cycle += 1;
        assert!(
            (300.0..1100.0).contains(&ev.isolation_center) && ev.isolation_center > prev_center,
            "scan {scan}: isolation center {} out of DIA order (previous {prev_center})",
            ev.isolation_center
        );
        prev_center = ev.isolation_center;
        assert_eq!(Some(ev.isolation_width), rec.get_f64("MS2 Isolation Width:"), "scan {scan}");
        let hcd: Option<f64> = rec.get_string("HCD Energy:").and_then(|t| t.trim().parse().ok());
        assert_eq!(Some(ev.collision_energy), hcd, "scan {scan}");
    }
    cycle_lengths.insert(windows_in_cycle);

    // Counts and first centers as reported by OpenTFRaw for the same file.
    assert_eq!((ms1, ms2), (2308, 150020));
    assert_eq!(cycle_lengths.into_iter().collect::<Vec<_>>(), [65]);
    assert_eq!(
        rf.scan_filter(2).as_deref(),
        Some("FTMS + c NSI Full ms2 351.41@hcd32.00 [200.00-2000.00]")
    );
    let centers: Vec<f64> = (2..=4).map(|s| rf.scan_event(s).unwrap().isolation_center).collect();
    for (got, want) in centers.iter().zip([351.4096, 353.4105, 355.4114]) {
        assert!((got - want).abs() < 5e-5, "{centers:?}");
    }
}
