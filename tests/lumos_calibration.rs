//! Calibration on a Fusion Lumos file, whose MS1 scan events have a 96-byte body.
//!
//! Original RAW from the public PRIDE study PXD031322 (Fusion Lumos). It is not
//! redistributed here; download it and point `THERMORAWFILE_LUMOS_RAW` at it:
//!
//! ```text
//! curl -O https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/OFL001513-YLL-GPF-15K-1.raw
//! THERMORAWFILE_LUMOS_RAW=$PWD/OFL001513-YLL-GPF-15K-1.raw cargo test --release --test lumos_calibration -- --ignored
//! ```
//!
//! The scan event stores the frequency→m/z calibration; the per-scan trailer stores the
//! same coefficients as "Conversion Parameter A/B/C". The two are written independently
//! by the instrument, so agreement on every MS1 scan, with the record found at event+160,
//! checks that offset for this layout (OpenTFRaw 6ec90a3).

use thermorawfile::RawFile;

#[test]
#[ignore = "requires public PXD031322 OFL001513-YLL-GPF-15K-1.raw; set THERMORAWFILE_LUMOS_RAW"]
fn lumos_ms1_calibration_matches_scan_trailers() {
    let path = std::env::var("THERMORAWFILE_LUMOS_RAW").expect("set THERMORAWFILE_LUMOS_RAW");
    let rf = RawFile::open(&path).unwrap();
    assert!(rf.has_scan_events(), "scan-event walk failed on this layout");

    let mut ms1 = 0;
    let mut mismatches = Vec::new();
    for scan in rf.first_scan..=rf.last_scan {
        if rf.scan_event(scan).unwrap().ms_order != 1 {
            continue;
        }
        ms1 += 1;
        let off = rf.scan_event_byte_offset(scan).unwrap();
        // The record must sit at the 96-byte-body offset itself, not be found by the
        // fallback scan: nparam 7 at event+160, and the +216 slot not a record at all.
        let u32at = |o: usize| u32::from_le_bytes(rf.bytes[o..o + 4].try_into().unwrap());
        assert_eq!(u32at(off + 160), 7, "scan {scan}: no nparam at event+160");
        assert!(!matches!(u32at(off + 216), 4 | 5 | 7), "scan {scan}: +216 also looks like a record");
        let cal = rf.calibration_at_event(off);
        let params = rf.scan_params(scan).expect("scan trailer");
        let rec = params.record();
        let want = (
            rec.get_f64("Conversion Parameter A:"),
            rec.get_f64("Conversion Parameter B:"),
            rec.get_f64("Conversion Parameter C:"),
        );
        let got = cal.map(|c| (Some(c.a), Some(c.b), Some(c.c)));
        if got != Some(want) {
            mismatches.push((scan, got, want));
        }
    }
    assert!(ms1 > 0, "no MS1 scans found");
    assert!(
        mismatches.is_empty(),
        "{} of {ms1} MS1 scans disagree with their trailer; first: {:?}",
        mismatches.len(),
        &mismatches[..mismatches.len().min(3)]
    );
}
