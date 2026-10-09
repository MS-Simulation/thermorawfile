//! Frequency↔m/z calibration round-trip (no data file needed).
//!
//! Uses the real coefficients decoded from a rev66 Orbitrap Astral MS1 scan
//! (PRIDE PXD074900): nparam=5, m/z = b/f² + c/f⁴.

use thermorawfile::Calibration;

#[test]
fn nparam5_inverse_round_trips() {
    let cal = Calibration {
        nparam: 5,
        a: 0.0,
        b: 1.6955368e8,
        c: 1.3376782e8,
    };
    for &mz in &[401.0, 500.0, 533.1506647, 650.0, 800.0, 899.0] {
        let f = cal.freq(mz).expect("reachable");
        let back = cal.mz(f);
        assert!(
            (back - mz).abs() < 1e-6,
            "m/z {mz} -> f {f} -> m/z {back} (Δ {})",
            (back - mz).abs()
        );
    }
}

#[test]
fn nparam4_inverse_round_trips() {
    let cal = Calibration {
        nparam: 4,
        a: -0.004,
        b: 3.34,
        c: 1.6955e8,
    };
    for &mz in &[300.0, 600.0, 1200.0] {
        let f = cal.freq(mz).expect("reachable");
        assert!((cal.mz(f) - mz).abs() < 1e-6);
    }
}

/// Both quadratic roots map back onto the target m/z, so `freq` must pick the physical one
/// (near the first-order `f = sqrt(b/mz)`), not the spurious root near `x = −b/c` (f ≈ 0.2).
/// Picking by round-trip error chose the spurious root for ~10% of m/z values on these
/// calibrations, which put the peak off the profile grid.
fn assert_physical_root(cal: Calibration, lo: f64, hi: f64) {
    let mut mz = lo;
    while mz <= hi {
        let f = cal.freq(mz).expect("reachable");
        let f0 = (cal.b / (mz - cal.a)).sqrt();
        assert!((f / f0 - 1.0).abs() < 1e-2, "m/z {mz}: f {f}, first order {f0}");
        assert!((cal.mz(f) - mz).abs() < 1e-6, "m/z {mz} -> f {f} -> {}", cal.mz(f));
        mz += 0.0137;
    }
}

#[test]
fn astral_picks_physical_root() {
    // Orbitrap Astral MS1 (LFQ_Astral_DIA_15min_50ng_Condition_A_REP1, scan 607), where
    // freq(744.4) used to return 0.1599 instead of 477.11.
    let cal = Calibration { nparam: 5, a: 0.0, b: 169452523.1179489, c: -4332928.760658183 };
    assert!((cal.freq(744.4).unwrap() - 477.1125).abs() < 1e-3);
    assert_physical_root(cal, 376.0, 990.0);
}

#[test]
fn fusion_lumos_picks_physical_root() {
    // Fusion Lumos MS1 (PRIDE PXD031322, scan 1).
    let cal = Calibration { nparam: 7, a: 0.0, b: 211782331.2992454, c: -270234696.3002101 };
    assert_physical_root(cal, 350.0, 1000.0);
}
