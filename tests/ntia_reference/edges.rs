//! Edge cases: refusals, warnings, the line-of-sight boundary and the limits of the ranges.

use p528::{P528Error, P528Model, Polarization, p528};

use crate::compare::{assert_all_agree, check, input};

/// A NaN or an infinity that NTIA's range checks catch gets NTIA's code: NaN passes NTIA's comparisons, so the first
/// range check that fails decides, as in NTIA's.
#[test]
fn non_finite_inputs_ntia_refuses_match() {
    let cases = [
        (f64::NEG_INFINITY, 10.0, 1000.0, 1000.0, 0, 50.0), // code 1
        (f64::NAN, 1.0, 1000.0, 1000.0, 0, 50.0),           // code 2
        (10.0, 10.0, f64::INFINITY, 1000.0, 0, 50.0),       // code 3
        (10.0, 1000.0, 10.0, f64::NAN, 0, 50.0),            // code 4
        (10.0, 10.0, 1000.0, f64::INFINITY, 0, 50.0),       // code 6
        (10.0, f64::NAN, 1000.0, 1000.0, 0, 0.5),           // code 7
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&(d, h1, h2, f, pol, p)| check(d, h1, h2, f, pol, p))
        .collect();
    assert_all_agree(&failures, cases.len());
}

/// A NaN that passes every range check, or an infinite distance, is refused with `NotFinite` by the single call and
/// by the model, where NTIA's reference computes on and returns NaN.
#[test]
fn non_finite_inputs_ntia_accepts_are_refused() {
    let nan = f64::NAN;
    let cases = [
        (nan, 10.0, 1000.0, 1000.0, 50.0),
        (f64::INFINITY, 10.0, 1000.0, 1000.0, 50.0),
        (10.0, nan, 1000.0, 1000.0, 50.0),
        (10.0, 10.0, nan, 1000.0, 50.0),
        (10.0, 10.0, 1000.0, nan, 50.0),
        (10.0, 10.0, 1000.0, 1000.0, nan),
    ];
    for (d, h1, h2, f, p) in cases {
        let single = p528(&input(d, h1, h2, f, 0, p));
        assert_eq!(
            single,
            Err(P528Error::NotFinite),
            "p528 d={d} h1={h1} h2={h2} f={f} p={p}"
        );
        let model =
            P528Model::create(h1, h2, f, Polarization::Horizontal).and_then(|m| m.loss(d, p));
        assert_eq!(
            model,
            Err(P528Error::NotFinite),
            "model d={d} h1={h1} h2={h2} f={f} p={p}"
        );
    }
}

/// Every input the reference refuses is refused by the port with the same code, including where several checks fail
/// at once (the first in NTIA's order wins).
#[test]
fn refusals_match() {
    let cases = [
        (-1.0, 10.0, 1000.0, 1000.0, 0, 50.0),
        (10.0, 1.4, 1000.0, 1000.0, 0, 50.0),
        (10.0, 10.0, 1.4, 1000.0, 0, 50.0),
        (10.0, 10.0, 80001.0, 1000.0, 0, 50.0),
        (10.0, 80001.0, 80002.0, 1000.0, 0, 50.0),
        (10.0, 1000.0, 10.0, 1000.0, 0, 50.0), // terminals swapped: refused, not reordered
        (10.0, 10.0, 1000.0, 99.9, 0, 50.0),
        (10.0, 10.0, 1000.0, 30000.1, 0, 50.0),
        (10.0, 10.0, 1000.0, 1000.0, 0, 0.9),
        (10.0, 10.0, 1000.0, 1000.0, 0, 99.1),
        (-1.0, 1.0, 0.5, 50.0, 0, 0.0),
        (10.0, 25000.0, 21000.0, 50.0, 0, 50.0),
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&(d, h1, h2, f, pol, p)| check(d, h1, h2, f, pol, p))
        .collect();
    assert_all_agree(&failures, cases.len());
}

/// Terminals at the same height and 0 km apart: NTIA returns success with every loss 0 and no mode; the port refuses
/// with NTIA's code for the case, 10.
#[test]
fn coincident_terminals_match() {
    let cases = [
        (0.0, 100.0, 100.0),
        (0.0, 1.5, 1.5),
        (0.0, 30000.0, 30000.0),
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&(d, h1, h2)| check(d, h1, h2, 1000.0, 0, 50.0))
        .collect();
    assert_all_agree(&failures, cases.len());
}

/// Heights above 20 km, which the reference accepts up to 80 km with a warning per terminal.
#[test]
fn height_warnings_match() {
    let mut failures = Vec::new();
    let mut n = 0;
    for (h1, h2) in [
        (10.0, 20000.0),
        (10.0, 20000.1),
        (10.0, 80000.0),
        (25000.0, 80000.0),
    ] {
        for d in [1.0, 100.0, 600.0, 1500.0] {
            failures.extend(check(d, h1, h2, 1000.0, 1, 50.0));
            n += 1;
        }
    }
    assert_all_agree(&failures, n);
}

/// Distances either side of the line-of-sight limit d_ML (the model switches 1 m inside it), at zero, and far beyond.
#[test]
fn line_of_sight_boundary_matches() {
    let mut failures = Vec::new();
    let mut n = 0;
    for (h1, h2, f) in [
        (1.5, 1.5, 100.0),
        (10.0, 3000.0, 1000.0),
        (100.0, 20000.0, 30000.0),
        (1.5, 20000.0, 5000.0),
    ] {
        // d_ML from the reference itself, then distances around it.
        let d_ml = ntia_p528_reference::p528(1.0, h1, h2, f, 0, 50.0)
            .path
            .d_ml_km;
        for offset in [
            -1.0, -0.01, -0.0011, -0.001, -0.0009, 0.0, 0.001, 0.01, 1.0, 3.0,
        ] {
            failures.extend(check(d_ml + offset, h1, h2, f, 0, 50.0));
            n += 1;
        }
        for d in [0.0, 1e-6, 1e-3, 5000.0, 10000.0] {
            failures.extend(check(d, h1, h2, f, 1, 10.0));
            n += 1;
        }
    }
    assert_all_agree(&failures, n);
}

/// The corners of the frequency, height and time ranges.
#[test]
fn range_corners_match() {
    let mut failures = Vec::new();
    let mut n = 0;
    for f in [100.0, 30000.0] {
        for (h1, h2) in [(1.5, 1.5), (1.5, 20000.0), (20000.0, 20000.0)] {
            for p in [1.0, 50.0, 99.0] {
                for pol in [0, 1] {
                    for d in [0.5, 50.0, 500.0] {
                        failures.extend(check(d, h1, h2, f, pol, p));
                        n += 1;
                    }
                }
            }
        }
    }
    assert_all_agree(&failures, n);
}
