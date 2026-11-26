//! NTIA's published examples (the reference repository's README), against the port and against the reference.

use crate::compare::{assert_all_agree, check, input};

/// `(d km, h1 m, h2 m, f MHz, polarisation, p %, A dB)`, as NTIA publishes them, rounded to 0.1 dB.
const EXAMPLES: [(f64, f64, f64, f64, i32, f64, f64); 4] = [
    (15.0, 10.0, 1000.0, 500.0, 0, 50.0, 110.0),
    (100.0, 100.0, 15000.0, 3600.0, 0, 90.0, 151.6),
    (1500.0, 15.0, 10000.0, 5700.0, 0, 10.0, 293.4),
    (30.0, 8.0, 20000.0, 22000.0, 1, 50.0, 151.1),
];

#[test]
fn published_examples() {
    let mut failures = Vec::new();
    for (d, h1, h2, f, pol, p, expected) in EXAMPLES {
        let port = p528::p528(&input(d, h1, h2, f, pol, p)).expect("an accepted input");
        if (port.basic_transmission_loss_db - expected).abs() > 0.05 {
            failures.push(format!(
                "{d} km: port {} dB, published {expected} dB",
                port.basic_transmission_loss_db
            ));
        }
        failures.extend(check(d, h1, h2, f, pol, p));
    }
    assert_all_agree(&failures, EXAMPLES.len());
}
