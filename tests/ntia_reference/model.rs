//! One `P528Model` reused across many distances and time percentages, as a batch caller uses it, against NTIA's
//! reference called afresh for each case. A result that depended on an earlier call would show here.

use p528::{P528Model, Polarization};

use crate::compare::{assert_all_agree, compare};

/// Height pairs and frequencies for a radar near the ground and targets up to 20 km, plus two high-low pairs.
const MODELS: [(f64, f64, f64, i32); 6] = [
    (10.0, 100.0, 1000.0, 0),
    (5.0, 3000.0, 3000.0, 1),
    (30.0, 5000.0, 1300.0, 0),
    (15.0, 20000.0, 9400.0, 1),
    (1000.0, 1000.0, 300.0, 0),
    (8000.0, 12000.0, 30000.0, 1),
];

/// Distances out of order and repeated, from inside the line-of-sight region to well beyond the horizon.
const DISTANCES: [f64; 12] = [
    80.0, 0.5, 400.0, 30.0, 0.5, 2500.0, 150.0, 5.0, 900.0, 80.0, 0.0, 3000.0,
];

#[test]
fn reused_model_agrees_with_reference() {
    let mut failures = Vec::new();
    let mut n = 0;
    for (h1, h2, f, pol) in MODELS {
        let polarization = if pol == 0 {
            Polarization::Horizontal
        } else {
            Polarization::Vertical
        };
        let model = P528Model::create(h1, h2, f, polarization).expect("an accepted model");
        for (k, d) in DISTANCES.into_iter().enumerate() {
            let p = [1.0, 10.0, 50.0, 90.0, 99.0, 37.5][k % 6];
            let label = format!("model h1={h1} h2={h2} f={f} pol={pol}: d={d} p={p}");
            let reference = ntia_p528_reference::p528(d, h1, h2, f, pol, p);
            failures.extend(compare(&label, model.loss(d, p), &reference));
            n += 1;
        }
    }
    assert_all_agree(&failures, n);
}

/// The parallel batch returns, in order, what the sequential calls return.
#[cfg(feature = "parallel")]
#[test]
fn parallel_batch_matches_sequential() {
    let model =
        P528Model::create(10.0, 5000.0, 3000.0, Polarization::Vertical).expect("an accepted model");
    let distances: Vec<f64> = (0..200).map(|i| f64::from(i) * 2.5).collect();
    let sequential: Vec<_> = distances.iter().map(|&d| model.loss(d, 50.0)).collect();
    assert_eq!(model.par_losses(&distances, 50.0), sequential);
}
