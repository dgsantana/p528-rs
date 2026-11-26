//! Generated cases across the ranges the reference accepts.

use crate::compare::{assert_all_agree, check};

/// A small deterministic generator (xorshift64*), so every run checks the same cases without a dependency.
struct Cases(u64);

impl Cases {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.unit()
    }
    /// Log-uniform: as many cases per decade, so short paths, low heights and VHF are as present as the rest.
    fn log_range(&mut self, low: f64, high: f64) -> f64 {
        self.range(low.ln(), high.ln()).exp()
    }
    /// One case: distance, ordered heights, frequency, polarisation; the time percentage is the caller's.
    fn geometry(&mut self) -> (f64, f64, f64, f64, i32) {
        let d = self.log_range(0.01, 2000.0);
        let a = self.log_range(1.5, 20000.0);
        let b = self.log_range(1.5, 20000.0);
        let f = self.log_range(100.0, 30000.0);
        let pol = i32::from(self.unit() < 0.5);
        (d, a.min(b), a.max(b), f, pol)
    }
}

/// Generated cases; the `P528_CASES` environment variable overrides the count.
fn count() -> usize {
    std::env::var("P528_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2000)
}

#[test]
fn generated_cases_agree() {
    let mut cases = Cases(0x5052_3532_385F_5253);
    let n = count();
    let mut failures = Vec::new();
    for _ in 0..n {
        let (d, h1, h2, f, pol) = cases.geometry();
        let p = cases.range(1.0, 99.0);
        failures.extend(check(d, h1, h2, f, pol, p));
    }
    assert_all_agree(&failures, n);
}

/// The tabulated time percentages, where the interpolation in the Nakagami-Rice and variability tables changes
/// segment, at a spread of geometries.
#[test]
fn tabulated_percentages_agree() {
    let percentages = [
        1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 85.0, 90.0, 95.0,
        98.0, 99.0,
    ];
    let mut cases = Cases(0x7461_626C_6573);
    let mut failures = Vec::new();
    let mut n = 0;
    for p in percentages {
        for _ in 0..20 {
            let (d, h1, h2, f, pol) = cases.geometry();
            failures.extend(check(d, h1, h2, f, pol, p));
            n += 1;
        }
    }
    assert_all_agree(&failures, n);
}
