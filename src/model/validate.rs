//! Input validation (NTIA `ValidateInputs.cpp`), split into the checks on the terminals and frequency and the checks
//! on the distance and time percentage, so a model built for one pair of heights checks each part once.

use crate::{P528Error, Warnings};

/// Checks the terminal heights and the frequency in NTIA's order, returning the height warnings.
pub fn terminals(h1_m: f64, h2_m: f64, f_mhz: f64) -> Result<Warnings, P528Error> {
    if h1_m.is_nan() || h2_m.is_nan() || f_mhz.is_nan() {
        return Err(P528Error::NotFinite);
    }
    if !(1.5..=80000.0).contains(&h1_m) {
        return Err(P528Error::H1);
    }
    if !(1.5..=80000.0).contains(&h2_m) {
        return Err(P528Error::H2);
    }
    let mut warnings = Warnings::empty();
    if h1_m > 20000.0 {
        warnings |= Warnings::HEIGHT_LIMIT_H1;
    }
    if h2_m > 20000.0 {
        warnings |= Warnings::HEIGHT_LIMIT_H2;
    }
    if h1_m > h2_m {
        return Err(P528Error::TerminalGeometry);
    }
    if f_mhz < 100.0 {
        return Err(P528Error::FrequencyLow);
    }
    if f_mhz > 30000.0 {
        return Err(P528Error::FrequencyHigh);
    }
    // NTIA checks the polarisation here; `Polarization` cannot hold an invalid value.
    Ok(warnings)
}

/// Checks the distance, in NTIA's order its first check.
pub fn distance(d_km: f64) -> Result<(), P528Error> {
    if d_km.is_nan() || d_km.is_infinite() {
        return Err(P528Error::NotFinite);
    }
    if d_km < 0.0 {
        return Err(P528Error::Distance);
    }
    Ok(())
}

/// Checks the time percentage and, last, coincident terminals, in NTIA's order.
pub fn percentage(p: f64, coincident: bool) -> Result<(), P528Error> {
    if p.is_nan() {
        return Err(P528Error::NotFinite);
    }
    if p < 1.0 {
        return Err(P528Error::PercentLow);
    }
    if p > 99.0 {
        return Err(P528Error::PercentHigh);
    }
    if coincident {
        return Err(P528Error::CoincidentTerminals);
    }
    Ok(())
}
