//! Input validation (NTIA `ValidateInputs.cpp`), split into the checks on the distance, on the terminals and frequency,
//! and on the time percentage, so a model built for one pair of heights checks each part once.
//!
//! The range checks are NTIA's comparisons, so a NaN passes them as it passes NTIA's, and an input NTIA refuses gets
//! NTIA's code. The checks for NaN and infinity, which NTIA does not have, come after the range checks.

use crate::{P528Error, Warnings};

/// NTIA's distance check, its first.
pub fn distance(d_km: f64) -> Result<(), P528Error> {
    if d_km < 0.0 {
        return Err(P528Error::Distance);
    }
    Ok(())
}

/// NTIA's checks on the terminal heights and the frequency, in NTIA's order, returning the height warnings.
pub fn terminals(h1_m: f64, h2_m: f64, f_mhz: f64) -> Result<Warnings, P528Error> {
    if h1_m < 1.5 || h1_m > 80000.0 {
        return Err(P528Error::H1);
    }
    if h2_m < 1.5 || h2_m > 80000.0 {
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

/// NTIA's checks on the time percentage.
pub fn percentage(p: f64) -> Result<(), P528Error> {
    if p < 1.0 {
        return Err(P528Error::PercentLow);
    }
    if p > 99.0 {
        return Err(P528Error::PercentHigh);
    }
    Ok(())
}

/// The port's own check, after NTIA's range checks: no input is NaN, and the distance is finite. (The other inputs
/// are finite once NTIA's range checks pass.)
pub fn finite(values: &[f64], d_km: Option<f64>) -> Result<(), P528Error> {
    if values.iter().any(|v| v.is_nan()) || d_km.is_some_and(|d| !d.is_finite()) {
        return Err(P528Error::NotFinite);
    }
    Ok(())
}

/// NTIA's last check: the terminals are at the same height and 0 km apart.
pub fn coincident(h1_m: f64, h2_m: f64, d_km: f64) -> Result<(), P528Error> {
    if h1_m == h2_m && d_km == 0.0 {
        return Err(P528Error::CoincidentTerminals);
    }
    Ok(())
}
