//! Input validation (NTIA `ValidateInputs.cpp`).

use crate::{P528Error, P528Input, Warnings};

/// What validation lets through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Valid {
    /// The model can run.
    Run,
    /// The terminals are at the same height and 0 km apart (`ERROR_HEIGHT_AND_DISTANCE`). NTIA returns success with
    /// every loss 0.
    Coincident,
}

/// Checks `input` in NTIA's order, setting height warnings on the way as NTIA does (`ValidateInputs`). NaN, and an
/// infinite distance, are refused first; NTIA has no such check.
pub fn validate(input: &P528Input, warnings: &mut Warnings) -> Result<Valid, P528Error> {
    let P528Input {
        distance_km: d_km,
        h1_m,
        h2_m,
        frequency_mhz: f_mhz,
        time_percent: p,
        ..
    } = *input;

    if [d_km, h1_m, h2_m, f_mhz, p].iter().any(|v| v.is_nan()) || d_km.is_infinite() {
        return Err(P528Error::NotFinite);
    }

    if d_km < 0.0 {
        return Err(P528Error::Distance);
    }
    if !(1.5..=80000.0).contains(&h1_m) {
        return Err(P528Error::H1);
    }
    if !(1.5..=80000.0).contains(&h2_m) {
        return Err(P528Error::H2);
    }
    if h1_m > 20000.0 {
        *warnings |= Warnings::HEIGHT_LIMIT_H1;
    }
    if h2_m > 20000.0 {
        *warnings |= Warnings::HEIGHT_LIMIT_H2;
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
    if p < 1.0 {
        return Err(P528Error::PercentLow);
    }
    if p > 99.0 {
        return Err(P528Error::PercentHigh);
    }
    if h1_m == h2_m && d_km == 0.0 {
        return Ok(Valid::Coincident);
    }
    Ok(Valid::Run)
}
