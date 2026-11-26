//! Inputs the model refuses.

use std::fmt;

/// Why the model refused an input or could not compute a result.
///
/// The variants with an NTIA return code are checked in NTIA's order, so an input NTIA refuses gives the same code.
/// The others cover what NTIA's reference lets through and answers with a value that is not a loss: NaN, or 0 dB for
/// coincident terminals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum P528Error {
    /// An input is NaN, or the distance is infinite. NTIA's reference has no such check: a NaN passes its
    /// validation and an infinite distance its range checks, and the computation then returns NaN.
    NotFinite,
    /// The distance is below 0 km (`ERROR_VALIDATION__D_KM`, 1).
    Distance,
    /// The low terminal is below 1.5 m or above 80,000 m (`ERROR_VALIDATION__H_1`, 2).
    H1,
    /// The high terminal is below 1.5 m or above 80,000 m (`ERROR_VALIDATION__H_2`, 3).
    H2,
    /// The low terminal is above the high terminal (`ERROR_VALIDATION__TERM_GEO`, 4). The model does not reorder
    /// them; the caller passes the lower one as `h1_m`.
    TerminalGeometry,
    /// The frequency is below 100 MHz (`ERROR_VALIDATION__F_MHZ_LOW`, 5).
    FrequencyLow,
    /// The frequency is above 30,000 MHz (`ERROR_VALIDATION__F_MHZ_HIGH`, 6).
    FrequencyHigh,
    /// The time percentage is below 1 (`ERROR_VALIDATION__PERCENT_LOW`, 7).
    PercentLow,
    /// The time percentage is above 99 (`ERROR_VALIDATION__PERCENT_HIGH`, 8).
    PercentHigh,
    /// The terminals are at the same height and 0 km apart (`ERROR_HEIGHT_AND_DISTANCE`, 10). NTIA defines this code
    /// but returns success with every loss 0.
    CoincidentTerminals,
    /// The path leaves the reference atmosphere (P.835): on long transhorizon paths, from about 2,400 km, the
    /// troposcatter common volume is above 100 km, the top of the atmosphere's profiles. NTIA's reference returns NaN
    /// with a success code.
    OutsideAtmosphere,
}

impl P528Error {
    /// NTIA's return code for this error; `None` for the errors NTIA does not have.
    pub fn ntia_code(self) -> Option<i32> {
        match self {
            Self::NotFinite => None,
            Self::Distance => Some(1),
            Self::H1 => Some(2),
            Self::H2 => Some(3),
            Self::TerminalGeometry => Some(4),
            Self::FrequencyLow => Some(5),
            Self::FrequencyHigh => Some(6),
            Self::PercentLow => Some(7),
            Self::PercentHigh => Some(8),
            Self::CoincidentTerminals => Some(10),
            Self::OutsideAtmosphere => None,
        }
    }
}

impl fmt::Display for P528Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::NotFinite => "an input is NaN, or the distance is infinite",
            Self::Distance => "path distance must be at least 0 km",
            Self::H1 => "low terminal height must be 1.5 m to 80,000 m",
            Self::H2 => "high terminal height must be 1.5 m to 80,000 m",
            Self::TerminalGeometry => "low terminal height must not exceed high terminal height",
            Self::FrequencyLow => "frequency must be at least 100 MHz",
            Self::FrequencyHigh => "frequency must not exceed 30,000 MHz",
            Self::PercentLow => "time percentage must be at least 1",
            Self::PercentHigh => "time percentage must not exceed 99",
            Self::CoincidentTerminals => "the terminals are at the same height and 0 km apart",
            Self::OutsideAtmosphere => "the path leaves the reference atmosphere",
        })
    }
}

impl std::error::Error for P528Error {}
