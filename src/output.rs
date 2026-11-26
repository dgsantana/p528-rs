//! The model's results and warnings.

use bitflags::bitflags;

/// Mode of propagation that sets the loss (NTIA's `propagation_mode`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PropagationMode {
    /// Line of sight (`PROP_MODE__LOS`, 1).
    LineOfSight,
    /// Diffraction (`PROP_MODE__DIFFRACTION`, 2).
    Diffraction,
    /// Troposcatter (`PROP_MODE__SCATTERING`, 3).
    Troposcatter,
}

impl PropagationMode {
    /// NTIA's numeric value for this mode.
    pub fn ntia_code(self) -> i32 {
        match self {
            Self::LineOfSight => 1,
            Self::Diffraction => 2,
            Self::Troposcatter => 3,
        }
    }
}

bitflags! {
    /// Conditions NTIA's reference flags on a result it still returns. The bit values are NTIA's.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
    pub struct Warnings: u32 {
        /// The diffraction and troposcatter models may not be physically consistent with each other on this path
        /// (`WARNING__DFRAC_TROPO_REGION`).
        const DIFFRACTION_TROPO_REGION = 0x01;
        /// The low terminal is above the 20 km the in-force Recommendation covers; the result is informative, not an
        /// official P.528 prediction (`WARNING__HEIGHT_LIMIT_H_1`).
        const HEIGHT_LIMIT_H1 = 0x02;
        /// The high terminal is above the 20 km the in-force Recommendation covers; the result is informative, not an
        /// official P.528 prediction (`WARNING__HEIGHT_LIMIT_H_2`).
        const HEIGHT_LIMIT_H2 = 0x04;
    }
}

/// One prediction's results (NTIA's `Result`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct P528Result {
    /// Basic transmission loss, dB.
    pub basic_transmission_loss_db: f64,
    /// Free-space basic transmission loss, dB.
    pub free_space_loss_db: f64,
    /// Median atmospheric (gaseous) absorption loss, dB.
    pub absorption_loss_db: f64,
    /// Path distance the loss is for, km. On a line-of-sight path this is the distance the ray geometry converged on,
    /// which can differ from the input by up to 1 m.
    pub distance_km: f64,
    /// Elevation angle of the ray at the low terminal, rad.
    pub theta_h1_rad: f64,
    /// Mode of propagation.
    pub mode: PropagationMode,
    /// Warnings.
    pub warnings: Warnings,
}
