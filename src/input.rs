//! The model's inputs.

/// Polarisation of the antennas (NTIA's `T_pol`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Polarization {
    /// Horizontal (`POLARIZATION__HORIZONTAL`, 0).
    Horizontal,
    /// Vertical (`POLARIZATION__VERTICAL`, 1).
    Vertical,
}

/// One prediction's inputs. The limits are those NTIA's reference enforces; outside them [`crate::p528`] returns a
/// [`crate::P528Error`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct P528Input {
    /// Great-circle path distance between the terminals, km. At least 0.
    pub distance_km: f64,
    /// Height of the low terminal, m. 1.5 to 80,000 m, and not above `h2_m`; above 20,000 m the result carries
    /// [`crate::Warnings::HEIGHT_LIMIT_H1`].
    pub h1_m: f64,
    /// Height of the high terminal, m. 1.5 to 80,000 m; above 20,000 m the result carries
    /// [`crate::Warnings::HEIGHT_LIMIT_H2`].
    pub h2_m: f64,
    /// Frequency, MHz. 100 to 30,000 MHz.
    pub frequency_mhz: f64,
    /// Time percentage. 1 to 99.
    pub time_percent: f64,
    /// Polarisation.
    pub polarization: Polarization,
}
