//! Terminal geometry (NTIA `TerminalGeometry.cpp`, P.528 Annex 2, Step 1 and section 4).

use std::f64::consts::PI;

use super::{A_0_KM, A_E_KM};
use crate::p676::slant_path_attenuation;
use crate::p835::OutsideAtmosphere;

/// One terminal's geometry (`Terminal`).
#[derive(Debug, Default, Clone, Copy)]
pub struct Terminal {
    /// Real terminal height, km.
    pub h_r_km: f64,
    /// Effective terminal height, km.
    pub h_e_km: f64,
    /// Difference between real and effective height, km.
    pub delta_h_km: f64,
    /// Ray-traced horizon distance, km.
    pub d_r_km: f64,
    /// Total ray path length to the horizon, km.
    pub a_km: f64,
    /// Incident angle of the grazing ray at the terminal, rad.
    pub theta_rad: f64,
    /// Median atmospheric absorption loss to the horizon, dB.
    pub a_a_db: f64,
}

impl Terminal {
    /// Geometry of a terminal at `h_r_km` for frequency `f_mhz`, from a ray traced from the ground at grazing
    /// incidence up to the terminal (`TerminalGeometry`). Fails only for a height outside the reference atmosphere.
    pub fn create(h_r_km: f64, f_mhz: f64) -> Result<Self, OutsideAtmosphere> {
        let theta_tx_rad = 0.0;
        let slant = slant_path_attenuation(f_mhz / 1000.0, 0.0, h_r_km, PI / 2.0 - theta_tx_rad)?;

        // Arc distance: [Thayer, equation 2], rearranged.
        let central_angle = (PI / 2.0 - slant.angle_rad) - theta_tx_rad + slant.bending_rad;
        let d_r_km = A_0_KM * central_angle;

        let phi_rad = d_r_km / A_E_KM; // [Eqn 4-1]
        let h_e_km = (A_E_KM / phi_rad.cos()) - A_E_KM; // [Eqn 4-2]

        Ok(Self {
            h_r_km,
            h_e_km,
            delta_h_km: h_r_km - h_e_km, // [Eqn 4-3]
            d_r_km,
            a_km: slant.a_km,
            theta_rad: PI / 2.0 - slant.angle_rad,
            a_a_db: slant.a_gas_db,
        })
    }
}
