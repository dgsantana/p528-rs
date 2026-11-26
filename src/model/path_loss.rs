//! Line-of-sight path loss from the two-ray model and the ground reflection coefficients (NTIA `GetPathLoss.cpp`,
//! `ReflectionCoefficients.cpp`, P.528 Annex 2, sections 8 and 9).

use std::f64::consts::PI;

use super::ray_optics::LineOfSightParams;
use super::{EPSILON_R, SIGMA, min};
use crate::Polarization;

/// Path distances shared by the line-of-sight steps (`Path`).
#[derive(Debug, Default, Clone, Copy)]
pub struct Path {
    /// Maximum line-of-sight distance, km.
    pub d_ml_km: f64,
    /// `d_0`, the distance up to which the two-ray model is used, km.
    pub d_0_km: f64,
    /// Distance at which smooth-earth diffraction loss is 0 dB, km.
    pub d_d_km: f64,
}

/// Magnitude and phase of the ground reflection coefficient at grazing angle `psi_rad` (`ReflectionCoefficients`).
fn reflection_coefficients(psi_rad: f64, f_mhz: f64, polarization: Polarization) -> (f64, f64) {
    let (sin_psi, cos_psi) = if psi_rad <= 0.0 {
        (0.0, 1.0)
    } else if psi_rad >= PI / 2.0 {
        (1.0, 0.0)
    } else {
        (psi_rad.sin(), psi_rad.cos())
    };

    let x = (18000.0 * SIGMA) / f_mhz; // [Eqn 9-1]
    let y = EPSILON_R - cos_psi.powf(2.0); // [Eqn 9-2]
    let t = (y.powf(2.0) + x.powf(2.0)).sqrt() + y; // [Eqn 9-3]
    let p = (t * 0.5).sqrt(); // [Eqn 9-4]
    let q = x / (2.0 * p); // [Eqn 9-5]

    let pq = p.powf(2.0) + q.powf(2.0);
    let (b, a, alpha, beta) = match polarization {
        Polarization::Horizontal => (
            1.0 / pq,                // [Eqn 9-6]
            (2.0 * p) / pq,          // [Eqn 9-7]
            (-q).atan2(sin_psi - p), // [Eqn 9-9]
            q.atan2(sin_psi + p),    // [Eqn 9-10]
        ),
        Polarization::Vertical => (
            (EPSILON_R.powf(2.0) + x.powf(2.0)) / pq, // [Eqn 9-6]
            (2.0 * (p * EPSILON_R + q * x)) / pq,     // [Eqn 9-7]
            ((EPSILON_R * sin_psi) - q).atan2(EPSILON_R * sin_psi - p), // [Eqn 9-9]
            ((x * sin_psi) + q).atan2(EPSILON_R * sin_psi + p), // [Eqn 9-10]
        ),
    };

    // [Eqn 9-8]
    let r_g = ((1.0 + (b * sin_psi.powf(2.0)) - (a * sin_psi))
        / (1.0 + (b * sin_psi.powf(2.0)) + (a * sin_psi)))
        .sqrt();

    (r_g, alpha - beta) // [Eqn 9-11]
}

/// Line-of-sight loss at grazing angle `psi_rad`, written to `params.a_los_db`; returns `R_Tg`, the effective
/// reflection coefficient (`GetPathLoss`). Beyond `d_0` the loss is the straight line from `a_d_0_db` at `d_0` to
/// `a_dml_db` at `d_ML`.
#[allow(clippy::too_many_arguments)]
pub fn path_loss(
    psi_rad: f64,
    path: &Path,
    f_mhz: f64,
    psi_limit: f64,
    a_dml_db: f64,
    a_d_0_db: f64,
    polarization: Polarization,
    params: &mut LineOfSightParams,
) -> f64 {
    let (r_g, phi_g) = reflection_coefficients(psi_rad, f_mhz, polarization);

    let d_v = if psi_rad.tan() >= 0.1 {
        1.0
    } else {
        let r_1 = params.big_d_km[0] / psi_rad.cos(); // [Eqn 8-3]
        let r_2 = params.big_d_km[1] / psi_rad.cos(); // [Eqn 8-3]
        let r_r = (r_1 * r_2) / params.r_12_km; // [Eqn 8-4]

        let term_1 =
            (2.0 * r_r * (1.0 + psi_rad.sin().powf(2.0))) / (params.a_a_km * psi_rad.sin());
        let term_2 = (2.0 * r_r / params.a_a_km).powf(2.0);
        (1.0 + term_1 + term_2).powf(-0.5) // [Eqn 8-5]
    };

    let f_r = min(params.r_0_km / params.r_12_km, 1.0); // Ray-length factor, [Eqn 8-6]

    let r_tg = r_g * d_v * f_r; // [Eqn 8-7]

    if params.d_km > path.d_0_km {
        // [Eqn 8-1]
        params.a_los_db = ((params.d_km - path.d_0_km) * (a_dml_db - a_d_0_db)
            / (path.d_ml_km - path.d_0_km))
            + a_d_0_db;
    } else {
        let lambda_km = 0.2997925 / f_mhz; // [Eqn 8-2]

        if psi_rad > psi_limit {
            // Ignore the phase lag; Step 8-2.
            params.a_los_db = 0.0;
        } else {
            // Total phase lag of the ground-reflected ray relative to the direct ray, [Eqn 8-8].
            let phi_tg = (2.0 * PI * params.delta_r_km / lambda_km) + phi_g;

            // |1 + R_Tg e^(-j phi_Tg)|, [Eqns 8-9 and 8-10].
            let w_rl = min((1.0 + r_tg * phi_tg.cos()).hypot(-r_tg * phi_tg.sin()), 1.0);
            let w_r0 = w_rl.powf(2.0); // [Eqn 8-11]

            params.a_los_db = 10.0 * w_r0.log10(); // [Eqn 8-12]
        }
    }

    r_tg
}
