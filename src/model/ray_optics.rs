//! Ray optics over the adjusted earth (NTIA `RayOptics.cpp`, P.528 Annex 2, section 7).

use super::terminal::Terminal;
use super::{A_0_KM, A_E_KM, max};

/// Line-of-sight geometry for one grazing angle, and the line-of-sight loss once computed (`LineOfSightParams`).
#[derive(Debug, Default, Clone, Copy)]
pub struct LineOfSightParams {
    /// Terminal radii from the adjusted earth's centre, km.
    pub z_km: [f64; 2],
    /// Path distance between the terminals, km.
    pub d_km: f64,
    /// Direct ray length, km.
    pub r_0_km: f64,
    /// Indirect (reflected) ray length, km.
    pub r_12_km: f64,
    /// Distances from each terminal to the reflection point along the tangent plane, km.
    pub big_d_km: [f64; 2],
    /// Take-off angle from the low terminal to the high terminal, rad.
    pub theta_h1_rad: f64,
    /// Take-off angle from the high terminal to the low terminal, rad.
    pub theta_h2_rad: f64,
    /// Central angles from each terminal to the reflection point, rad.
    pub theta: [f64; 2],
    /// Adjusted earth radius, km.
    pub a_a_km: f64,
    /// Path length difference between the reflected and direct rays, km.
    pub delta_r_km: f64,
    /// Line-of-sight loss, dB. Set by `path_loss`, not by `ray_optics`.
    pub a_los_db: f64,
}

/// Geometry of the direct and reflected rays for grazing angle `psi` (`RayOptics`). Leaves `a_los_db` as it is in
/// `params`, as NTIA's function leaves that field untouched.
pub fn ray_optics(
    terminal_1: &Terminal,
    terminal_2: &Terminal,
    psi: f64,
    params: &mut LineOfSightParams,
) {
    let z = (A_0_KM / A_E_KM) - 1.0; // [Eqn 7-1]
    let k_a = 1.0 / (1.0 + z * psi.cos()); // [Eqn 7-2]
    params.a_a_km = A_0_KM * k_a; // [Eqn 7-3]

    let delta_h_a1_km = terminal_1.delta_h_km * (params.a_a_km - A_0_KM) / (A_E_KM - A_0_KM); // [Eqn 7-4]
    let delta_h_a2_km = terminal_2.delta_h_km * (params.a_a_km - A_0_KM) / (A_E_KM - A_0_KM); // [Eqn 7-4]

    let h_km = [
        terminal_1.h_r_km - delta_h_a1_km, // [Eqn 7-5]
        terminal_2.h_r_km - delta_h_a2_km, // [Eqn 7-5]
    ];

    let mut h_prime_km = [0.0; 2];
    for i in 0..2 {
        params.z_km[i] = params.a_a_km + h_km[i]; // [Eqn 7-6]
        params.theta[i] = (params.a_a_km * psi.cos() / params.z_km[i]).acos() - psi; // [Eqn 7-7]
        params.big_d_km[i] = params.z_km[i] * params.theta[i].sin(); // [Eqn 7-8]

        // [Eqn 7-9]
        h_prime_km[i] = if psi > 1.56 {
            h_km[i]
        } else {
            params.big_d_km[i] * psi.tan()
        };
    }

    let delta_z = (params.z_km[0] - params.z_km[1]).abs(); // [Eqn 7-10]

    params.d_km = max(params.a_a_km * (params.theta[0] + params.theta[1]), 0.0); // [Eqn 7-11]

    let d_sum = params.big_d_km[0] + params.big_d_km[1];
    let alpha = ((h_prime_km[1] - h_prime_km[0]) / d_sum).atan(); // [Eqn 7-12]
    params.r_0_km = max(delta_z, d_sum / alpha.cos()); // [Eqn 7-13]
    params.r_12_km = d_sum / psi.cos(); // [Eqn 7-14]

    params.delta_r_km = 4.0 * h_prime_km[0] * h_prime_km[1] / (params.r_0_km + params.r_12_km); // [Eqn 7-15]

    params.theta_h1_rad = alpha - params.theta[0]; // [Eqn 7-16]
    params.theta_h2_rad = -(alpha + params.theta[1]); // [Eqn 7-17]
}
