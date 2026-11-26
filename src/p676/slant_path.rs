//! Slant-path attenuation by atmospheric gases (NTIA `SlantPath.cpp`).

use std::f64::consts::PI;

use super::ray_trace::ray_trace;
use super::specific_attenuation::{global_wet_pressure, refractive_index};
use crate::model::A_0_KM;
use crate::p835::{global_pressure, global_temperature};

/// The result of a slant-path ray trace (`SlantPathAttenuationResult`).
#[derive(Debug, Default, Clone, Copy)]
pub struct SlantPath {
    /// Median gaseous absorption, dB.
    pub a_gas_db: f64,
    /// Bending angle, rad.
    pub bending_rad: f64,
    /// Ray length, km.
    pub a_km: f64,
    /// Incident angle, rad.
    pub angle_rad: f64,
    /// Excess atmospheric path length, km.
    pub delta_l_km: f64,
}

/// Gaseous attenuation from height `h_1_km` to `h_2_km` for a ray leaving `h_1_km` at `beta_1_rad` from zenith
/// (`SlantPathAttenuation`). A ray pointing below the horizon (`beta_1_rad > π/2`) is traced from its grazing height
/// in both directions (P.676, section 2.2.2).
pub fn slant_path_attenuation(f_ghz: f64, h_1_km: f64, h_2_km: f64, beta_1_rad: f64) -> SlantPath {
    if beta_1_rad > PI / 2.0 {
        // Refractive index at h_1.
        let n_1 = refractive_index(
            global_pressure(h_1_km),
            global_temperature(h_1_km),
            global_wet_pressure(h_1_km),
        );

        // Binary search for the grazing height h_G, starting half way between h_1 and the surface.
        let mut h_g_km = h_1_km;
        let mut delta = h_1_km / 2.0;
        let mut diff = 100.0;
        loop {
            if diff > 0.0 {
                h_g_km -= delta;
            } else {
                h_g_km += delta;
            }
            delta /= 2.0;

            let n_g = refractive_index(
                global_pressure(h_g_km),
                global_temperature(h_g_km),
                global_wet_pressure(h_g_km),
            );
            let grazing_term = n_g * (A_0_KM + h_g_km);
            let start_term = n_1 * (A_0_KM + h_1_km) * beta_1_rad.sin();
            diff = grazing_term - start_term;
            if diff.abs() <= 0.001 || diff.is_nan() {
                break;
            }
        }

        // Converged on h_G: trace in both directions from the grazing angle.
        let beta_graze_rad = PI / 2.0;
        let result_1 = ray_trace(f_ghz, h_g_km, h_1_km, beta_graze_rad);
        let result_2 = ray_trace(f_ghz, h_g_km, h_2_km, beta_graze_rad);

        SlantPath {
            angle_rad: result_2.angle_rad,
            a_gas_db: result_1.a_gas_db + result_2.a_gas_db,
            a_km: result_1.a_km + result_2.a_km,
            bending_rad: result_1.bending_rad + result_2.bending_rad,
            delta_l_km: result_1.delta_l_km + result_2.delta_l_km,
        }
    } else {
        ray_trace(f_ghz, h_1_km, h_2_km, beta_1_rad)
    }
}
