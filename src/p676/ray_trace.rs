//! The layered ray trace through the reference atmosphere (NTIA `RayTrace.cpp`, P.676 Annex 1, section 2.2).

use super::slant_path::SlantPath;
use super::specific_attenuation::{global_wet_pressure, refractive_index, specific_attenuation};
use crate::model::{A_0_KM, min};
use crate::p835::{global_pressure, global_temperature};

/// Thickness of layer `i`, km (P.676 equation 14).
fn layer_thickness(m: f64, i: i32) -> f64 {
    m * (f64::from(i - 1) / 100.0).exp()
}

/// Refractive index and specific attenuation (dB/km) at height `h_km` (`GetLayerProperties`).
fn layer_properties(f_ghz: f64, h_km: f64) -> (f64, f64) {
    let t_kelvin = global_temperature(h_km);
    let p_hpa = global_pressure(h_km);
    let e_hpa = global_wet_pressure(h_km);
    (
        refractive_index(p_hpa, t_kelvin, e_hpa),
        specific_attenuation(f_ghz, t_kelvin, e_hpa, p_hpa),
    )
}

/// Layer index for height `h_km` before rounding (P.676 equations 16a and 16b).
fn layer_index(h_km: f64) -> f64 {
    100.0 * (1e4 * h_km * ((1.0_f64 / 100.0).exp() - 1.0) + 1.0).ln() + 1.0
}

/// Traces a ray from `h_1_km` to `h_2_km`, leaving at `beta_1_rad` from zenith, summing gaseous attenuation, ray
/// length, bending and excess path length layer by layer (`RayTrace`).
pub fn ray_trace(f_ghz: f64, h_1_km: f64, h_2_km: f64, beta_1_rad: f64) -> SlantPath {
    // P.676 equations 16a to 16c. NTIA converts the rounded doubles to `int`.
    let i_lower = layer_index(h_1_km).floor() as i32;
    let i_upper = layer_index(h_2_km).ceil() as i32;
    let m = (((2.0_f64 / 100.0).exp() - (1.0_f64 / 100.0).exp())
        / ((f64::from(i_upper) / 100.0).exp() - (f64::from(i_lower) / 100.0).exp()))
        * (h_2_km - h_1_km);
    let step = (1.0_f64 / 100.0).exp() - 1.0;
    let lower = (f64::from(i_lower - 1) / 100.0).exp();

    let mut result = SlantPath::default();
    let mut alpha_i_rad = beta_1_rad;

    // Starting layer, at h_1. (NTIA writes its height as h_1 + m * (lower - lower) / step; that differs from h_1 only
    // when m is not finite, and then no layer is traced.)
    let mut delta_i_km = layer_thickness(m, i_lower);
    let (mut n_i, mut gamma_i) = layer_properties(f_ghz, h_1_km + delta_i_km / 2.0);
    let mut r_i_km = A_0_KM + h_1_km;

    // Bottom-layer properties, for the alpha and beta calculations.
    let r_1_km = r_i_km;
    let n_1 = n_i;

    // Summation of P.676 equation 13.
    for i in i_lower..i_upper {
        let delta_ii_km = layer_thickness(m, i + 1);
        let h_ii_km = h_1_km + m * (((f64::from(i) / 100.0).exp() - lower) / step);
        let (n_ii, gamma_ii) = layer_properties(f_ghz, h_ii_km + delta_ii_km / 2.0);
        let r_ii_km = A_0_KM + h_ii_km;
        delta_i_km = layer_thickness(m, i);

        // P.676 equation 19b.
        let beta_i_rad = min(1.0, (n_1 * r_1_km) / (n_i * r_i_km) * beta_1_rad.sin()).asin();
        // Entry angle into the layer interface, P.676 equation 18a.
        alpha_i_rad = min(1.0, (n_1 * r_1_km) / (n_i * r_ii_km) * beta_1_rad.sin()).asin();
        // Path length through layer i, P.676 equation 17.
        let a_i_km = -r_i_km * beta_i_rad.cos()
            + (r_i_km.powf(2.0) * beta_i_rad.cos().powf(2.0)
                + 2.0 * r_i_km * delta_i_km
                + delta_i_km.powf(2.0))
            .sqrt();

        result.a_km += a_i_km;
        result.a_gas_db += a_i_km * gamma_i;
        result.delta_l_km += a_i_km * (n_i - 1.0); // P.676 equation 23

        let beta_ii_rad = (n_i / n_ii * alpha_i_rad.sin()).asin();

        // Bending angle, P.676 equation 22a; the summation stops at i_max - 1.
        if i != i_upper - 1 {
            result.bending_rad += beta_ii_rad - alpha_i_rad;
        }

        n_i = n_ii;
        gamma_i = gamma_ii;
        r_i_km = r_ii_km;
    }

    result.angle_rad = alpha_i_rad;
    result
}
