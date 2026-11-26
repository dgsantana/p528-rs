//! Specific attenuation and refractive index of an atmospheric layer (NTIA `SpecificAttenuation.cpp`,
//! `Refractivity.cpp`, `LineShapeFactor.cpp`, `NonresonantDebyeAttenuation.cpp`, `RefractiveIndex.cpp`,
//! `GlobalWetPressure.cpp`).

use super::spectroscopy::{OXYGEN, WATER_VAPOUR};
use crate::model::max;
use crate::p835::{
    Profile, global_pressure, global_temperature, global_water_vapour_density,
    water_vapour_density_to_pressure,
};

/// Surface water vapour density of the global reference atmosphere, g/m³ (`RHO_0__M_KG`).
const RHO_0: f64 = 7.5;

/// Water vapour partial pressure, hPa, at height `h_km` (`GlobalWetPressure`).
pub fn global_wet_pressure(h_km: f64) -> Profile {
    let t_kelvin = global_temperature(h_km)?;
    let p_hpa = global_pressure(h_km)?;
    let rho_g_m3 = max(
        global_water_vapour_density(h_km, RHO_0)?,
        2.0 * 10.0_f64.powf(-6.0) * 216.7 * p_hpa / t_kelvin,
    );
    Ok(water_vapour_density_to_pressure(rho_g_m3, t_kelvin))
}

/// Refractive index of the reference atmosphere at height `h_km`.
pub fn refractive_index_at(h_km: f64) -> Profile {
    Ok(refractive_index(
        global_pressure(h_km)?,
        global_temperature(h_km)?,
        global_wet_pressure(h_km)?,
    ))
}

/// Refractive index from dry pressure, temperature and water vapour pressure (`RefractiveIndex`).
pub fn refractive_index(p_hpa: f64, t_kelvin: f64, e_hpa: f64) -> f64 {
    let n_dry = 77.6 * p_hpa / t_kelvin;
    let n_wet = 72.0 * e_hpa / t_kelvin + 3.75e5 * e_hpa / t_kelvin.powf(2.0);
    let n = n_dry + n_wet;
    1.0 + n * 10.0_f64.powf(-6.0)
}

/// Specific attenuation, dB/km, of oxygen and water vapour together (`SpecificAttenuation`, P.676 equation 1).
pub fn specific_attenuation(f_ghz: f64, t_kelvin: f64, e_hpa: f64, p_hpa: f64) -> f64 {
    let gamma_o = 0.1820 * f_ghz * oxygen_refractivity(f_ghz, t_kelvin, e_hpa, p_hpa);
    let gamma_w = 0.1820 * f_ghz * water_vapour_refractivity(f_ghz, t_kelvin, e_hpa, p_hpa);
    gamma_o + gamma_w
}

/// Line-shape factor (`LineShapeFactor`, P.676 equation 5).
fn line_shape_factor(f_ghz: f64, f_i_ghz: f64, delta_f_ghz: f64, delta: f64) -> f64 {
    let term1 = f_ghz / f_i_ghz;
    let term2 = (delta_f_ghz - delta * (f_i_ghz - f_ghz))
        / ((f_i_ghz - f_ghz).powf(2.0) + delta_f_ghz.powf(2.0));
    let term3 = (delta_f_ghz - delta * (f_i_ghz + f_ghz))
        / ((f_i_ghz + f_ghz).powf(2.0) + delta_f_ghz.powf(2.0));
    term1 * (term2 + term3)
}

/// Dry continuum from the non-resonant Debye spectrum of oxygen and pressure-induced nitrogen absorption
/// (`NonresonantDebyeAttenuation`, P.676 equations 8 and 9).
fn nonresonant_debye_attenuation(f_ghz: f64, e_hpa: f64, p_hpa: f64, theta: f64) -> f64 {
    let d = 5.6e-4 * (p_hpa + e_hpa) * theta.powf(0.8);
    let frac_1 = 6.14e-5 / (d * (1.0 + (f_ghz / d).powf(2.0)));
    let frac_2 = (1.4e-12 * p_hpa * theta.powf(1.5)) / (1.0 + 1.9e-5 * f_ghz.powf(1.5));
    f_ghz * p_hpa * theta.powf(2.0) * (frac_1 + frac_2)
}

/// Imaginary part of the oxygen refractivity (`OxygenRefractivity`, P.676 equations 2a, 3, 6a, 6b, 7).
fn oxygen_refractivity(f_ghz: f64, t_kelvin: f64, e_hpa: f64, p_hpa: f64) -> f64 {
    let theta = 300.0 / t_kelvin;
    let mut n = 0.0;
    for line in &OXYGEN {
        let [f_0, a_1, a_2, a_3, a_4, a_5, a_6] = *line;
        let s_i = a_1 * 1e-7 * p_hpa * theta.powf(3.0) * (a_2 * (1.0 - theta)).exp();
        let delta_f_ghz = a_3 * 1e-4 * (p_hpa * theta.powf(0.8 - a_4) + 1.1 * e_hpa * theta);
        // Zeeman splitting of the oxygen lines.
        let delta_f_ghz = (delta_f_ghz.powf(2.0) + 2.25e-6).sqrt();
        // Interference between oxygen lines.
        let delta = (a_5 + a_6 * theta) * 1e-4 * (p_hpa + e_hpa) * theta.powf(0.8);
        n += s_i * line_shape_factor(f_ghz, f_0, delta_f_ghz, delta);
    }
    n + nonresonant_debye_attenuation(f_ghz, e_hpa, p_hpa, theta)
}

/// Imaginary part of the water vapour refractivity (`WaterVapourRefractivity`, P.676 equations 2b, 3, 6a, 6b).
fn water_vapour_refractivity(f_ghz: f64, t_kelvin: f64, e_hpa: f64, p_hpa: f64) -> f64 {
    let theta = 300.0 / t_kelvin;
    let mut n_w = 0.0;
    for line in &WATER_VAPOUR {
        let [f_0, b_1, b_2, b_3, b_4, b_5, b_6] = *line;
        let s_i = 0.1 * b_1 * e_hpa * theta.powf(3.5) * (b_2 * (1.0 - theta)).exp();
        let delta_f_ghz = 1e-4 * b_3 * (p_hpa * theta.powf(b_4) + b_5 * e_hpa * theta.powf(b_6));
        // Doppler broadening of the water vapour lines.
        let term1 = 0.217 * delta_f_ghz.powf(2.0) + (2.1316e-12 * f_0.powf(2.0) / theta);
        let delta_f_ghz = 0.535 * delta_f_ghz + term1.sqrt();
        n_w += s_i * line_shape_factor(f_ghz, f_0, delta_f_ghz, 0.0);
    }
    n_w
}
