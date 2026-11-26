//! The layered ray trace through the reference atmosphere (NTIA `RayTrace.cpp`, P.676 Annex 1, section 2.2).
//!
//! The layer heights depend only on the two end heights, not on the ray's angle. [`SlantPathLayers`] computes the
//! layers' properties, where the time goes, once for a pair of heights and traces any angle through them; the
//! arithmetic is the same as tracing directly, so the results are identical.

use std::f64::consts::PI;

use super::slant_path::{SlantPath, slant_path_attenuation};
use super::specific_attenuation::{global_wet_pressure, refractive_index, specific_attenuation};
use crate::model::{A_0_KM, min};
use crate::p835::{OutsideAtmosphere, global_pressure, global_temperature};

/// Thickness of layer `i`, km (P.676 equation 14).
fn layer_thickness(m: f64, i: i32) -> f64 {
    m * (f64::from(i - 1) / 100.0).exp()
}

/// Refractive index and specific attenuation (dB/km) at height `h_km` (`GetLayerProperties`).
fn layer_properties(f_ghz: f64, h_km: f64) -> Result<(f64, f64), OutsideAtmosphere> {
    let t_kelvin = global_temperature(h_km)?;
    let p_hpa = global_pressure(h_km)?;
    let e_hpa = global_wet_pressure(h_km)?;
    Ok((
        refractive_index(p_hpa, t_kelvin, e_hpa),
        specific_attenuation(f_ghz, t_kelvin, e_hpa, p_hpa),
    ))
}

/// Layer index for height `h_km` before rounding (P.676 equations 16a and 16b).
fn layer_index(h_km: f64) -> f64 {
    100.0 * (1e4 * h_km * ((1.0_f64 / 100.0).exp() - 1.0) + 1.0).ln() + 1.0
}

/// The layers between two heights (P.676 equations 14 and 16).
#[derive(Debug, Clone, Copy)]
struct LayerGrid {
    h_1_km: f64,
    i_lower: i32,
    i_upper: i32,
    m: f64,
    step: f64,
    lower: f64,
}

impl LayerGrid {
    fn new(h_1_km: f64, h_2_km: f64) -> Self {
        // P.676 equations 16a to 16c. NTIA converts the rounded doubles to `int`.
        let i_lower = layer_index(h_1_km).floor() as i32;
        let i_upper = layer_index(h_2_km).ceil() as i32;
        let m = (((2.0_f64 / 100.0).exp() - (1.0_f64 / 100.0).exp())
            / ((f64::from(i_upper) / 100.0).exp() - (f64::from(i_lower) / 100.0).exp()))
            * (h_2_km - h_1_km);
        Self {
            h_1_km,
            i_lower,
            i_upper,
            m,
            step: (1.0_f64 / 100.0).exp() - 1.0,
            lower: (f64::from(i_lower - 1) / 100.0).exp(),
        }
    }

    /// Height of the bottom of layer `j`, km. NTIA writes the first layer's as `h_1 + m * (lower - lower) / step`,
    /// which differs from `h_1` only when `m` is not finite, and then no layer is traced.
    fn bottom(&self, j: i32) -> f64 {
        if j == self.i_lower {
            self.h_1_km
        } else {
            self.h_1_km + self.m * (((f64::from(j - 1) / 100.0).exp() - self.lower) / self.step)
        }
    }

    /// Height of the middle of layer `j`, km, where its properties are taken.
    fn midpoint(&self, j: i32) -> f64 {
        self.bottom(j) + layer_thickness(self.m, j) / 2.0
    }

    /// Traces a ray leaving the bottom at `beta_1_rad` from zenith through the layers, given their properties
    /// (refractive index, specific attenuation) by index (`RayTrace`).
    fn trace(
        &self,
        beta_1_rad: f64,
        mut properties: impl FnMut(i32) -> Result<(f64, f64), OutsideAtmosphere>,
    ) -> Result<SlantPath, OutsideAtmosphere> {
        let mut result = SlantPath::default();
        let mut alpha_i_rad = beta_1_rad;

        // No layer to trace: NTIA reads the first layer's properties (at a height that is not finite here) but uses
        // none of them.
        if self.i_lower >= self.i_upper {
            result.angle_rad = alpha_i_rad;
            return Ok(result);
        }

        let (mut n_i, mut gamma_i) = properties(self.i_lower)?;
        let mut r_i_km = A_0_KM + self.h_1_km;

        // Bottom-layer properties, for the alpha and beta calculations.
        let r_1_km = r_i_km;
        let n_1 = n_i;

        // Summation of P.676 equation 13.
        for i in self.i_lower..self.i_upper {
            let (n_ii, gamma_ii) = properties(i + 1)?;
            let r_ii_km = A_0_KM + self.bottom(i + 1);
            let delta_i_km = layer_thickness(self.m, i);

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
            if i != self.i_upper - 1 {
                result.bending_rad += beta_ii_rad - alpha_i_rad;
            }

            n_i = n_ii;
            gamma_i = gamma_ii;
            r_i_km = r_ii_km;
        }

        result.angle_rad = alpha_i_rad;
        Ok(result)
    }
}

/// Traces a ray from `h_1_km` to `h_2_km`, leaving at `beta_1_rad` from zenith, summing gaseous attenuation, ray
/// length, bending and excess path length layer by layer (`RayTrace`). Fails when a layer lies outside the reference
/// atmosphere.
pub fn ray_trace(
    f_ghz: f64,
    h_1_km: f64,
    h_2_km: f64,
    beta_1_rad: f64,
) -> Result<SlantPath, OutsideAtmosphere> {
    let grid = LayerGrid::new(h_1_km, h_2_km);
    grid.trace(beta_1_rad, |j| layer_properties(f_ghz, grid.midpoint(j)))
}

/// Slant-path attenuation between two fixed heights for any angle, with the layers' properties computed once.
#[derive(Debug, Clone)]
pub struct SlantPathLayers {
    f_ghz: f64,
    h_1_km: f64,
    h_2_km: f64,
    grid: LayerGrid,
    /// Refractive index and specific attenuation of layers `i_lower..=i_upper`.
    properties: Vec<(f64, f64)>,
}

impl SlantPathLayers {
    /// The layers between `h_1_km` and `h_2_km` for frequency `f_ghz`. Fails when a layer lies outside the
    /// reference atmosphere.
    pub fn create(f_ghz: f64, h_1_km: f64, h_2_km: f64) -> Result<Self, OutsideAtmosphere> {
        let grid = LayerGrid::new(h_1_km, h_2_km);
        // Empty when no layer is traced (`trace` reads none then).
        let properties = (grid.i_lower..=grid.i_upper)
            .filter(|_| grid.i_lower < grid.i_upper)
            .map(|j| layer_properties(f_ghz, grid.midpoint(j)))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            f_ghz,
            h_1_km,
            h_2_km,
            grid,
            properties,
        })
    }

    /// The same as `slant_path_attenuation(f_ghz, h_1_km, h_2_km, beta_1_rad)`. A ray above the horizon reuses the
    /// stored layers; one below it starts from its own grazing height, so it is traced in full.
    pub fn attenuation(&self, beta_1_rad: f64) -> Result<SlantPath, OutsideAtmosphere> {
        if beta_1_rad > PI / 2.0 {
            return slant_path_attenuation(self.f_ghz, self.h_1_km, self.h_2_km, beta_1_rad);
        }
        let first = self.grid.i_lower;
        self.grid
            .trace(beta_1_rad, |j| Ok(self.properties[(j - first) as usize]))
    }
}
