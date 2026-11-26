//! The step-by-step method of Recommendation ITU-R P.528-5 Annex 2, ported from NTIA/ITS's reference implementation.
//! Each submodule maps to one or two of NTIA's source files, named in its documentation.

mod diffraction;
mod line_of_sight;
mod nakagami_rice;
mod p528;
mod path_loss;
mod ray_optics;
mod terminal;
mod transhorizon;
mod troposcatter;
mod validate;
mod variability;

pub use p528::p528;

/// Mean earth radius, km (`a_0__km`).
pub const A_0_KM: f64 = 6371.0;
/// Effective earth radius, km (`a_e__km`).
pub const A_E_KM: f64 = 9257.0;
/// Surface refractivity, N-units (`N_s`).
const N_S: f64 = 341.0;
/// Relative permittivity of the ground (`epsilon_r`).
const EPSILON_R: f64 = 15.0;
/// Conductivity of the ground, S/m (`sigma`).
const SIGMA: f64 = 0.005;
/// One third (`THIRD`).
const THIRD: f64 = 1.0 / 3.0;

/// NTIA's `MIN` macro: `y` when `x < y` is false, so a NaN in either place gives `y`.
pub fn min(x: f64, y: f64) -> f64 {
    if x < y { x } else { y }
}

/// NTIA's `MAX` macro: `y` when `x > y` is false, so a NaN in either place gives `y`.
pub fn max(x: f64, y: f64) -> f64 {
    if x > y { x } else { y }
}

/// Straight-line interpolation through `(x1, y1)` and `(x2, y2)` at `x` (`LinearInterpolation`).
fn linear_interpolation(x1: f64, y1: f64, x2: f64, y2: f64, x: f64) -> f64 {
    (y1 * (x2 - x) + y2 * (x - x1)) / (x2 - x1)
}
