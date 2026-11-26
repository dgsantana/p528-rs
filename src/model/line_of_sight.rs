//! The line-of-sight region (NTIA `LineOfSight.cpp`, P.528 Annex 2, Step 4 and sections 6, 8 and 13).

use std::f64::consts::PI;

use super::max;
use super::nakagami_rice::{find_k_for_y_pi_at_99_percent, nakagami_rice};
use super::path_loss::{Path, path_loss};
use super::ray_optics::{LineOfSightParams, ray_optics};
use super::terminal::Terminal;
use super::variability::{combine_distributions, long_term_variability};
use crate::Polarization;
use crate::p676::SlantPathLayers;
use crate::p835::OutsideAtmosphere;

/// What `line_of_sight` reports back to the caller (the fields of NTIA's `Result` it writes, and `K_LOS`).
#[derive(Debug, Default, Clone, Copy)]
pub struct LineOfSightResult {
    /// Path distance the ray geometry converged on, km.
    pub d_km: f64,
    /// Basic transmission loss, dB.
    pub a_db: f64,
    /// Free-space loss, dB.
    pub a_fs_db: f64,
    /// Atmospheric absorption, dB.
    pub a_a_db: f64,
    /// Elevation angle of the ray at the low terminal, rad.
    pub theta_h1_rad: f64,
    /// `K_LOS`, the Nakagami-Rice K for the line-of-sight path, dB.
    pub k_los_db: f64,
}

/// Ray geometry for grazing angle `psi`, with `a_los_db` zero (NTIA computes it into an uninitialised temporary and
/// does not read that field).
fn geometry(terminal_1: &Terminal, terminal_2: &Terminal, psi: f64) -> LineOfSightParams {
    let mut params = LineOfSightParams::default();
    ray_optics(terminal_1, terminal_2, psi, &mut params);
    params
}

/// Grazing angle at which the path distance is `d_km`, to within 1 m, by bisection (`FindPsiAtDistance`).
#[allow(
    clippy::neg_cmp_op_on_partial_ord,
    reason = "a NaN ends the loop, as in NTIA's `while`"
)]
fn find_psi_at_distance(d_km: f64, terminal_1: &Terminal, terminal_2: &Terminal) -> f64 {
    if d_km == 0.0 {
        return PI / 2.0;
    }

    // Start at the mid-point.
    let mut psi = PI / 2.0;
    let mut delta_psi = -PI / 4.0;
    loop {
        psi += delta_psi;
        let d_psi_km = geometry(terminal_1, terminal_2, psi).d_km;

        delta_psi = if d_psi_km > d_km {
            delta_psi.abs() / 2.0
        } else {
            -delta_psi.abs() / 2.0
        };

        if !((d_km - d_psi_km).abs() > 1e-3 && delta_psi.abs() > 1e-12) {
            return psi;
        }
    }
}

/// Bisection on the grazing angle until the path length difference is within `terminate` of `delta_r_km`, returning
/// the ray geometry it stopped at (`FindPsiAtDeltaR` and `FindDistanceAtDeltaR`, which differ only in what they
/// return).
///
/// NTIA's loop has no other exit. Once the step reaches zero no further iteration changes anything, so this returns
/// there instead of spinning; where NTIA's loop terminates, the result is identical.
#[allow(
    clippy::neg_cmp_op_on_partial_ord,
    reason = "a NaN ends the loop, as in NTIA's `while`"
)]
fn find_at_delta_r(
    delta_r_km: f64,
    terminal_1: &Terminal,
    terminal_2: &Terminal,
    terminate: f64,
) -> (f64, LineOfSightParams) {
    let mut psi = PI / 2.0;
    let mut delta_psi = -PI / 4.0;
    loop {
        psi += delta_psi;
        let params = geometry(terminal_1, terminal_2, psi);

        delta_psi = if params.delta_r_km > delta_r_km {
            -delta_psi.abs() / 2.0
        } else {
            delta_psi.abs() / 2.0
        };

        if !((params.delta_r_km - delta_r_km).abs() > terminate) || delta_psi == 0.0 {
            return (psi, params);
        }
    }
}

/// The part of NTIA's `LineOfSight` that depends only on the terminals and the frequency: the two-ray limits, `d_0`
/// and the loss at `d_0`.
#[derive(Debug, Clone, Copy)]
pub struct LineOfSight {
    /// Wavelength, km ([Eqn 6-1]).
    lambda_km: f64,
    /// Grazing angle above which the phase lag is ignored.
    psi_limit: f64,
    /// Line-of-sight loss at `d_0`, dB.
    a_d_0_db: f64,
    /// Diffraction loss at `d_ML`, negated as NTIA passes it, dB.
    a_dml_db: f64,
}

impl LineOfSight {
    /// Computes the distance-independent steps and writes `path.d_0_km`. `a_dml_db` is the diffraction loss at
    /// `d_ML`, negated as NTIA passes it.
    pub fn new(
        path: &mut Path,
        terminal_1: &Terminal,
        terminal_2: &Terminal,
        f_mhz: f64,
        a_dml_db: f64,
        polarization: Polarization,
    ) -> Self {
        // 0.2997925 is the speed of light in gigametres per second.
        let lambda_km = 0.2997925 / f_mhz; // [Eqn 6-1]
        let terminate = lambda_km / 1e6;

        // psi_limit, where the free-space model gives way to the two-ray model; lambda / 2 starts the lobe closest
        // to d_ML.
        let (psi_limit, _) = find_at_delta_r(lambda_km / 2.0, terminal_1, terminal_2, terminate);

        // "[d_y6__km] is the largest distance at which a free-space value is obtained in a two-ray model of
        // reflection from a smooth earth with a reflection coefficient of -1" [ES-83-3, page 44].
        let d_y6_km = find_at_delta_r(lambda_km / 6.0, terminal_1, terminal_2, terminate)
            .1
            .d_km;

        // d_0. IF-73's values (d_d in IF-77) were too small when both antennas are low; this heuristic fixes that
        // [Eqns 8-2 and 8-3].
        path.d_0_km = if terminal_1.d_r_km >= path.d_d_km || path.d_d_km >= path.d_ml_km {
            if terminal_1.d_r_km > d_y6_km || d_y6_km > path.d_ml_km {
                terminal_1.d_r_km
            } else {
                d_y6_km
            }
        } else if path.d_d_km < d_y6_km && d_y6_km < path.d_ml_km {
            d_y6_km
        } else {
            path.d_d_km
        };

        // Walk d_0 forward 1 m at a time, as far as it goes without leaving the line-of-sight region.
        let mut d_temp_km = path.d_0_km;
        loop {
            let psi = find_psi_at_distance(d_temp_km, terminal_1, terminal_2);
            let los_result = geometry(terminal_1, terminal_2, psi);

            if los_result.d_km >= path.d_0_km || (d_temp_km + 0.001) >= path.d_ml_km {
                path.d_0_km = los_result.d_km;
                break;
            }
            d_temp_km += 0.001;
        }

        // Loss at d_0.
        let psi_d0 = find_psi_at_distance(path.d_0_km, terminal_1, terminal_2);
        let mut at_d0 = geometry(terminal_1, terminal_2, psi_d0);
        path_loss(
            psi_d0,
            path,
            f_mhz,
            psi_limit,
            a_dml_db,
            0.0,
            polarization,
            &mut at_d0,
        );

        Self {
            lambda_km,
            psi_limit,
            a_d_0_db: at_d0.a_los_db,
            a_dml_db,
        }
    }

    /// Loss over a line-of-sight path of `d_km` at time percentage `p` (the rest of NTIA's `LineOfSight`).
    #[allow(clippy::too_many_arguments)]
    pub fn at_distance(
        &self,
        path: &Path,
        terminal_1: &Terminal,
        terminal_2: &Terminal,
        absorption: &SlantPathLayers,
        f_mhz: f64,
        p: f64,
        d_km: f64,
        polarization: Polarization,
    ) -> Result<LineOfSightResult, OutsideAtmosphere> {
        let lambda_km = self.lambda_km;

        let psi = find_psi_at_distance(d_km, terminal_1, terminal_2);
        let mut los_params = geometry(terminal_1, terminal_2, psi);
        let r_tg = path_loss(
            psi,
            path,
            f_mhz,
            self.psi_limit,
            self.a_dml_db,
            self.a_d_0_db,
            polarization,
            &mut los_params,
        );

        // Atmospheric absorption.
        let slant = absorption.attenuation(PI / 2.0 - los_params.theta_h1_rad)?;
        let a_a_db = slant.a_gas_db;

        // Free-space loss, [Eqn 6-4].
        let a_fs_db = 20.0 * los_params.r_0_km.log10() + 20.0 * f_mhz.log10() + 32.45;

        // Variability. [Eqn 13-1]
        let f_theta_h = if los_params.theta_h1_rad <= 0.0 {
            1.0
        } else if los_params.theta_h1_rad >= 1.0 {
            0.0
        } else {
            max(
                0.5 - (1.0 / PI) * (20.0 * (32.0 * los_params.theta_h1_rad).log10()).atan(),
                0.0,
            )
        };

        let (y_e_db, _) = long_term_variability(
            terminal_1.d_r_km,
            terminal_2.d_r_km,
            d_km,
            f_mhz,
            p,
            f_theta_h,
            los_params.a_los_db,
        );
        let (y_e_50_db, a_y) = long_term_variability(
            terminal_1.d_r_km,
            terminal_2.d_r_km,
            d_km,
            f_mhz,
            50.0,
            f_theta_h,
            los_params.a_los_db,
        );

        // [Eqn 13-2]
        let f_ay = if a_y <= 0.0 {
            1.0
        } else if a_y >= 9.0 {
            0.1
        } else {
            (1.1 + (0.9 * ((a_y / 9.0) * PI).cos())) / 2.0
        };

        // [Eqn 175]
        let f_delta_r = if los_params.delta_r_km >= (lambda_km / 2.0) {
            1.0
        } else if los_params.delta_r_km <= lambda_km / 6.0 {
            0.1
        } else {
            0.5 * (1.1
                - (0.9
                    * (((3.0 * PI) / lambda_km) * (los_params.delta_r_km - (lambda_km / 6.0)))
                        .cos()))
        };

        let r_s = r_tg * f_delta_r * f_ay; // [Eqn 13-4]

        let y_pi_99_db = 10.0 * (f_mhz * slant.a_km.powf(3.0)).log10() - 84.26; // [Eqn 13-5]
        let k_t = find_k_for_y_pi_at_99_percent(y_pi_99_db);

        let w_a = 10.0_f64.powf(k_t / 10.0); // [Eqn 13-6]
        let w_r = r_s.powf(2.0) + 0.01_f64.powf(2.0); // [Eqn 13-7]
        let w = w_r + w_a; // [Eqn 13-8]

        // [Eqn 13-9]
        let k_los_db = if w <= 0.0 {
            -40.0
        } else {
            let k = 10.0 * w.log10();
            if k < -40.0 { -40.0 } else { k }
        };

        let y_pi_50_db = 0.0; // zero mean
        let y_pi_db = nakagami_rice(k_los_db, p);

        let y_total_db = -combine_distributions(y_e_50_db, y_e_db, y_pi_50_db, y_pi_db, p);

        Ok(LineOfSightResult {
            d_km: los_params.d_km,
            a_db: a_fs_db + a_a_db - los_params.a_los_db + y_total_db,
            a_fs_db,
            a_a_db,
            theta_h1_rad: los_params.theta_h1_rad,
            k_los_db,
        })
    }
}
