//! The main entry point (NTIA `P528.cpp`, P.528-5 Annex 2, section 3), split into the steps that depend only on the
//! terminals and frequency ([`P528Model::create`]) and those that depend on the distance and time percentage
//! ([`P528Model::loss`]).

use std::f64::consts::PI;

use super::diffraction::smooth_earth_diffraction;
use super::line_of_sight::LineOfSight;
use super::nakagami_rice::nakagami_rice;
use super::path_loss::Path;
use super::terminal::Terminal;
use super::transhorizon::{Case, Transhorizon, transhorizon_search};
use super::troposcatter::troposcatter;
use super::variability::{combine_distributions, long_term_variability};
use super::{A_E_KM, THIRD, validate};
use crate::p676::{SlantPathLayers, slant_path_attenuation};
use crate::p835::OutsideAtmosphere;
use crate::{P528Error, P528Input, P528Result, Polarization, PropagationMode, Warnings};

impl From<OutsideAtmosphere> for P528Error {
    fn from(_: OutsideAtmosphere) -> Self {
        Self::OutsideAtmosphere
    }
}

/// P.528 for one pair of terminal heights, one frequency and one polarisation: the work that does not depend on the
/// distance or the time percentage, done once, for [`P528Model::loss`] to evaluate any number of distances.
///
/// Its results are the same as [`p528`]'s for the same inputs. It holds no interior mutability, so it is `Send` and
/// `Sync` and can be shared between threads.
#[derive(Debug, Clone)]
pub struct P528Model {
    h1_m: f64,
    h2_m: f64,
    f_mhz: f64,
    polarization: Polarization,
    /// Height warnings, on every result.
    warnings: Warnings,
    terminal_1: Terminal,
    terminal_2: Terminal,
    path: Path,
    line_of_sight: LineOfSight,
    /// Gaseous absorption between the two terminals, for the line-of-sight ray.
    absorption: SlantPathLayers,
    /// The diffraction line and the crossover to troposcatter.
    transhorizon: Transhorizon,
    /// Warnings from the transhorizon search, on transhorizon results only.
    transhorizon_warnings: Warnings,
    /// `K_LOS`, the Nakagami-Rice K of the line-of-sight model 1 km inside the horizon, dB.
    k_los_db: f64,
}

impl P528Model {
    /// The model for a low terminal at `h1_m`, a high terminal at `h2_m` (metres), frequency `frequency_mhz` and
    /// `polarization`.
    ///
    /// # Errors
    ///
    /// [`P528Error`] for a height or frequency outside the ranges NTIA's reference accepts, in NTIA's order, for a
    /// low terminal above the high one, and for NaN.
    pub fn create(
        h1_m: f64,
        h2_m: f64,
        frequency_mhz: f64,
        polarization: Polarization,
    ) -> Result<Self, P528Error> {
        let warnings = validate::terminals(h1_m, h2_m, frequency_mhz)?;
        let f_mhz = frequency_mhz;

        // Step 1: terminal geometries.
        let terminal_1 = Terminal::create(h1_m / 1000.0, f_mhz)?;
        let terminal_2 = Terminal::create(h2_m / 1000.0, f_mhz)?;

        // Step 2
        let mut path = Path {
            d_ml_km: terminal_1.d_r_km + terminal_2.d_r_km, // [Eqn 3-1]
            ..Path::default()
        };

        // Step 3: the smooth-earth diffraction line.
        // Step 3.1
        let d_3_km = path.d_ml_km + 0.5 * (A_E_KM.powf(2.0) / f_mhz).powf(THIRD); // [Eqn 3-2]
        let d_4_km = path.d_ml_km + 1.5 * (A_E_KM.powf(2.0) / f_mhz).powf(THIRD); // [Eqn 3-3]

        // Step 3.2
        let a_3_db = smooth_earth_diffraction(
            terminal_1.d_r_km,
            terminal_2.d_r_km,
            f_mhz,
            d_3_km,
            polarization,
        );
        let a_4_db = smooth_earth_diffraction(
            terminal_1.d_r_km,
            terminal_2.d_r_km,
            f_mhz,
            d_4_km,
            polarization,
        );

        // Step 3.3
        let m_d = (a_4_db - a_3_db) / (d_4_km - d_3_km); // [Eqn 3-4]
        let a_d0 = a_4_db - m_d * d_4_km; // [Eqn 3-5]

        // Step 3.4
        let a_dml_db = (m_d * path.d_ml_km) + a_d0; // [Eqn 3-6]
        path.d_d_km = -(a_d0 / m_d); // [Eqn 3-7]

        // Step 4, the distance-independent part. NTIA passes the diffraction loss at d_ML negated.
        let line_of_sight = LineOfSight::new(
            &mut path,
            &terminal_1,
            &terminal_2,
            f_mhz,
            -a_dml_db,
            polarization,
        );
        let absorption =
            SlantPathLayers::create(f_mhz / 1000.0, terminal_1.h_r_km, terminal_2.h_r_km)?;

        // K_LOS, from the line-of-sight model 1 km inside the horizon. It does not depend on the time percentage:
        // NTIA computes it from the geometry and the 50 % variability.
        let k_los_db = line_of_sight
            .at_distance(
                &path,
                &terminal_1,
                &terminal_2,
                &absorption,
                f_mhz,
                50.0,
                path.d_ml_km - 1.0,
                polarization,
            )?
            .k_los_db;

        // Step 6: search past the horizon for the crossover between diffraction and troposcatter.
        let mut transhorizon_warnings = Warnings::empty();
        let transhorizon = transhorizon_search(
            &path,
            &terminal_1,
            &terminal_2,
            f_mhz,
            a_dml_db,
            m_d,
            a_d0,
            &mut transhorizon_warnings,
        );

        Ok(Self {
            h1_m,
            h2_m,
            f_mhz,
            polarization,
            warnings,
            terminal_1,
            terminal_2,
            path,
            line_of_sight,
            absorption,
            transhorizon,
            transhorizon_warnings,
            k_los_db,
        })
    }

    /// Basic transmission loss over a path of `distance_km` for time percentage `time_percent`.
    ///
    /// # Errors
    ///
    /// [`P528Error`] for a distance or time percentage outside the ranges NTIA's reference accepts, for coincident
    /// terminals (same height, 0 km), for NaN or an infinite distance, and for a path that leaves the reference
    /// atmosphere.
    pub fn loss(&self, distance_km: f64, time_percent: f64) -> Result<P528Result, P528Error> {
        validate::distance(distance_km)?;
        validate::percentage(time_percent, self.h1_m == self.h2_m && distance_km == 0.0)?;
        self.compute(distance_km, time_percent)
    }

    /// Losses for every distance in `distances_km`, computed in parallel; the results are in the order of the
    /// distances.
    #[cfg(feature = "parallel")]
    pub fn par_losses(
        &self,
        distances_km: &[f64],
        time_percent: f64,
    ) -> Vec<Result<P528Result, P528Error>> {
        use rayon::prelude::*;
        distances_km
            .par_iter()
            .map(|&d| self.loss(d, time_percent))
            .collect()
    }

    /// The distance-dependent steps, for validated inputs.
    fn compute(&self, d_km: f64, p: f64) -> Result<P528Result, P528Error> {
        let (f_mhz, polarization) = (self.f_mhz, self.polarization);
        let (terminal_1, terminal_2) = (&self.terminal_1, &self.terminal_2);

        // Step 4: a line-of-sight path ends here.
        if self.path.d_ml_km - d_km > 0.001 {
            let los = self.line_of_sight.at_distance(
                &self.path,
                terminal_1,
                terminal_2,
                &self.absorption,
                f_mhz,
                p,
                d_km,
                polarization,
            )?;
            return Ok(P528Result {
                basic_transmission_loss_db: los.a_db,
                free_space_loss_db: los.a_fs_db,
                absorption_loss_db: los.a_a_db,
                distance_km: los.d_km,
                theta_h1_rad: los.theta_h1_rad,
                mode: PropagationMode::LineOfSight,
                warnings: self.warnings,
            });
        }

        let search = &self.transhorizon;

        // Step 7: terrain attenuation A_T.
        // Step 7.1
        let a_d_db = search.m_d * d_km + search.a_d0; // [Eqn 3-14]

        // Step 7.2
        let tropo = troposcatter(terminal_1, terminal_2, d_km, f_mhz);

        // Step 7.3
        let (a_t_db, mode) = if d_km < search.d_crx_km {
            // Always diffraction short of d_crx.
            (a_d_db, PropagationMode::Diffraction)
        } else {
            match search.case {
                Case::LowerLoss if tropo.a_s_db <= a_d_db => {
                    (tropo.a_s_db, PropagationMode::Troposcatter)
                }
                Case::LowerLoss => (a_d_db, PropagationMode::Diffraction),
                Case::Troposcatter => (tropo.a_s_db, PropagationMode::Troposcatter),
            }
        };

        // Variability. f_theta_h is unity for transhorizon paths.
        let f_theta_h = 1.0;

        // The p % and 50 % values of the long-term variability distribution.
        let (y_e_db, _) = long_term_variability(
            terminal_1.d_r_km,
            terminal_2.d_r_km,
            d_km,
            f_mhz,
            p,
            f_theta_h,
            -a_t_db,
        );
        let (y_e_50_db, _) = long_term_variability(
            terminal_1.d_r_km,
            terminal_2.d_r_km,
            d_km,
            f_mhz,
            50.0,
            f_theta_h,
            -a_t_db,
        );

        // The p % and 50 % values of the Nakagami-Rice distribution.
        let k_los = self.k_los_db;
        let angle = 0.02617993878; // 1.5 degrees
        let k_t_db = if tropo.theta_s >= angle {
            20.0
        } else if tropo.theta_s <= 0.0 {
            k_los
        } else {
            (tropo.theta_s * (20.0 - k_los) / angle) + k_los
        };

        let y_pi_50_db = 0.0; // zero mean
        let y_pi_db = nakagami_rice(k_t_db, p);

        // Combine the long-term and Nakagami-Rice distributions.
        let y_total_db = combine_distributions(y_e_50_db, y_e_db, y_pi_50_db, y_pi_db, p);

        // Atmospheric absorption for the transhorizon path, [Eqn 3-17].
        let vertical = slant_path_attenuation(f_mhz / 1000.0, 0.0, tropo.h_v_km, PI / 2.0)?;
        let a_a_db = terminal_1.a_a_db + terminal_2.a_a_db + 2.0 * vertical.a_gas_db;

        // Free-space loss.
        let r_fs_km = terminal_1.a_km + terminal_2.a_km + 2.0 * vertical.a_km; // [Eqn 3-18]
        let a_fs_db = 20.0 * f_mhz.log10() + 20.0 * r_fs_km.log10() + 32.45; // [Eqn 3-19]

        Ok(P528Result {
            basic_transmission_loss_db: a_fs_db + a_a_db + a_t_db - y_total_db, // [Eqn 3-20]
            free_space_loss_db: a_fs_db,
            absorption_loss_db: a_a_db,
            distance_km: d_km,
            theta_h1_rad: -terminal_1.theta_rad,
            mode,
            warnings: self.warnings | self.transhorizon_warnings,
        })
    }
}

/// Basic transmission loss between a low and a high terminal over a smooth earth, by the step-by-step method of
/// Recommendation ITU-R P.528-5 Annex 2, as NTIA/ITS's reference implementation computes it (`P528`).
///
/// For many distances between the same two heights, build a [`P528Model`] once and call [`P528Model::loss`]: this
/// function builds one per call.
///
/// # Errors
///
/// [`P528Error`] for an input outside the ranges NTIA's reference accepts, in the order NTIA checks them, for
/// coincident terminals, for NaN or an infinite distance, and for a path that leaves the reference atmosphere.
pub fn p528(input: &P528Input) -> Result<P528Result, P528Error> {
    let P528Input {
        distance_km,
        h1_m,
        h2_m,
        frequency_mhz,
        time_percent,
        polarization,
    } = *input;
    // NTIA's order: the distance, the terminals and frequency, the time percentage, then coincident terminals; all
    // checked before the model's work.
    validate::distance(distance_km)?;
    validate::terminals(h1_m, h2_m, frequency_mhz)?;
    validate::percentage(time_percent, h1_m == h2_m && distance_km == 0.0)?;
    let model = P528Model::create(h1_m, h2_m, frequency_mhz, polarization)?;
    model.loss(distance_km, time_percent)
}
