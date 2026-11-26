//! The main entry point (NTIA `P528.cpp`, P.528-5 Annex 2, section 3).

use std::f64::consts::PI;

use super::diffraction::smooth_earth_diffraction;
use super::line_of_sight::line_of_sight;
use super::nakagami_rice::nakagami_rice;
use super::path_loss::Path;
use super::ray_optics::LineOfSightParams;
use super::terminal::Terminal;
use super::transhorizon::{Case, transhorizon_search};
use super::troposcatter::troposcatter;
use super::validate::{Valid, validate};
use super::variability::{combine_distributions, long_term_variability};
use super::{A_E_KM, THIRD};
use crate::p676::slant_path_attenuation;
use crate::{P528Error, P528Input, P528Result, PropagationMode, Warnings};

/// Basic transmission loss between a low and a high terminal over a smooth earth, by the step-by-step method of
/// Recommendation ITU-R P.528-5 Annex 2, as NTIA/ITS's reference implementation computes it (`P528`).
///
/// Pure and deterministic: no global state and no allocation, so it can be called from many threads at once.
///
/// # Errors
///
/// [`P528Error`] for an input outside the ranges NTIA's reference accepts, in the order NTIA checks them, and for NaN.
pub fn p528(input: &P528Input) -> Result<P528Result, P528Error> {
    let mut warnings = Warnings::empty();
    if validate(input, &mut warnings)? == Valid::Coincident {
        return Ok(P528Result {
            basic_transmission_loss_db: 0.0,
            free_space_loss_db: 0.0,
            absorption_loss_db: 0.0,
            distance_km: 0.0,
            theta_h1_rad: 0.0,
            mode: PropagationMode::NotSet,
            warnings,
        });
    }

    let P528Input {
        distance_km: d_km,
        h1_m,
        h2_m,
        frequency_mhz: f_mhz,
        time_percent: p,
        polarization,
    } = *input;

    // Step 1: terminal geometries.
    let terminal_1 = Terminal::new(h1_m / 1000.0, f_mhz);
    let terminal_2 = Terminal::new(h2_m / 1000.0, f_mhz);

    // Step 2
    let mut path = Path {
        d_ml_km: terminal_1.d_r_km + terminal_2.d_r_km,
        ..Path::default()
    }; // [Eqn 3-1]

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

    let mut los_params = LineOfSightParams::default();

    // Step 4: a line-of-sight path ends here.
    if path.d_ml_km - d_km > 0.001 {
        let los = line_of_sight(
            &mut path,
            &terminal_1,
            &terminal_2,
            &mut los_params,
            f_mhz,
            -a_dml_db,
            p,
            d_km,
            polarization,
        );
        return Ok(P528Result {
            basic_transmission_loss_db: los.a_db,
            free_space_loss_db: los.a_fs_db,
            absorption_loss_db: los.a_a_db,
            distance_km: los.d_km,
            theta_h1_rad: los.theta_h1_rad,
            mode: PropagationMode::LineOfSight,
            warnings,
        });
    }

    // K_LOS, from the line-of-sight model 1 km inside the horizon.
    let d_inside_km = path.d_ml_km - 1.0;
    let k_los = line_of_sight(
        &mut path,
        &terminal_1,
        &terminal_2,
        &mut los_params,
        f_mhz,
        -a_dml_db,
        p,
        d_inside_km,
        polarization,
    )
    .k_los_db;

    // Step 6: search past the horizon for the crossover between diffraction and troposcatter.
    let search = transhorizon_search(
        &path,
        &terminal_1,
        &terminal_2,
        f_mhz,
        a_dml_db,
        m_d,
        a_d0,
        &mut warnings,
    );

    // Step 7: terrain attenuation A_T.
    // Step 7.1
    let a_d_db = search.m_d * d_km + search.a_d0; // [Eqn 3-14]

    // Step 7.2
    let tropo = troposcatter(&terminal_1, &terminal_2, d_km, f_mhz);

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
    let vertical = slant_path_attenuation(f_mhz / 1000.0, 0.0, tropo.h_v_km, PI / 2.0);
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
        warnings,
    })
}
