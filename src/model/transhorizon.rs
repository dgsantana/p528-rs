//! Search past the horizon for the crossover from diffraction to troposcatter (NTIA `TranshorizonSearch.cpp`,
//! P.528 Annex 2, Step 6).

use super::path_loss::Path;
use super::terminal::Terminal;
use super::troposcatter::troposcatter;
use crate::Warnings;

/// How the transhorizon loss is chosen beyond the crossover distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Case {
    /// The lower of diffraction and troposcatter (`CASE_1`). Also the outcome when the search finds no crossover:
    /// NTIA then sets `CONST_MODE__DIFFRACTION`, which has `CASE_1`'s value.
    LowerLoss,
    /// Troposcatter, with the diffraction line adjusted to meet it (`CASE_2`).
    Troposcatter,
}

/// The diffraction line and crossover the search settles on.
#[derive(Debug, Clone, Copy)]
pub struct Transhorizon {
    /// Slope of the diffraction line, dB/km.
    pub m_d: f64,
    /// Intercept of the diffraction line, dB.
    pub a_d0: f64,
    /// Crossover distance, km.
    pub d_crx_km: f64,
    /// How the loss is chosen beyond `d_crx_km`.
    pub case: Case,
}

/// Steps 1 km at a time from 3 km beyond `d_ML`, for at most 100 km, until the troposcatter line is no steeper than the
/// diffraction line (`TranshorizonSearch`). Sets `Warnings::DIFFRACTION_TROPO_REGION` when it never is.
#[allow(clippy::too_many_arguments)]
pub fn transhorizon_search(
    path: &Path,
    terminal_1: &Terminal,
    terminal_2: &Terminal,
    f_mhz: f64,
    a_dml_db: f64,
    m_d: f64,
    a_d0: f64,
    warnings: &mut Warnings,
) -> Transhorizon {
    let mut k = 0;

    // Step 6.1: d' and d" [Eqns 3-8 and 3-9].
    let mut d_search_km = [path.d_ml_km + 3.0, path.d_ml_km + 2.0];
    let mut a_s_db = [0.0_f64; 2];

    let search_limit = 100; // 100 km beyond the starting point

    for _ in 0..search_limit {
        a_s_db[1] = a_s_db[0];

        // Step 6.2
        a_s_db[0] = troposcatter(terminal_1, terminal_2, d_search_km[0], f_mhz).a_s_db;

        // Below 20 dB the troposcatter result is outside the valid part of the model; and two points are needed to
        // draw a line.
        let valid = a_s_db[0] >= 20.0;
        if valid {
            k += 1;
        }
        if !valid || k <= 1 {
            d_search_km[1] = d_search_km[0];
            d_search_km[0] += 1.0;
            continue;
        }

        // Step 6.3
        let m_s = (a_s_db[0] - a_s_db[1]) / (d_search_km[0] - d_search_km[1]); // [Eqn 3-10]

        if m_s <= m_d {
            // Step 6.6
            let a_d_db = m_d * d_search_km[1] + a_d0; // [Eqn 3-11]

            return if a_s_db[1] >= a_d_db {
                Transhorizon {
                    m_d,
                    a_d0,
                    d_crx_km: d_search_km[0],
                    case: Case::LowerLoss,
                }
            } else {
                // Adjust the diffraction line to the troposcatter model.
                let m_d = (a_s_db[1] - a_dml_db) / (d_search_km[1] - path.d_ml_km); // [Eqn 3-12]
                let a_d0 = a_s_db[1] - (m_d * d_search_km[1]); // [Eqn 3-13]
                Transhorizon {
                    m_d,
                    a_d0,
                    d_crx_km: d_search_km[0],
                    case: Case::Troposcatter,
                }
            };
        }

        d_search_km[1] = d_search_km[0];
        d_search_km[0] += 1.0;
    }

    // M_s was always greater than M_d: default to the diffraction-only transhorizon model.
    *warnings |= Warnings::DIFFRACTION_TROPO_REGION;
    Transhorizon {
        m_d,
        a_d0,
        d_crx_km: d_search_km[1],
        case: Case::LowerLoss,
    }
}
