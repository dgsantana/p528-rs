//! Tropospheric scatter (NTIA `Troposcatter.cpp`, P.528 Annex 2, section 11).

use super::terminal::Terminal;
use super::{A_0_KM, A_E_KM, N_S, min};

/// Troposcatter geometry and loss at one distance (the fields of NTIA's `TroposcatterParams` the model reads).
#[derive(Debug, Default, Clone, Copy)]
pub struct TroposcatterParams {
    /// Height of the common-volume cross-over point, km.
    pub h_v_km: f64,
    /// Scattering angle, rad.
    pub theta_s: f64,
    /// Troposcatter loss, dB.
    pub a_s_db: f64,
}

/// Troposcatter loss over a path of `d_km` (`Troposcatter`). Within the line-of-sight distance everything is zero.
pub fn troposcatter(
    terminal_1: &Terminal,
    terminal_2: &Terminal,
    d_km: f64,
    f_mhz: f64,
) -> TroposcatterParams {
    let d_s_km = d_km - terminal_1.d_r_km - terminal_2.d_r_km; // [Eqn 11-2]

    if d_s_km <= 0.0 {
        return TroposcatterParams::default();
    }

    // Geometric parameters.
    let d_z_km = 0.5 * d_s_km; // [Eqn 11-6]

    let a_m = 1.0 / A_0_KM; // [Eqn 11-7]
    let d_n = a_m - (1.0 / A_E_KM); // [Eqn 11-8]
    let gamma_e_km = (N_S * 1e-6) / d_n; // [Eqn 11-9]

    let z_a_km = 1.0 / (2.0 * A_E_KM) * (d_z_km / 2.0).powf(2.0); // [Eqn 11-10]
    let z_b_km = 1.0 / (2.0 * A_E_KM) * d_z_km.powf(2.0); // [Eqn 11-11]

    let q_o = a_m - d_n; // [Eqn 11-12]

    let q_a = a_m - d_n / min(35.0, z_a_km / gamma_e_km).exp(); // [Eqn 11-13]
    let q_b = a_m - d_n / min(35.0, z_b_km / gamma_e_km).exp(); // [Eqn 11-13]

    let big_z_a_km = (7.0 * q_o + 6.0 * q_a - q_b) * (d_z_km.powf(2.0) / 96.0); // [Eqn 11-14]
    let big_z_b_km = (q_o + 2.0 * q_a) * (d_z_km.powf(2.0) / 6.0); // [Eqn 11-15]

    let big_q_a = a_m - d_n / min(35.0, big_z_a_km / gamma_e_km).exp(); // [Eqn 11-16]
    let big_q_b = a_m - d_n / min(35.0, big_z_b_km / gamma_e_km).exp(); // [Eqn 11-16]

    let h_v_km = (q_o + 2.0 * big_q_a) * (d_z_km.powf(2.0) / 6.0); // [Eqn 11-17]
    let theta_a = (q_o + 4.0 * big_q_a + big_q_b) * d_z_km / 6.0; // [Eqn 11-18]
    let theta_s = 2.0 * theta_a; // [Eqn 11-19]

    // Scattering efficiency term.
    let epsilon_1 = 5.67e-6 * N_S.powf(2.0) - 0.00232 * N_S + 0.031; // [Eqn 11-20]
    let epsilon_2 = 0.0002 * N_S.powf(2.0) - 0.06 * N_S + 6.6; // [Eqn 11-21]

    let gamma = 0.1424 * (1.0 + epsilon_1 / min(35.0, (h_v_km / 4.0).powf(6.0)).exp()); // [Eqn 11-22]

    // [Eqn 11-23]
    let s_e_db = 83.1 - epsilon_2 / (1.0 + 0.07716 * h_v_km.powf(2.0))
        + 20.0 * ((0.1424 / gamma).powf(2.0) * (gamma * h_v_km).exp()).log10();

    // Scattering volume term.
    let x_a1_km2 = terminal_1.h_e_km.powf(2.0)
        + 4.0
            * (A_E_KM + terminal_1.h_e_km)
            * A_E_KM
            * (terminal_1.d_r_km / (A_E_KM * 2.0)).sin().powf(2.0); // [Eqn 11-24]
    let x_a2_km2 = terminal_2.h_e_km.powf(2.0)
        + 4.0
            * (A_E_KM + terminal_2.h_e_km)
            * A_E_KM
            * (terminal_2.d_r_km / (A_E_KM * 2.0)).sin().powf(2.0); // [Eqn 11-24]

    let ell_1_km = x_a1_km2.sqrt() + d_z_km; // [Eqn 11-25]
    let ell_2_km = x_a2_km2.sqrt() + d_z_km; // [Eqn 11-25]
    let ell_km = ell_1_km + ell_2_km; // [Eqn 11-26]

    let s = (ell_1_km - ell_2_km) / ell_km; // [Eqn 11-27]
    let eta = gamma * theta_s * ell_km / 2.0; // [Eqn 11-28]

    let kappa = f_mhz / 0.0477; // [Eqn 11-29]

    let rho_1_km = 2.0 * kappa * theta_s * terminal_1.h_e_km; // [Eqn 11-30]
    let rho_2_km = 2.0 * kappa * theta_s * terminal_2.h_e_km; // [Eqn 11-30]

    let sqrt2 = 2.0_f64.sqrt();

    let a = (1.0 - s.powf(2.0)).powf(2.0); // [Eqn 11-36]

    let x_v1 = (1.0 + s).powf(2.0) * eta; // [Eqn 11-32]
    let x_v2 = (1.0 - s).powf(2.0) * eta; // [Eqn 11-33]

    let q_1 = x_v1.powf(2.0) + rho_1_km.powf(2.0); // [Eqn 11-34]
    let q_2 = x_v2.powf(2.0) + rho_2_km.powf(2.0); // [Eqn 11-35]

    // [Eqn 11-37]
    let b_s = 6.0
        + 8.0 * s.powf(2.0)
        + 8.0 * (1.0 - s) * x_v1.powf(2.0) * rho_1_km.powf(2.0) / q_1.powf(2.0)
        + 8.0 * (1.0 + s) * x_v2.powf(2.0) * rho_2_km.powf(2.0) / q_2.powf(2.0)
        + 2.0
            * (1.0 - s.powf(2.0))
            * (1.0 + 2.0 * x_v1.powf(2.0) / q_1)
            * (1.0 + 2.0 * x_v2.powf(2.0) / q_2);

    // [Eqn 11-38]
    let c_s = 12.0
        * ((rho_1_km + sqrt2) / rho_1_km).powf(2.0)
        * ((rho_2_km + sqrt2) / rho_2_km).powf(2.0)
        * (rho_1_km + rho_2_km)
        / (rho_1_km + rho_2_km + 2.0 * sqrt2);

    let temp =
        (a * eta.powf(2.0) + b_s * eta) * q_1 * q_2 / (rho_1_km.powf(2.0) * rho_2_km.powf(2.0));

    let s_v_db = 10.0 * (temp + c_s).log10();

    TroposcatterParams {
        h_v_km,
        theta_s,
        a_s_db: s_e_db + s_v_db + 10.0 * (kappa * theta_s.powf(3.0) / ell_km).log10(),
    }
}
