//! Long-term variability and the combination of distributions (NTIA `LongTermVariability.cpp`,
//! `CombineDistributions.cpp`, `InverseComplementaryCumulativeDistributionFunction.cpp`, P.528 Annex 2, section 14).

use super::nakagami_rice::{P, upper_bound};
use super::{THIRD, linear_interpolation, max};

/// Inverse of the complementary cumulative normal distribution, by the rational approximation of Abramowitz and
/// Stegun 26.2.23 (`InverseComplementaryCumulativeDistributionFunction`).
fn inverse_ccdf(q: f64) -> f64 {
    let c_0 = 2.515516;
    let c_1 = 0.802853;
    let c_2 = 0.010328;
    let d_1 = 1.432788;
    let d_2 = 0.189269;
    let d_3 = 0.001308;

    let x = if q > 0.5 { 1.0 - q } else { q };
    let t_x = (-2.0 * x.ln()).sqrt();
    let zeta_x = ((c_2 * t_x + c_1) * t_x + c_0) / (((d_3 * t_x + d_2) * t_x + d_1) * t_x + 1.0);
    let q_q = t_x - zeta_x;

    if q > 0.5 { -q_q } else { q_q }
}

/// Combines two distributions given by their medians (`a_m`, `b_m`) and their values at `p` (`a_p`, `b_p`)
/// (`CombineDistributions`).
pub fn combine_distributions(a_m: f64, a_p: f64, b_m: f64, b_p: f64, p: f64) -> f64 {
    let c_m = a_m + b_m;
    let y_1 = a_p - a_m;
    let y_2 = b_p - b_m;
    let y_3 = (y_1.powf(2.0) + y_2.powf(2.0)).sqrt();

    if p < 50.0 { c_m + y_3 } else { c_m - y_3 }
}

/// Long-term variability `Y_e` at time percentage `p`, dB, and `A_Y` (`LongTermVariability`). `a_t` is the loss
/// relative to free space the variability is applied to, dB.
#[allow(clippy::too_many_arguments)]
pub fn long_term_variability(
    d_r1_km: f64,
    d_r2_km: f64,
    d_km: f64,
    f_mhz: f64,
    p: f64,
    f_theta_h: f64,
    a_t: f64,
) -> (f64, f64) {
    let d_qs_km = 65.0 * (100.0 / f_mhz).powf(THIRD); // [Eqn 14-1]
    let d_lq_km = d_r1_km + d_r2_km; // [Eqn 14-2]
    let d_q_km = d_lq_km + d_qs_km; // [Eqn 14-3]

    // [Eqn 14-4]
    let d_e_km = if d_km <= d_q_km {
        (130.0 * d_km) / d_q_km
    } else {
        130.0 + d_km - d_q_km
    };

    // [Eqns 14-5 and 14-6]
    let (g_10, g_90) = if f_mhz > 1600.0 {
        (1.05, 1.05)
    } else {
        (
            (0.21 * (5.22 * (f_mhz / 200.0).log10()).sin()) + 1.28,
            (0.18 * (5.22 * (f_mhz / 200.0).log10()).sin()) + 1.23,
        )
    };

    // Constants from Tech Note 101, Vol 2. Column 1: Table III.4, Row A* (page III-50). Column 2: Table III.3,
    // Row A* (page III-49). Column 3: Table III.5, Row Continental Temperate (page III-51).
    let c_1 = [2.93e-4, 5.25e-4, 1.59e-5];
    let c_2 = [3.78e-8, 1.57e-6, 1.56e-11];
    let c_3 = [1.02e-7, 4.70e-7, 2.77e-8];
    let n_1 = [2.00, 1.97, 2.32];
    let n_2 = [2.88, 2.31, 4.08];
    let n_3 = [3.15, 2.90, 3.25];
    let f_inf = [3.2, 5.4, 0.0];
    let f_m = [8.2, 10.0, 3.9];

    // [Y_0(90), Y_0(10), V(50)]
    let z_db: [f64; 3] = std::array::from_fn(|i| {
        let f_2 = f_inf[i] + ((f_m[i] - f_inf[i]) * (-c_2[i] * d_e_km.powf(n_2[i])).exp());
        (c_1[i] * d_e_km.powf(n_1[i]) - f_2) * (-c_3[i] * d_e_km.powf(n_3[i])).exp() + f_2
    });

    let y_p_db = if p == 50.0 {
        z_db[2]
    } else if p > 50.0 {
        let z_90 = inverse_ccdf(90.0 / 100.0);
        let z_p = inverse_ccdf(p / 100.0);
        let c_p = z_p / z_90;
        let y = c_p * (-z_db[0] * g_90);
        y + z_db[2]
    } else {
        let c_p = if p >= 10.0 {
            let z_10 = inverse_ccdf(10.0 / 100.0);
            let z_p = inverse_ccdf(p / 100.0);
            z_p / z_10
        } else {
            // Values for p < 10: [15], Table 10, page 34, Climate 6.
            let ps = [1.0, 2.0, 5.0, 10.0];
            let c_ps = [1.9507, 1.7166, 1.3265, 1.0000];
            let dist = upper_bound(&P, p);
            linear_interpolation(ps[dist - 1], c_ps[dist - 1], ps[dist], c_ps[dist], p)
        };
        let y = c_p * (z_db[1] * g_10);
        y + z_db[2]
    };

    let y_10_db = (z_db[1] * g_10) + z_db[2]; // [Eqn 14-20]
    let y_ei_db = f_theta_h * y_p_db; // [Eqn 14-21]
    let y_ei_10_db = f_theta_h * y_10_db; // [Eqn 14-22]

    // A_Y "is used to prevent available signal powers from exceeding levels expected for free-space propagation by an
    // unrealistic amount when the variability about L_b(50) is large and L_b(50) is near its free-space level"
    // [ES-83-3, p3-4].
    let a_yi = (a_t + y_ei_10_db) - 3.0; // [Eqn 14-23]
    let a_y = max(a_yi, 0.0); // [Eqn 14-24]
    let mut y_e_db = y_ei_db - a_y; // [Eqn 14-25]

    // For percentages below 10 %, "prevent available signal powers from exceeding levels expected from free-space
    // levels by unrealistic amounts" [Gierhart 1970].
    if p < 10.0 {
        let c_y = [-5.0, -4.5, -3.7, 0.0];
        let dist = upper_bound(&P, p);
        let c_yi = linear_interpolation(P[dist - 1], c_y[dist - 1], P[dist], c_y[dist], p);

        y_e_db += a_t;
        if y_e_db > -c_yi {
            y_e_db = -c_yi;
        }
        y_e_db -= a_t;
    }

    (y_e_db, a_y)
}
