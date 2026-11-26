//! Smooth-earth diffraction (NTIA `SmoothEarthDiffraction.cpp`, P.528 Annex 2, section 10).
//!
//! References: [Vogler 1964] *Calculation of Groundwave Attenuation in the Far Diffraction Region*; [FAA-ES-83-3]
//! *The IF-77 Electromagnetic Wave Propagation Model*, Gierhart and Johnson, 1983.

use super::{EPSILON_R, SIGMA, THIRD};
use crate::Polarization;

/// Distance function, dB ([Vogler 1964, equation 13]).
fn distance_function(x_km: f64) -> f64 {
    0.05751 * x_km - 10.0 * x_km.log10()
}

/// Height function, dB.
fn height_function(x_km: f64, k: f64) -> f64 {
    let y_db = 40.0 * x_km.log10() - 117.0; // [FAA-ES-83-3, equation 73]
    let g_x_db = distance_function(x_km);

    if x_km <= 200.0 {
        let x_t_km = 450.0 / -k.log10().powf(3.0); // [Eqn 109]

        // [Eqn 110]
        if x_km >= x_t_km {
            if y_db.abs() < 117.0 { y_db } else { -117.0 }
        } else {
            20.0 * k.log10() - 15.0 + (0.000025 * x_km.powf(2.0) / k)
        }
    } else if x_km > 2000.0 {
        // [Vogler 1964] F_x ~= G_x for large x (Figure 7).
        g_x_db
    } else {
        // Blend y with G_x for 200 < x < 2000: [FAA-ES-83-3, equations 72 and 75].
        let w = 0.0134 * x_km * (-0.005 * x_km).exp();
        w * y_db + (1.0 - w) * g_x_db
    }
}

/// Smooth-earth diffraction loss, dB, at path distance `d_0_km` for terminals with horizon distances `d_1_km` and
/// `d_2_km` (`SmoothEarthDiffraction`).
pub fn smooth_earth_diffraction(
    d_1_km: f64,
    d_2_km: f64,
    f_mhz: f64,
    d_0_km: f64,
    polarization: Polarization,
) -> f64 {
    let s = 18000.0 * SIGMA / f_mhz;

    let k = match polarization {
        Polarization::Horizontal => {
            0.01778 * f_mhz.powf(-THIRD) * ((EPSILON_R - 1.0).powf(2.0) + s.powf(2.0)).powf(-0.25)
        }
        Polarization::Vertical => {
            0.01778
                * f_mhz.powf(-THIRD)
                * ((EPSILON_R.powf(2.0) + s.powf(2.0))
                    / ((EPSILON_R - 1.0).powf(2.0) + s.powf(2.0)).powf(0.5))
                .powf(0.5)
        }
    };

    let b_0 = 1.607;

    // [Vogler 1964, equation 2] with C_0 = 1, the 4/3 earth.
    let x_0_km = (b_0 - k) * f_mhz.powf(THIRD) * d_0_km;
    let x_1_km = (b_0 - k) * f_mhz.powf(THIRD) * d_1_km;
    let x_2_km = (b_0 - k) * f_mhz.powf(THIRD) * d_2_km;

    let g_x_db = distance_function(x_0_km);
    let f_x1_db = height_function(x_1_km, k);
    let f_x2_db = height_function(x_2_km, k);

    // [Vogler 1964, equation 1] with C_1(K, b^0) = 20, the approximate value for all K (Figure 5).
    g_x_db - f_x1_db - f_x2_db - 20.0
}
