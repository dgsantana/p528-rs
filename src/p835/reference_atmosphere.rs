//! The mean annual global reference atmosphere (NTIA `MeanAnnualGlobalReferenceAtmosphere.cpp`, `Conversions.cpp`).
//!
//! As in NTIA's code, a height outside a profile's range returns the error code as the value (`-1` below, `-2` above).
//! P.528 limits terminals to 80 km, inside every profile, so these values are never reached through the public API.

/// NTIA's `ERROR_HEIGHT_TOO_SMALL`, returned as a value.
const HEIGHT_TOO_SMALL: f64 = -1.0;
/// NTIA's `ERROR_HEIGHT_TOO_LARGE`, returned as a value.
const HEIGHT_TOO_LARGE: f64 = -2.0;

fn geopotential_height(h_km: f64) -> f64 {
    (6356.766 * h_km) / (6356.766 + h_km)
}

/// Temperature, kelvin, at geometric height `h_km` (`GlobalTemperature`).
pub fn global_temperature(h_km: f64) -> f64 {
    if h_km < 0.0 {
        return HEIGHT_TOO_SMALL;
    }
    if h_km > 100.0 {
        return HEIGHT_TOO_LARGE;
    }
    if h_km < 86.0 {
        temperature_regime_1(geopotential_height(h_km))
    } else {
        temperature_regime_2(h_km)
    }
}

fn temperature_regime_1(h_prime_km: f64) -> f64 {
    if h_prime_km < 0.0 {
        HEIGHT_TOO_SMALL
    } else if h_prime_km <= 11.0 {
        288.15 - 6.5 * h_prime_km
    } else if h_prime_km <= 20.0 {
        216.65
    } else if h_prime_km <= 32.0 {
        216.65 + (h_prime_km - 20.0)
    } else if h_prime_km <= 47.0 {
        228.65 + 2.8 * (h_prime_km - 32.0)
    } else if h_prime_km <= 51.0 {
        270.65
    } else if h_prime_km <= 71.0 {
        270.65 - 2.8 * (h_prime_km - 51.0)
    } else if h_prime_km <= 84.852 {
        214.65 - 2.0 * (h_prime_km - 71.0)
    } else {
        HEIGHT_TOO_LARGE
    }
}

fn temperature_regime_2(h_km: f64) -> f64 {
    if h_km < 86.0 {
        HEIGHT_TOO_SMALL
    } else if h_km <= 91.0 {
        186.8673
    } else if h_km <= 100.0 {
        263.1905 - 76.3232 * (1.0 - ((h_km - 91.0) / 19.9429).powf(2.0)).sqrt()
    } else {
        HEIGHT_TOO_LARGE
    }
}

/// Dry pressure, hPa, at geometric height `h_km` (`GlobalPressure`).
pub fn global_pressure(h_km: f64) -> f64 {
    if h_km < 0.0 {
        return HEIGHT_TOO_SMALL;
    }
    if h_km > 100.0 {
        return HEIGHT_TOO_LARGE;
    }
    if h_km < 86.0 {
        pressure_regime_1(geopotential_height(h_km))
    } else {
        pressure_regime_2(h_km)
    }
}

fn pressure_regime_1(h_prime_km: f64) -> f64 {
    if h_prime_km < 0.0 {
        HEIGHT_TOO_SMALL
    } else if h_prime_km <= 11.0 {
        1013.25 * (288.15 / (288.15 - 6.5 * h_prime_km)).powf(-34.1632 / 6.5)
    } else if h_prime_km <= 20.0 {
        226.3226 * (-34.1632 * (h_prime_km - 11.0) / 216.65).exp()
    } else if h_prime_km <= 32.0 {
        54.74980 * (216.65 / (216.65 + (h_prime_km - 20.0))).powf(34.1632)
    } else if h_prime_km <= 47.0 {
        8.680422 * (228.65 / (228.65 + 2.8 * (h_prime_km - 32.0))).powf(34.1632 / 2.8)
    } else if h_prime_km <= 51.0 {
        1.109106 * (-34.1632 * (h_prime_km - 47.0) / 270.65).exp()
    } else if h_prime_km <= 71.0 {
        0.6694167 * (270.65 / (270.65 - 2.8 * (h_prime_km - 51.0))).powf(-34.1632 / 2.8)
    } else if h_prime_km <= 84.852 {
        0.03956649 * (214.65 / (214.65 - 2.0 * (h_prime_km - 71.0))).powf(-34.1632 / 2.0)
    } else {
        HEIGHT_TOO_LARGE
    }
}

fn pressure_regime_2(h_km: f64) -> f64 {
    if h_km < 86.0 {
        return HEIGHT_TOO_SMALL;
    }
    if h_km > 100.0 {
        return HEIGHT_TOO_LARGE;
    }
    let a_0 = 95.571899;
    let a_1 = -4.011801;
    let a_2 = 6.424731e-2;
    let a_3 = -4.789660e-4;
    let a_4 = 1.340543e-6;
    (a_0 + a_1 * h_km + a_2 * h_km.powf(2.0) + a_3 * h_km.powf(3.0) + a_4 * h_km.powf(4.0)).exp()
}

/// Water vapour density, g/mÂ³, at geometric height `h_km` for a surface density `rho_0` (`GlobalWaterVapourDensity`).
pub fn global_water_vapour_density(h_km: f64, rho_0: f64) -> f64 {
    if h_km < 0.0 {
        return HEIGHT_TOO_SMALL;
    }
    if h_km > 100.0 {
        return HEIGHT_TOO_LARGE;
    }
    let h_0_km = 2.0; // scale height
    rho_0 * (-h_km / h_0_km).exp()
}

/// Water vapour partial pressure, hPa, from density `rho` (g/m³) and temperature (`WaterVapourDensityToPressure`).
pub fn water_vapour_density_to_pressure(rho: f64, t_kelvin: f64) -> f64 {
    (rho * t_kelvin) / 216.7
}
