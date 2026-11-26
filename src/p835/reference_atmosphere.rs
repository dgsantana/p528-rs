//! The mean annual global reference atmosphere (NTIA `MeanAnnualGlobalReferenceAtmosphere.cpp`, `Conversions.cpp`).
//!
//! NTIA's functions return their error codes (`-1` below a profile, `-2` above it) as if they were temperatures and
//! pressures, and the computation carries on with them. Here a height outside a profile is an error instead.

/// A height outside the reference atmosphere's profiles (below 0 km or above 100 km, or in the gap NTIA's profiles
/// leave just below 86 km).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OutsideAtmosphere;

/// Result of a profile lookup.
pub type Profile = Result<f64, OutsideAtmosphere>;

fn geopotential_height(h_km: f64) -> f64 {
    (6356.766 * h_km) / (6356.766 + h_km)
}

/// NTIA's range check, common to every profile: 0 km to 100 km.
fn in_range(h_km: f64) -> Result<(), OutsideAtmosphere> {
    if (0.0..=100.0).contains(&h_km) {
        Ok(())
    } else {
        Err(OutsideAtmosphere)
    }
}

/// Temperature, kelvin, at geometric height `h_km` (`GlobalTemperature`).
pub fn global_temperature(h_km: f64) -> Profile {
    in_range(h_km)?;
    if h_km < 86.0 {
        temperature_regime_1(geopotential_height(h_km))
    } else {
        temperature_regime_2(h_km)
    }
}

fn temperature_regime_1(h_prime_km: f64) -> Profile {
    Ok(if h_prime_km < 0.0 {
        return Err(OutsideAtmosphere);
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
        return Err(OutsideAtmosphere);
    })
}

fn temperature_regime_2(h_km: f64) -> Profile {
    if h_km <= 91.0 {
        Ok(186.8673)
    } else {
        Ok(263.1905 - 76.3232 * (1.0 - ((h_km - 91.0) / 19.9429).powf(2.0)).sqrt())
    }
}

/// Dry pressure, hPa, at geometric height `h_km` (`GlobalPressure`).
pub fn global_pressure(h_km: f64) -> Profile {
    in_range(h_km)?;
    if h_km < 86.0 {
        pressure_regime_1(geopotential_height(h_km))
    } else {
        Ok(pressure_regime_2(h_km))
    }
}

fn pressure_regime_1(h_prime_km: f64) -> Profile {
    Ok(if h_prime_km < 0.0 {
        return Err(OutsideAtmosphere);
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
        return Err(OutsideAtmosphere);
    })
}

fn pressure_regime_2(h_km: f64) -> f64 {
    let a_0 = 95.571899;
    let a_1 = -4.011801;
    let a_2 = 6.424731e-2;
    let a_3 = -4.789660e-4;
    let a_4 = 1.340543e-6;
    (a_0 + a_1 * h_km + a_2 * h_km.powf(2.0) + a_3 * h_km.powf(3.0) + a_4 * h_km.powf(4.0)).exp()
}

/// Water vapour density, g/m³, at geometric height `h_km` for a surface density `rho_0` (`GlobalWaterVapourDensity`).
pub fn global_water_vapour_density(h_km: f64, rho_0: f64) -> Profile {
    in_range(h_km)?;
    let h_0_km = 2.0; // scale height
    Ok(rho_0 * (-h_km / h_0_km).exp())
}

/// Water vapour partial pressure, hPa, from density `rho` (g/m³) and temperature (`WaterVapourDensityToPressure`).
pub fn water_vapour_density_to_pressure(rho: f64, t_kelvin: f64) -> f64 {
    (rho * t_kelvin) / 216.7
}
