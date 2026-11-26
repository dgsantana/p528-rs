//! Recommendation ITU-R P.835, the mean annual global reference atmosphere: the temperature, pressure and water vapour
//! profiles the gaseous-absorption ray trace reads. Only what NTIA's P.528 reference calls is ported.

mod reference_atmosphere;

pub use reference_atmosphere::{
    OutsideAtmosphere, Profile, global_pressure, global_temperature, global_water_vapour_density,
    water_vapour_density_to_pressure,
};
