//! Recommendation ITU-R P.676 gaseous absorption along a slant path, ray-traced through the P.835 reference
//! atmosphere, as NTIA's P.528 reference computes it. Only what P.528 calls is ported.

mod ray_trace;
mod slant_path;
mod specific_attenuation;
mod spectroscopy;

pub use ray_trace::SlantPathLayers;
pub use slant_path::slant_path_attenuation;
