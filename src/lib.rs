//! Recommendation ITU-R P.528: basic transmission loss for aeronautical paths, from a ground or airborne terminal to an
//! airborne one, over a smooth earth, 100 MHz to 30 GHz.
//!
//! A port of NTIA/ITS's reference implementation of P.528-5 Annex 2 (the step-by-step method), including the P.676
//! gaseous absorption and the P.835 reference atmosphere that reference uses. Its results are checked against that
//! reference, compiled from NTIA's C++ sources, in this crate's tests.
//!
//! ```
//! use p528::{P528Input, Polarization, p528};
//!
//! let result = p528(&P528Input {
//!     distance_km: 15.0,
//!     h1_m: 10.0,
//!     h2_m: 1000.0,
//!     frequency_mhz: 500.0,
//!     time_percent: 50.0,
//!     polarization: Polarization::Horizontal,
//! })?;
//! assert!((result.basic_transmission_loss_db - 110.0).abs() < 0.05);
//! # Ok::<(), p528::P528Error>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod error;
mod input;
mod model;
mod output;
mod p676;
mod p835;

pub use error::P528Error;
pub use input::{P528Input, Polarization};
pub use model::{P528Model, p528};
pub use output::{P528Result, PropagationMode, Warnings};
