//! Comparison of one case between the port and NTIA's reference.

use ntia_p528_reference::Reference;
use p528::{P528Error, P528Input, P528Result, Polarization, p528};

/// Agreement asked of every loss, dB.
const TOLERANCE_DB: f64 = 0.01;
/// Agreement asked of the distance the result is for, km (1 cm).
const TOLERANCE_KM: f64 = 1e-5;
/// Agreement asked of the elevation angle, rad.
const TOLERANCE_RAD: f64 = 1e-9;

/// The input as NTIA's call takes it.
pub fn input(d_km: f64, h1_m: f64, h2_m: f64, f_mhz: f64, pol: i32, p: f64) -> P528Input {
    P528Input {
        distance_km: d_km,
        h1_m,
        h2_m,
        frequency_mhz: f_mhz,
        time_percent: p,
        polarization: if pol == 0 {
            Polarization::Horizontal
        } else {
            Polarization::Vertical
        },
    }
}

/// Runs both on one case; `None` when they agree, else the reason.
pub fn check(d_km: f64, h1_m: f64, h2_m: f64, f_mhz: f64, pol: i32, p: f64) -> Option<String> {
    let reference = ntia_p528_reference::p528(d_km, h1_m, h2_m, f_mhz, pol, p);
    let port = p528(&input(d_km, h1_m, h2_m, f_mhz, pol, p));
    let label = format!("d={d_km} km h1={h1_m} m h2={h2_m} m f={f_mhz} MHz pol={pol} p={p}");
    compare(&label, port, &reference)
}

/// What the port must return for a case, from NTIA's answer.
enum Expected {
    /// A loss within tolerance of NTIA's.
    Loss,
    /// A refusal with this NTIA code.
    Refusal(i32),
    /// `P528Error::OutsideAtmosphere`, where NTIA returns NaN.
    OutsideAtmosphere,
}

/// Compares one result of the port with NTIA's for the same case; `None` when they agree.
///
/// Where NTIA returns success with a value that is not a loss, the port must refuse instead: NaN (the path leaves the
/// reference atmosphere) with [`P528Error::OutsideAtmosphere`], and 0 dB with no mode (coincident terminals) with
/// [`P528Error::CoincidentTerminals`].
pub fn compare(
    label: &str,
    port: Result<P528Result, P528Error>,
    reference: &Reference,
) -> Option<String> {
    let r = &reference.output;
    let expected = if reference.code != 0 && reference.code != 11 {
        Expected::Refusal(reference.code)
    } else if r.propagation_mode == 0 {
        Expected::Refusal(10)
    } else if r.a_db.is_nan() {
        Expected::OutsideAtmosphere
    } else {
        Expected::Loss
    };
    match (port, expected) {
        (Ok(port), Expected::Loss) => {
            let close = (port.basic_transmission_loss_db - r.a_db).abs() <= TOLERANCE_DB
                && (port.free_space_loss_db - r.a_fs_db).abs() <= TOLERANCE_DB
                && (port.absorption_loss_db - r.a_a_db).abs() <= TOLERANCE_DB
                && (port.distance_km - r.d_km).abs() <= TOLERANCE_KM
                && (port.theta_h1_rad - r.theta_h1_rad).abs() <= TOLERANCE_RAD;
            let same = port.mode.ntia_code() == r.propagation_mode
                && port.warnings.bits() as i32 == r.warnings;
            (!close || !same).then(|| {
                format!(
                    "{label}: port A {:.4} fs {:.4} a {:.4} d {:.6} th {:.3e} mode {} w{:#x};                      NTIA A {:.4} fs {:.4} a {:.4} d {:.6} th {:.3e} mode {} w{:#x}",
                    port.basic_transmission_loss_db,
                    port.free_space_loss_db,
                    port.absorption_loss_db,
                    port.distance_km,
                    port.theta_h1_rad,
                    port.mode.ntia_code(),
                    port.warnings.bits(),
                    r.a_db,
                    r.a_fs_db,
                    r.a_a_db,
                    r.d_km,
                    r.theta_h1_rad,
                    r.propagation_mode,
                    r.warnings
                )
            })
        }
        (Err(P528Error::OutsideAtmosphere), Expected::OutsideAtmosphere) => None,
        (Err(error), Expected::Refusal(code)) => (error.ntia_code() != Some(code)).then(|| {
            format!("{label}: port refuses with {error:?}, NTIA's equivalent code is {code}")
        }),
        (port, _) => Some(format!(
            "{label}: port returns {port:?}, NTIA returns code {} with {r:?}",
            reference.code
        )),
    }
}

/// Fails with the disagreements listed.
pub fn assert_all_agree(failures: &[String], cases: usize) {
    assert!(
        failures.is_empty(),
        "{} of {cases} cases disagree:\n{}",
        failures.len(),
        failures
            .iter()
            .take(40)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
