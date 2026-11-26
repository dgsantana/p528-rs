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

/// Within `tolerance`, or NaN in both: NTIA's reference returns NaN on some accepted inputs (from about 2,400 km,
/// where the troposcatter common volume rises above the 100 km the reference atmosphere covers), and the port must
/// return the same.
fn near(port: f64, reference: f64, tolerance: f64) -> bool {
    (port - reference).abs() <= tolerance || (port.is_nan() && reference.is_nan())
}

fn compare(
    label: &str,
    port: Result<P528Result, P528Error>,
    reference: &Reference,
) -> Option<String> {
    let reference_ok = reference.code == 0 || reference.code == 11;
    let r = &reference.output;
    match port {
        Ok(port) if reference_ok => {
            let close = near(port.basic_transmission_loss_db, r.a_db, TOLERANCE_DB)
                && near(port.free_space_loss_db, r.a_fs_db, TOLERANCE_DB)
                && near(port.absorption_loss_db, r.a_a_db, TOLERANCE_DB)
                && near(port.distance_km, r.d_km, TOLERANCE_KM)
                && near(port.theta_h1_rad, r.theta_h1_rad, TOLERANCE_RAD);
            let same = port.mode.ntia_code() == r.propagation_mode
                && port.warnings.bits() as i32 == r.warnings;
            (!close || !same).then(|| {
                format!(
                    "{label}: port A {:.4} fs {:.4} a {:.4} d {:.6} th {:.3e} mode {} w{:#x}; \
                     NTIA A {:.4} fs {:.4} a {:.4} d {:.6} th {:.3e} mode {} w{:#x}",
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
        Err(error) if !reference_ok => (error.ntia_code() != Some(reference.code)).then(|| {
            format!(
                "{label}: port refuses with {error:?}, NTIA with code {}",
                reference.code
            )
        }),
        Ok(port) => Some(format!(
            "{label}: port returns {port:?}, NTIA refuses with code {}",
            reference.code
        )),
        Err(error) => Some(format!(
            "{label}: port refuses with {error:?}, NTIA returns {r:?}"
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
