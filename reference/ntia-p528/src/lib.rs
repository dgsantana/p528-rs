//! Bindings to NTIA's C++ P.528, the reference the library's tests compare with. Test-only.
//!
//! The structures mirror `ntia/include/p528.h` field for field, so `P528_Ex` can fill them.

use std::ffi::c_int;

/// NTIA's `Result`.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub struct Output {
    /// Mode of propagation: 0 not set, 1 line of sight, 2 diffraction, 3 troposcatter.
    pub propagation_mode: c_int,
    /// Warning bits.
    pub warnings: c_int,
    /// Path distance used in the calculation, km.
    pub d_km: f64,
    /// Basic transmission loss, dB.
    pub a_db: f64,
    /// Free-space basic transmission loss, dB.
    pub a_fs_db: f64,
    /// Median atmospheric absorption loss, dB.
    pub a_a_db: f64,
    /// Elevation angle of the ray at the low terminal, rad.
    pub theta_h1_rad: f64,
}

/// NTIA's `Terminal`.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub struct Terminal {
    /// Real terminal height, km.
    pub h_r_km: f64,
    /// Effective terminal height, km.
    pub h_e_km: f64,
    /// `delta_h`, km.
    pub delta_h_km: f64,
    /// Ray-traced horizon distance, km.
    pub d_r_km: f64,
    /// Ray path length to the horizon, km.
    pub a_km: f64,
    /// Central angle to the smooth-earth horizon, rad.
    pub phi_rad: f64,
    /// Incident angle of the grazing ray at the terminal, rad.
    pub theta_rad: f64,
    /// Median atmospheric absorption loss, dB.
    pub a_a_db: f64,
}

/// NTIA's `Path`.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub struct PathParams {
    /// Maximum line-of-sight distance, km.
    pub d_ml_km: f64,
    /// `d_0`, km.
    pub d_0_km: f64,
    /// Distance at which smooth-earth diffraction is 0 dB, km.
    pub d_d_km: f64,
}

/// NTIA's `TroposcatterParams`.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub struct Troposcatter {
    /// Scattering distance, km.
    pub d_s_km: f64,
    /// Half the scattering distance, km.
    pub d_z_km: f64,
    /// Height of the common-volume cross-over point, km.
    pub h_v_km: f64,
    /// Scattering angle, rad.
    pub theta_s: f64,
    /// Cross-over angle, rad.
    pub theta_a: f64,
    /// Troposcatter loss, dB.
    pub a_s_db: f64,
    /// Troposcatter loss of the previous search point, dB.
    pub a_s_prev_db: f64,
    /// Troposcatter line slope.
    pub m_s: f64,
}

/// NTIA's `LineOfSightParams`.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub struct LineOfSight {
    /// `z`, km.
    pub z_km: [f64; 2],
    /// Path distance, km.
    pub d_km: f64,
    /// Direct ray length, km.
    pub r_0_km: f64,
    /// Indirect ray length, km.
    pub r_12_km: f64,
    /// `D`, km.
    pub big_d_km: [f64; 2],
    /// Take-off angle at the low terminal, rad.
    pub theta_h1_rad: f64,
    /// Take-off angle at the high terminal, rad.
    pub theta_h2_rad: f64,
    /// `theta`, rad.
    pub theta: [f64; 2],
    /// Adjusted earth radius, km.
    pub a_a_km: f64,
    /// Ray length path difference, km.
    pub delta_r_km: f64,
    /// Line-of-sight loss, dB.
    pub a_los_db: f64,
}

/// Everything one reference call returns.
#[derive(Debug, Default, Clone, Copy)]
pub struct Reference {
    /// NTIA's return code: 0 success, 11 success with warnings, else an error.
    pub code: i32,
    /// NTIA's `Result`.
    pub output: Output,
    /// Low terminal.
    pub terminal_1: Terminal,
    /// High terminal.
    pub terminal_2: Terminal,
    /// Troposcatter parameters (transhorizon paths only).
    pub tropo: Troposcatter,
    /// Path parameters.
    pub path: PathParams,
    /// Line-of-sight parameters.
    pub los: LineOfSight,
}

unsafe extern "C" {
    #[link_name = "P528_Ex"]
    #[allow(clippy::too_many_arguments)]
    fn p528_ex(
        d_km: f64,
        h_1_meter: f64,
        h_2_meter: f64,
        f_mhz: f64,
        t_pol: c_int,
        p: f64,
        result: *mut Output,
        terminal_1: *mut Terminal,
        terminal_2: *mut Terminal,
        tropo: *mut Troposcatter,
        path: *mut PathParams,
        los_params: *mut LineOfSight,
    ) -> c_int;
}

/// NTIA's `P528_Ex`. `t_pol` is 0 horizontal, 1 vertical; `p` is the time percentage.
pub fn p528(d_km: f64, h_1_meter: f64, h_2_meter: f64, f_mhz: f64, t_pol: i32, p: f64) -> Reference {
    let mut reference = Reference::default();
    // SAFETY: plain values in; every out pointer is a valid, exclusively borrowed `#[repr(C)]` local whose layout
    // matches NTIA's structure.
    reference.code = unsafe {
        p528_ex(
            d_km,
            h_1_meter,
            h_2_meter,
            f_mhz,
            t_pol,
            p,
            &mut reference.output,
            &mut reference.terminal_1,
            &mut reference.terminal_2,
            &mut reference.tropo,
            &mut reference.path,
            &mut reference.los,
        )
    };
    reference
}
