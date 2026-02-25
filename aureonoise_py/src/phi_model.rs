//! Phi Head Geometry + Pinna Multipath
//!
//! Ellipsoid head model based on the golden ratio (phi), computing:
//! - ITD via Woodworth model with asymmetric head shape
//! - Pinna multipath (5-tap FIR with phi-scaled delays)
//! - Torso response (2-tap with LP filtering)
//! - Distance response (gain + lowpass rolloff)
//!
//! Ported from aureonoise::spatial::phi (C++ Max external).

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::constants::{PHI, INV_PHI, PI, TWO_PI};
use crate::math::{clamp, clamp01};

// ---------------------------------------------------------------------------
// Woodworth ITD with explicit radius
// ---------------------------------------------------------------------------

/// Woodworth ITD model: returns time-of-arrival difference (seconds)
/// for a sphere of given radius at angle theta.
#[inline]
fn woodworth_itd_radius(radius: f64, theta: f64) -> f64 {
    const C: f64 = 343.0; // speed of sound m/s
    let limit = 0.5 * PI;
    let th = clamp(theta, -limit, limit);
    (radius / C) * (th + th.sin())
}

// ---------------------------------------------------------------------------
// Geometry structs
// ---------------------------------------------------------------------------

/// Configuration for the ellipsoid head geometry.
#[derive(Clone, Debug)]
pub struct GeometryConfig {
    pub head_b: f64,       // lateral semi-axis (m)
    pub phi_ratio: f64,    // golden ratio
    pub nose_ratio: f64,   // nose extrusion relative to lateral axis
    pub front_base: f64,   // base ellipsoid scale on front axis
    pub back_scale: f64,   // back semi-axis relative to lateral
}

impl Default for GeometryConfig {
    fn default() -> Self {
        Self {
            head_b: 0.0875,
            phi_ratio: PHI,
            nose_ratio: 0.20,
            front_base: 0.83,
            back_scale: 1.05,
        }
    }
}

/// Computed ellipsoid geometry (axes in metres).
#[derive(Clone, Debug)]
pub struct Geometry {
    pub axis_lat: f64,
    pub axis_vert: f64,
    pub axis_front: f64,
    pub axis_front_base: f64,
    pub axis_back: f64,
    pub nose_len: f64,
}

impl Default for Geometry {
    fn default() -> Self {
        Self {
            axis_lat: 0.0875,
            axis_vert: 0.1416,
            axis_front: 0.0,
            axis_front_base: 0.0,
            axis_back: 0.0,
            nose_len: 0.0,
        }
    }
}

/// Head computation result: radii, lateral component, ITD.
#[derive(Clone, Debug)]
pub struct HeadResult {
    pub horizontal_radius: f64,
    pub direction_radius: f64,
    pub lateral: f64,
    pub itd_samples: f64,
}

impl Default for HeadResult {
    fn default() -> Self {
        Self {
            horizontal_radius: 0.0,
            direction_radius: 0.0,
            lateral: 0.0,
            itd_samples: 0.0,
        }
    }
}

/// Pinna tuning parameters.
#[derive(Clone, Debug)]
pub struct PinnaTuning {
    pub elev_norm: f64,
    pub elev_notch: f64,
    pub notch_shift: bool,
    pub elev_plus: bool,
}

impl Default for PinnaTuning {
    fn default() -> Self {
        Self {
            elev_norm: 0.0,
            elev_notch: 0.0,
            notch_shift: false,
            elev_plus: false,
        }
    }
}

/// Pinna response: amplitude, HF strength, peak/notch freqs, 5-tap FIR.
#[derive(Clone, Debug)]
pub struct PinnaResponse {
    pub amp_scale: f64,
    pub hf_strength: f64,
    pub peak_freqs: [f64; 3],
    pub notch_freqs: [f64; 3],
    pub tap_count: usize,
    pub tap_delay_sec: [f64; 5],
    pub tap_gains: [f64; 5],
}

impl Default for PinnaResponse {
    fn default() -> Self {
        Self {
            amp_scale: 1.0,
            hf_strength: 0.0,
            peak_freqs: [0.0; 3],
            notch_freqs: [0.0; 3],
            tap_count: 0,
            tap_delay_sec: [0.0; 5],
            tap_gains: [0.0; 5],
        }
    }
}

/// Torso tuning parameters.
#[derive(Clone, Debug)]
pub struct TorsoTuning {
    pub torso_mix: f64,
    pub torso_ms: f64,
    pub torso_hp_hz: f64,
}

impl Default for TorsoTuning {
    fn default() -> Self {
        Self {
            torso_mix: 0.0,
            torso_ms: 1.0,
            torso_hp_hz: 600.0,
        }
    }
}

/// Torso response: 2-tap delay with LP filtering.
#[derive(Clone, Debug)]
pub struct TorsoResponse {
    pub tap_count: usize,
    pub delay_samples: [f64; 2],
    pub gains: [f64; 2],
    pub lp_alpha: [f64; 2],
    pub pan_abs: f64,
    pub active: bool,
}

impl Default for TorsoResponse {
    fn default() -> Self {
        Self {
            tap_count: 0,
            delay_samples: [0.0; 2],
            gains: [0.0; 2],
            lp_alpha: [0.0; 2],
            pan_abs: 0.0,
            active: false,
        }
    }
}

/// Distance response: gain attenuation + lowpass rolloff.
#[derive(Clone, Debug)]
pub struct DistanceResponse {
    pub direct_gain: f64,
    pub hf_weight: f64,
    pub lowpass_alpha: f64,
}

impl Default for DistanceResponse {
    fn default() -> Self {
        Self {
            direct_gain: 1.0,
            hf_weight: 1.0,
            lowpass_alpha: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Pure functions (no state, composable)
// ---------------------------------------------------------------------------

/// Build geometry from configuration.
pub fn build_geometry(cfg: &GeometryConfig) -> Geometry {
    let axis_lat = cfg.head_b;
    let axis_vert = cfg.phi_ratio * cfg.head_b;
    let axis_front_base = cfg.front_base * cfg.head_b;
    let nose_len = cfg.nose_ratio * cfg.head_b;
    let axis_front = axis_front_base + nose_len;
    let axis_back = cfg.back_scale * cfg.head_b;

    Geometry {
        axis_lat,
        axis_vert,
        axis_front,
        axis_front_base,
        axis_back,
        nose_len,
    }
}

/// Compute the radius of the ellipsoid in a given direction (unit vector).
/// Front/back asymmetry: uz >= 0 uses axis_front, uz < 0 uses axis_back.
pub fn ellipsoid_radius(geom: &Geometry, ux: f64, uy: f64, uz: f64) -> f64 {
    let eps = 1.0e-12;
    let al = geom.axis_lat.max(eps);
    let av = geom.axis_vert.max(eps);
    let af = geom.axis_front.max(eps);
    let ab = geom.axis_back.max(eps);

    let axis_z = if uz >= 0.0 { af } else { ab };

    let denom = (ux * ux) / (al * al)
        + (uy * uy) / (av * av)
        + (uz * uz) / (axis_z * axis_z);

    if denom <= eps {
        return al;
    }
    1.0 / denom.sqrt()
}

/// Compute head result: horizontal/direction radii, lateral component, ITD in samples.
pub fn compute_head_result(
    geom: &Geometry,
    sample_rate: f64,
    pan: f64,
    pitch_rad: f64,
    phi_distance: f64,
) -> HeadResult {
    let mut res = HeadResult::default();

    let theta = clamp(pan, -1.0, 1.0) * (PI * 0.5);
    let dir_x = theta.sin();
    let dir_y = 0.0_f64;
    let dir_z = theta.cos();

    // Pitch rotation (around X axis, negative pitch_rad)
    let c = (-pitch_rad).cos();
    let s = (-pitch_rad).sin();
    let ux = dir_x;
    let uy = c * dir_y - s * dir_z;
    let uz = s * dir_y + c * dir_z;

    let horiz_norm = (ux * ux + uz * uz).sqrt();
    let theta_h = if horiz_norm > 1e-9 {
        ux.atan2(uz)
    } else {
        0.0
    };
    res.lateral = if horiz_norm > 1e-9 {
        ux / horiz_norm
    } else {
        0.0
    };

    res.horizontal_radius = ellipsoid_radius(geom, ux, 0.0, uz);
    res.direction_radius = ellipsoid_radius(geom, ux, uy, uz);

    let sr = sample_rate.max(48000.0);
    let mut itd = woodworth_itd_radius(res.horizontal_radius, theta_h);
    let itd_scale = 0.75 + 0.25 * (1.0 - phi_distance);
    itd *= itd_scale;
    res.itd_samples = clamp(itd * sr, -sr * 0.002, sr * 0.002);

    res
}

/// Design pinna multipath response (5-tap FIR with phi-scaled delays).
pub fn design_pinna_response(
    geom: &Geometry,
    pan: f64,
    tuning: &PinnaTuning,
) -> PinnaResponse {
    let mut resp = PinnaResponse::default();
    resp.tap_count = 5;

    let c = 343.0_f64; // speed of sound
    let pan_abs = pan.abs();

    // Pinna anatomical dimensions (metres)
    let _hp = 0.45 * geom.axis_vert;
    let _wp = 0.55 * _hp;
    let l_ec = 0.027;   // ear canal length
    let r_ec = 0.0035;  // ear canal radius
    let l_eff = l_ec + 0.6 * r_ec;
    let d_con = 0.0185;  // concha depth
    let d_cym = 0.0110;  // cymba depth

    // Resonant peak frequencies (quarter-wave resonances)
    resp.peak_freqs[0] = c / (4.0 * l_eff);
    resp.peak_freqs[1] = c / (4.0 * d_con);
    resp.peak_freqs[2] = c / (4.0 * d_cym);

    // Notch frequencies from path-length differences
    let base_dd = 0.019 + 0.0035 * pan_abs + 0.0025 * tuning.elev_norm;
    let azimuth_scale = if pan_abs < 0.4 {
        let t = pan_abs / 0.4;
        0.90 + 0.10 * t
    } else {
        1.0
    };

    let mut dd1 = base_dd * azimuth_scale;
    let mut dd2 = dd1 * 0.78;
    let mut dd3 = dd1 * 1.22;

    if tuning.notch_shift {
        dd1 *= 0.90;
        dd2 *= 0.90;
        dd3 *= 0.90;
    }
    if tuning.elev_plus {
        dd1 *= 0.95;
        dd2 *= 0.95;
        dd3 *= 0.95;
    }

    resp.notch_freqs[0] = c / (2.0 * dd1);
    resp.notch_freqs[1] = c / (2.0 * dd2);
    resp.notch_freqs[2] = c / (2.0 * dd3);

    // Amplitude scaling and HF strength
    resp.amp_scale = 0.35 + 0.65 * (1.0 - (-((pan_abs / 0.12).powi(2))).exp());
    resp.hf_strength = (0.25 + 0.55 * tuning.elev_notch)
        * (0.4 + 0.6 * (1.0 - (-((pan_abs / 0.12).powi(2))).exp()));
    if tuning.elev_plus {
        resp.hf_strength *= 1.35;
    }

    // Tap gains and delays
    let gain_base = resp.amp_scale * (0.55 + 0.45 * tuning.elev_notch);

    // Tap 0: direct path
    resp.tap_delay_sec[0] = 0.0;
    resp.tap_gains[0] = gain_base * (0.12 + 0.10 * (1.0 - pan_abs));

    let sign: f64 = if pan >= 0.0 { 1.0 } else { -1.0 };
    let mut delta = dd1;

    for k in 1..5_usize {
        if k > 1 {
            delta *= PHI;
        }
        resp.tap_delay_sec[k] = delta / c;

        let mut gain = gain_base * INV_PHI.powi(k as i32);

        if pan_abs < 0.4 {
            gain *= 0.7 + 0.3 * (pan_abs / 0.4);
            if k == 1 || k == 2 {
                gain = -gain;
            }
        } else if k % 2 == 1 {
            gain = -gain;
        }

        if sign < 0.0 && k % 2 == 1 {
            gain = -gain;
        }

        resp.tap_gains[k] = gain;
    }

    resp
}

/// Compute torso reflection response (2-tap with LP filtering).
pub fn compute_torso_response(
    _geom: &Geometry,
    sample_rate: f64,
    pan: f64,
    _pitch_rad: f64,
    tuning: &TorsoTuning,
) -> TorsoResponse {
    let mut resp = TorsoResponse::default();

    if tuning.torso_mix <= 1e-6 {
        resp.active = false;
        return resp;
    }

    let sr = sample_rate.max(48000.0);
    let pan_abs = pan.abs();
    let base_ms = clamp(tuning.torso_ms, 0.5, 2.0);
    let primary_delay = (base_ms * 1e-3) * (0.85 + 0.15 * (1.0 - pan_abs));
    let secondary_delay = primary_delay + 0.00011;

    resp.active = true;
    resp.tap_count = 2;
    resp.pan_abs = pan_abs;
    resp.delay_samples[0] = primary_delay * sr;
    resp.delay_samples[1] = secondary_delay * sr;

    let mut mix = clamp(tuning.torso_mix, 0.0, 0.5);
    if pan_abs < 0.3 {
        let t = pan_abs / 0.3;
        mix *= 0.25 + 0.75 * t;
    }
    resp.gains[0] = mix * 0.65;
    resp.gains[1] = mix * 0.45;

    let lp_fc = clamp(tuning.torso_hp_hz, 500.0, 3000.0);
    let alpha = (-TWO_PI * lp_fc / sr).exp();
    resp.lp_alpha[0] = alpha;
    resp.lp_alpha[1] = alpha;

    resp
}

/// Compute distance-dependent gain and lowpass response.
pub fn compute_distance_response(sample_rate: f64, phi_distance: f64) -> DistanceResponse {
    let mut resp = DistanceResponse::default();

    let dist = clamp01(phi_distance);
    let distance_m = 0.8 + 3.2 * dist;
    resp.direct_gain = clamp(0.8 / distance_m, 0.25, 1.0);

    let fc_dist = 14000.0 - 6000.0 * dist;
    let sr = sample_rate.max(48000.0);
    let alpha = (-TWO_PI * fc_dist.max(800.0) / sr).exp();
    resp.lowpass_alpha = clamp(alpha, 0.0, 0.9999);
    resp.hf_weight = 1.0 - 0.30 * dist;

    resp
}

// ---------------------------------------------------------------------------
// PyO3 wrapper: PhiModel
// ---------------------------------------------------------------------------

/// Phi Head Geometry model.
///
/// Wraps the ellipsoid head geometry, pinna multipath, torso, and distance
/// computations behind a stateful object that caches the geometry.
#[pyclass]
#[derive(Clone)]
pub struct PhiModel {
    #[allow(dead_code)]
    cfg: GeometryConfig,
    geom: Geometry,
}

#[pymethods]
impl PhiModel {
    /// Create a new PhiModel with given head parameters.
    #[new]
    #[pyo3(signature = (head_b = 0.0875, phi_ratio = 1.6180339887498948482, nose_ratio = 0.20, front_base = 0.83, back_scale = 1.05))]
    pub fn new(
        head_b: f64,
        phi_ratio: f64,
        nose_ratio: f64,
        front_base: f64,
        back_scale: f64,
    ) -> Self {
        let cfg = GeometryConfig {
            head_b,
            phi_ratio,
            nose_ratio,
            front_base,
            back_scale,
        };
        let geom = build_geometry(&cfg);
        Self { cfg, geom }
    }

    /// Return the computed geometry as a Python dict.
    pub fn geometry<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new_bound(py);
        d.set_item("axis_lat", self.geom.axis_lat)?;
        d.set_item("axis_vert", self.geom.axis_vert)?;
        d.set_item("axis_front", self.geom.axis_front)?;
        d.set_item("axis_front_base", self.geom.axis_front_base)?;
        d.set_item("axis_back", self.geom.axis_back)?;
        d.set_item("nose_len", self.geom.nose_len)?;
        Ok(d)
    }

    /// Compute head result (radii, lateral, ITD) for given spatial parameters.
    #[pyo3(signature = (sr = 48000.0, pan = 0.0, pitch_rad = 0.0, phi_distance = 0.5))]
    pub fn compute_head<'py>(
        &self,
        py: Python<'py>,
        sr: f64,
        pan: f64,
        pitch_rad: f64,
        phi_distance: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let res = compute_head_result(&self.geom, sr, pan, pitch_rad, phi_distance);
        let d = PyDict::new_bound(py);
        d.set_item("horizontal_radius", res.horizontal_radius)?;
        d.set_item("direction_radius", res.direction_radius)?;
        d.set_item("lateral", res.lateral)?;
        d.set_item("itd_samples", res.itd_samples)?;
        Ok(d)
    }

    /// Design pinna multipath response.
    #[pyo3(signature = (pan = 0.0, elev_norm = 0.0, elev_notch = 0.0, notch_shift = false, elev_plus = false))]
    pub fn design_pinna<'py>(
        &self,
        py: Python<'py>,
        pan: f64,
        elev_norm: f64,
        elev_notch: f64,
        notch_shift: bool,
        elev_plus: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let tuning = PinnaTuning {
            elev_norm,
            elev_notch,
            notch_shift,
            elev_plus,
        };
        let resp = design_pinna_response(&self.geom, pan, &tuning);
        let d = PyDict::new_bound(py);
        d.set_item("amp_scale", resp.amp_scale)?;
        d.set_item("hf_strength", resp.hf_strength)?;
        d.set_item("peak_freqs", resp.peak_freqs.to_vec())?;
        d.set_item("notch_freqs", resp.notch_freqs.to_vec())?;
        d.set_item("tap_count", resp.tap_count)?;
        d.set_item("tap_delay_sec", resp.tap_delay_sec.to_vec())?;
        d.set_item("tap_gains", resp.tap_gains.to_vec())?;
        Ok(d)
    }

    /// Compute torso reflection response.
    #[pyo3(signature = (sr = 48000.0, pan = 0.0, pitch_rad = 0.0, torso_mix = 0.0, torso_ms = 1.0, torso_hp_hz = 600.0))]
    pub fn compute_torso<'py>(
        &self,
        py: Python<'py>,
        sr: f64,
        pan: f64,
        pitch_rad: f64,
        torso_mix: f64,
        torso_ms: f64,
        torso_hp_hz: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let tuning = TorsoTuning {
            torso_mix,
            torso_ms,
            torso_hp_hz,
        };
        let resp = compute_torso_response(&self.geom, sr, pan, pitch_rad, &tuning);
        let d = PyDict::new_bound(py);
        d.set_item("active", resp.active)?;
        d.set_item("tap_count", resp.tap_count)?;
        d.set_item("delay_samples", resp.delay_samples.to_vec())?;
        d.set_item("gains", resp.gains.to_vec())?;
        d.set_item("lp_alpha", resp.lp_alpha.to_vec())?;
        d.set_item("pan_abs", resp.pan_abs)?;
        Ok(d)
    }

    /// Compute distance-dependent response.
    #[pyo3(signature = (sr = 48000.0, phi_distance = 0.5))]
    pub fn compute_distance<'py>(
        &self,
        py: Python<'py>,
        sr: f64,
        phi_distance: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let resp = compute_distance_response(sr, phi_distance);
        let d = PyDict::new_bound(py);
        d.set_item("direct_gain", resp.direct_gain)?;
        d.set_item("hf_weight", resp.hf_weight)?;
        d.set_item("lowpass_alpha", resp.lowpass_alpha)?;
        Ok(d)
    }
}

// ---------------------------------------------------------------------------
// Air absorption
// ---------------------------------------------------------------------------

/// Compute frequency-dependent air absorption coefficient (dB/m).
/// Based on ISO 9613-1 simplified model.
/// Returns attenuation factor (linear) for a given frequency and distance.
#[inline]
pub fn compute_air_absorption(freq_hz: f64, distance_m: f64) -> f64 {
    let f_khz = freq_hz * 0.001;
    let alpha_db_per_m = 0.0002 + 0.0002 * f_khz + 0.0006 * f_khz * f_khz;
    let atten_db = alpha_db_per_m * distance_m;
    // Convert dB attenuation to linear gain
    10.0_f64.powf(-atten_db / 20.0)
}

// ---------------------------------------------------------------------------
// Unit tests (Rust side)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_geometry_phi_ratio() {
        let cfg = GeometryConfig::default();
        let geom = build_geometry(&cfg);
        let ratio = geom.axis_vert / geom.axis_lat;
        assert!(
            (ratio - PHI).abs() < 0.001,
            "vertical/lateral ratio {ratio} should be phi"
        );
    }

    #[test]
    fn test_ellipsoid_radius_front_back_asymmetry() {
        let geom = build_geometry(&GeometryConfig::default());
        let r_front = ellipsoid_radius(&geom, 0.0, 0.0, 1.0);
        let r_back = ellipsoid_radius(&geom, 0.0, 0.0, -1.0);
        // back_scale > front_base, so back radius > front radius
        assert!(
            r_back > r_front,
            "back radius {r_back} should exceed front {r_front}"
        );
    }

    #[test]
    fn test_woodworth_itd_zero_at_zero() {
        let itd = woodworth_itd_radius(0.0875, 0.0);
        assert!(itd.abs() < 1e-12, "ITD at theta=0 should be 0");
    }

    #[test]
    fn test_itd_physiological_range() {
        let geom = build_geometry(&GeometryConfig::default());
        let res = compute_head_result(&geom, 48000.0, 1.0, 0.0, 0.5);
        let itd_us = res.itd_samples.abs() / 48000.0 * 1e6;
        assert!(
            itd_us > 200.0 && itd_us < 900.0,
            "ITD {itd_us:.0}us outside physiological range"
        );
    }

    #[test]
    fn test_itd_zero_at_center() {
        let geom = build_geometry(&GeometryConfig::default());
        let res = compute_head_result(&geom, 48000.0, 0.0, 0.0, 0.5);
        assert!(
            res.itd_samples.abs() < 1.0,
            "ITD at center should be ~0, got {}",
            res.itd_samples
        );
    }

    #[test]
    fn test_pinna_5_taps() {
        let geom = build_geometry(&GeometryConfig::default());
        let tuning = PinnaTuning::default();
        let resp = design_pinna_response(&geom, 0.5, &tuning);
        assert_eq!(resp.tap_count, 5);
    }

    #[test]
    fn test_pinna_peak_ranges() {
        let geom = build_geometry(&GeometryConfig::default());
        let tuning = PinnaTuning {
            elev_norm: 0.5,
            elev_notch: 0.5,
            ..Default::default()
        };
        let resp = design_pinna_response(&geom, 0.5, &tuning);
        assert!(resp.peak_freqs[0] > 2500.0 && resp.peak_freqs[0] < 4000.0);
        assert!(resp.peak_freqs[1] > 3500.0 && resp.peak_freqs[1] < 6000.0);
        assert!(resp.peak_freqs[2] > 6000.0 && resp.peak_freqs[2] < 10000.0);
    }

    #[test]
    fn test_torso_inactive_at_zero_mix() {
        let geom = build_geometry(&GeometryConfig::default());
        let tuning = TorsoTuning::default(); // torso_mix = 0
        let resp = compute_torso_response(&geom, 48000.0, 0.5, 0.0, &tuning);
        assert!(!resp.active);
    }

    #[test]
    fn test_torso_active_with_mix() {
        let geom = build_geometry(&GeometryConfig::default());
        let tuning = TorsoTuning {
            torso_mix: 0.3,
            torso_ms: 1.0,
            torso_hp_hz: 600.0,
        };
        let resp = compute_torso_response(&geom, 48000.0, 0.5, 0.0, &tuning);
        assert!(resp.active);
        assert_eq!(resp.tap_count, 2);
    }

    #[test]
    fn test_distance_attenuation() {
        let near = compute_distance_response(48000.0, 0.0);
        let far = compute_distance_response(48000.0, 1.0);
        assert!(
            near.direct_gain > far.direct_gain,
            "near gain {} should exceed far gain {}",
            near.direct_gain,
            far.direct_gain
        );
    }

    #[test]
    fn test_distance_hf_rolloff() {
        let near = compute_distance_response(48000.0, 0.0);
        let far = compute_distance_response(48000.0, 1.0);
        assert!(
            near.hf_weight > far.hf_weight,
            "near hf_weight {} should exceed far {}",
            near.hf_weight,
            far.hf_weight
        );
    }
}
