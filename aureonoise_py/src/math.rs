//! aureonoise - Math utilities
//! φ-based mapping functions and DSP helpers

use crate::constants::*;

/// Clamp value to [min, max]
#[inline]
pub fn clamp(x: f64, min: f64, max: f64) -> f64 {
    x.max(min).min(max)
}

/// Clamp value to [0, 1]
#[inline]
pub fn clamp01(x: f64) -> f64 {
    clamp(x, 0.0, 1.0)
}

/// Convert decibels to linear gain
#[inline]
pub fn db_to_lin(db: f64) -> f64 {
    10.0_f64.powf(db / 20.0)
}

/// Convert linear gain to decibels
#[inline]
pub fn lin_to_db(lin: f64) -> f64 {
    20.0 * lin.max(TINY).log10()
}

/// Hann window function for phase in [0, 1]
#[inline]
pub fn hann01(phase: f64) -> f64 {
    let p = clamp01(phase);
    0.5 - 0.5 * (TWO_PI * p).cos()
}

/// Soft saturation using tanh
#[inline]
pub fn soft_tanh(x: f64) -> f64 {
    x.tanh()
}

/// φ-based exponential mapping
/// Maps u in [0,1] to [vmin, vmax] using φ as the exponential base
/// This creates a perceptually uniform distribution based on the golden ratio
#[inline]
pub fn map_phi_range(vmin: f64, vmax: f64, u: f64) -> f64 {
    let u = clamp01(u);
    let vmin = vmin.max(1.0e-12);
    let vmax = vmax.max(vmin * 1.000001);
    let k = (vmax / vmin).ln() / PHI.ln();
    vmin * PHI.powf(k * u)
}

/// Map sample rate crush amount to hold samples
/// Uses φ-based mapping for perceptually uniform distribution
#[inline]
pub fn map_sr_hold_base(amt: f64, u: f64) -> i32 {
    let amt = clamp01(amt);
    let steps = map_phi_range(1.0, 64.0, amt);
    let n = steps.round().max(1.0) as i32;
    let jitter = (clamp01(u) * 0.999 * n as f64).floor() as i32;
    (n - jitter).max(1)
}

/// Map bit crush amount to bit depth
#[inline]
pub fn map_bits(amt: f64) -> i32 {
    let amt = clamp01(amt);
    let bits = (16.0 - amt * 12.0).round() as i32;
    bits.clamp(4, 16)
}

/// Quantize sample to given bit depth
#[inline]
pub fn quantize_bits(x: f64, bits: i32) -> f64 {
    let levels = ((1 << (bits - 1)) - 1) as f64;
    if levels <= 0.0 {
        return 0.0;
    }
    (x * levels).round() / levels
}

/// Equal-power pan law
/// Returns (gain_left, gain_right) for pan in [-1, 1]
#[inline]
pub fn pan_equal_power(pan: f64, width: f64) -> (f64, f64) {
    let p = clamp(pan * width, -1.0, 1.0);
    let theta = (p + 1.0) * (PI * 0.25);
    (theta.cos(), theta.sin())
}

/// Woodworth ITD model
/// Returns ITD factor based on azimuth angle
#[inline]
pub fn woodworth_itd(theta: f64) -> f64 {
    let limit = 0.5 * PI;
    let theta = theta.clamp(-limit, limit);
    theta + theta.sin()
}

/// Map ITD in samples based on pan position
#[inline]
pub fn map_itd_samples(sr: f64, itd_us: f64, pan: f64, u: f64) -> f64 {
    let sr = if sr > 0.0 { sr } else { DEFAULT_SR };
    let max_sec = clamp(itd_us, 0.0, 1600.0) * 1.0e-6;
    if max_sec <= 0.0 {
        return 0.0;
    }

    let head_radius = 0.0875; // ~17.5 cm diameter
    let speed_of_sound = 343.0;
    let theta = clamp(pan, -1.0, 1.0) * (PI * 0.5);
    let geom = woodworth_itd(theta);
    let geom_max = woodworth_itd(0.5 * PI);
    
    let mut itd_sec = (head_radius / speed_of_sound) * geom;
    let itd_sec_max = (head_radius / speed_of_sound) * geom_max;
    
    if itd_sec_max > 1.0e-9 {
        itd_sec = (max_sec / itd_sec_max) * itd_sec;
    } else {
        itd_sec = 0.0;
    }

    let jitter = (2.0 * clamp01(u) - 1.0) * 0.08;
    itd_sec = clamp(itd_sec + jitter * max_sec, -max_sec, max_sec);
    itd_sec * sr
}

/// Map ILD in dB based on pan position
#[inline]
pub fn map_ild_db(ild_db_max: f64, pan: f64, u: f64) -> f64 {
    let max_db = clamp(ild_db_max, 0.0, 30.0);
    if max_db <= 1.0e-9 {
        return 0.0;
    }

    let theta = clamp(pan, -1.0, 1.0) * (PI * 0.5);
    let lateral = theta.sin();
    let lateral_mag = lateral.abs().powf(0.9);
    let random = map_phi_range(0.35, 1.0, clamp01((2.0 * u - 1.0).abs()));
    let mut ild = max_db * lateral_mag * random;
    let front_focus = 0.65 + 0.35 * lateral_mag;
    ild *= front_focus;
    
    if lateral >= 0.0 { ild } else { -ild }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clamp() {
        assert_eq!(clamp(0.5, 0.0, 1.0), 0.5);
        assert_eq!(clamp(-0.5, 0.0, 1.0), 0.0);
        assert_eq!(clamp(1.5, 0.0, 1.0), 1.0);
    }

    #[test]
    fn test_map_phi_range() {
        let v = map_phi_range(1.0, 10.0, 0.0);
        assert!((v - 1.0).abs() < 0.001);
        
        let v = map_phi_range(1.0, 10.0, 1.0);
        assert!((v - 10.0).abs() < 0.001);
    }

    #[test]
    fn test_pan_equal_power() {
        let (l, r) = pan_equal_power(0.0, 1.0);
        assert!((l - r).abs() < 0.001); // center = equal
        
        let (l, r) = pan_equal_power(-1.0, 1.0);
        assert!(l > r); // left pan = more left
        
        let (l, r) = pan_equal_power(1.0, 1.0);
        assert!(r > l); // right pan = more right
    }
}
