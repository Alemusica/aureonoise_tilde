//! aureonoise - Near-field DVF (Distance Variation Function)
//! 2nd-order shelving filter for close-distance ILD boost.
//! At 30cm, ILD increases 8-15 dB at high frequencies (DVF research).
//!
//! When a sound source is closer than ~0.5m, the head's acoustic shadow
//! produces a frequency-dependent ILD that grows sharply with 1/r².
//! This module implements a per-ear high-shelf biquad (Audio EQ Cookbook)
//! whose gain tracks the inverse-square near-field excess, producing
//! the intimate, "inside-the-head" spatial impression of ASMR content.

use crate::math::{clamp, clamp01};
use crate::constants::TWO_PI;

/// Near-field Distance Variation Function filter.
///
/// Per-ear 2nd-order high-shelf biquad.  The ipsilateral ear (closer to source)
/// receives a high-frequency boost derived from 1/r² near-field excess;
/// the contralateral ear receives a proportionally reduced boost.
pub struct DvfFilter {
    // --- Per-ear biquad state (2nd order shelving) ---
    // Left ear (ipsilateral boost when pan < 0)
    x1_l: f64,
    x2_l: f64,
    y1_l: f64,
    y2_l: f64,
    // Right ear
    x1_r: f64,
    x2_r: f64,
    y1_r: f64,
    y2_r: f64,
    // --- Biquad coefficients (recalculated on param change) ---
    b0: f64,
    b1: f64,
    b2: f64,
    a1: f64,
    a2: f64,
    // --- Cached parameters ---
    last_distance: f64,
    last_sr: f64,
    active: bool,
}

/// Shelving cutoff frequency (Hz).  Chosen at the lower boundary of
/// the head-shadow region: below ~625 Hz diffraction is negligible.
const DVF_FC: f64 = 625.0;

/// Butterworth Q for the shelving section.
const DVF_Q: f64 = 0.707;

/// Distance (m) beyond which the DVF effect is negligible.
const DVF_FAR_THRESHOLD: f64 = 1.5;

/// Minimum distance (m) to avoid singularity in 1/r².
const DVF_NEAR_CLAMP: f64 = 0.15;

/// Maximum distance (m) accepted by the filter.
const DVF_FAR_CLAMP: f64 = 4.0;

/// Maximum gain (dB) the shelf can reach — safety ceiling.
const DVF_MAX_GAIN_DB: f64 = 18.0;

impl DvfFilter {
    /// Create a new DVF filter, zero-initialised and inactive.
    pub fn new() -> Self {
        Self {
            x1_l: 0.0,
            x2_l: 0.0,
            y1_l: 0.0,
            y2_l: 0.0,
            x1_r: 0.0,
            x2_r: 0.0,
            y1_r: 0.0,
            y2_r: 0.0,
            b0: 1.0,
            b1: 0.0,
            b2: 0.0,
            a1: 0.0,
            a2: 0.0,
            last_distance: -1.0,
            last_sr: -1.0,
            active: false,
        }
    }

    /// Clear all filter state (delay elements).  Coefficients are preserved.
    pub fn reset(&mut self) {
        self.x1_l = 0.0;
        self.x2_l = 0.0;
        self.y1_l = 0.0;
        self.y2_l = 0.0;
        self.x1_r = 0.0;
        self.x2_r = 0.0;
        self.y1_r = 0.0;
        self.y2_r = 0.0;
    }

    /// Recalculate biquad coefficients from source distance (metres) and
    /// sample rate.  Only performs the computation when either parameter has
    /// actually changed.
    ///
    /// Gain formula (dB):
    ///   gain = clamp(15 * (1/d² − 1/1.5²), 0, 18)
    ///
    /// At 30 cm → ~15 * (11.11 − 0.44) = ~160 → clamped to 18 dB
    /// At 50 cm → ~15 * (4.00 − 0.44) = ~53  → clamped to 18 dB
    /// At 80 cm → ~15 * (1.56 − 0.44) = ~16.8 dB
    /// At 1.5 m → ~15 * (0.44 − 0.44) = 0 dB  (threshold)
    pub fn update_params(&mut self, distance: f64, sr: f64) {
        // Early-out: nothing changed.
        if (distance - self.last_distance).abs() < 1.0e-9
            && (sr - self.last_sr).abs() < 1.0e-9
        {
            return;
        }

        self.last_distance = distance;
        self.last_sr = sr;

        let d = clamp(distance, DVF_NEAR_CLAMP, DVF_FAR_CLAMP);
        let sr = if sr > 0.0 { sr } else { 44100.0 };

        // Near-field excess relative to far-field threshold.
        let inv_sq = 1.0 / (d * d);
        let inv_sq_ref = 1.0 / (DVF_FAR_THRESHOLD * DVF_FAR_THRESHOLD);
        let gain_db = clamp(15.0 * (inv_sq - inv_sq_ref), 0.0, DVF_MAX_GAIN_DB);

        self.active = gain_db > 0.01;

        if !self.active {
            // Unity passthrough coefficients.
            self.b0 = 1.0;
            self.b1 = 0.0;
            self.b2 = 0.0;
            self.a1 = 0.0;
            self.a2 = 0.0;
            return;
        }

        // ----------------------------------------------------------
        // High-shelf biquad (Audio EQ Cookbook, Robert Bristow-Johnson)
        // ----------------------------------------------------------
        let a_lin = 10.0_f64.powf(gain_db / 40.0); // sqrt of voltage gain
        let w0 = TWO_PI * DVF_FC / sr;
        let cos_w0 = w0.cos();
        let sin_w0 = w0.sin();
        let alpha = sin_w0 / (2.0 * DVF_Q);

        let two_sqrt_a_alpha = 2.0 * a_lin.sqrt() * alpha;

        let b0 = a_lin * ((a_lin + 1.0) + (a_lin - 1.0) * cos_w0 + two_sqrt_a_alpha);
        let b1 = -2.0 * a_lin * ((a_lin - 1.0) + (a_lin + 1.0) * cos_w0);
        let b2 = a_lin * ((a_lin + 1.0) + (a_lin - 1.0) * cos_w0 - two_sqrt_a_alpha);
        let a0 = (a_lin + 1.0) - (a_lin - 1.0) * cos_w0 + two_sqrt_a_alpha;
        let a1 = 2.0 * ((a_lin - 1.0) - (a_lin + 1.0) * cos_w0);
        let a2 = (a_lin + 1.0) - (a_lin - 1.0) * cos_w0 - two_sqrt_a_alpha;

        // Normalise by a0.
        let inv_a0 = 1.0 / a0;
        self.b0 = b0 * inv_a0;
        self.b1 = b1 * inv_a0;
        self.b2 = b2 * inv_a0;
        self.a1 = a1 * inv_a0;
        self.a2 = a2 * inv_a0;
    }

    /// Process one stereo sample pair through the DVF.
    ///
    /// `pan` ∈ [-1, 1]: negative = source on the left, positive = source on
    /// the right.  The ear closest to the source (ipsilateral) receives the
    /// full high-shelf boost; the far ear (contralateral) receives a reduced
    /// boost, blending toward unity.
    ///
    /// Returns `(out_left, out_right)`.
    #[inline]
    pub fn process(&mut self, left: f64, right: f64, pan: f64) -> (f64, f64) {
        if !self.active {
            return (left, right);
        }

        let pan = clamp(pan, -1.0, 1.0);

        // Ipsilateral mix: 1.0 = full shelf, 0.0 = bypass.
        // pan < 0 → source left  → left ear ipsi (mix_l=1), right ear contra.
        // pan > 0 → source right → right ear ipsi (mix_r=1), left ear contra.
        // At centre (pan=0) both ears get 0.5 shelf.
        let ipsi_l = clamp01(0.5 - 0.5 * pan); // 1 when pan=-1, 0.5 at centre, 0 when pan=+1
        let ipsi_r = clamp01(0.5 + 0.5 * pan); // 0 when pan=-1, 0.5 at centre, 1 when pan=+1

        // --- Left ear biquad ---
        let filtered_l = self.b0 * left + self.b1 * self.x1_l + self.b2 * self.x2_l
            - self.a1 * self.y1_l
            - self.a2 * self.y2_l;
        self.x2_l = self.x1_l;
        self.x1_l = left;
        self.y2_l = self.y1_l;
        self.y1_l = filtered_l;

        // --- Right ear biquad ---
        let filtered_r = self.b0 * right + self.b1 * self.x1_r + self.b2 * self.x2_r
            - self.a1 * self.y1_r
            - self.a2 * self.y2_r;
        self.x2_r = self.x1_r;
        self.x1_r = right;
        self.y2_r = self.y1_r;
        self.y1_r = filtered_r;

        // Blend between filtered (shelved) and dry based on ipsilateral weight.
        let out_l = ipsi_l * filtered_l + (1.0 - ipsi_l) * left;
        let out_r = ipsi_r * filtered_r + (1.0 - ipsi_r) * right;

        (out_l, out_r)
    }

    /// Whether the filter is actively modifying the signal.
    #[inline]
    pub fn is_active(&self) -> bool {
        self.active
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inactive_far_distance() {
        let mut dvf = DvfFilter::new();
        dvf.update_params(3.0, 48000.0);

        assert!(!dvf.is_active(), "DVF should be inactive at 3.0 m");

        // Signal must pass through unmodified.
        let (out_l, out_r) = dvf.process(0.75, -0.25, 0.3);
        assert!(
            (out_l - 0.75).abs() < 1e-15,
            "Left channel modified at far distance"
        );
        assert!(
            (out_r - (-0.25)).abs() < 1e-15,
            "Right channel modified at far distance"
        );
    }

    #[test]
    fn test_active_near_distance() {
        let mut dvf = DvfFilter::new();
        dvf.update_params(0.3, 48000.0);

        assert!(dvf.is_active(), "DVF should be active at 0.3 m");

        // Feed a burst of constant signal so the biquad state ramps up.
        let mut changed = false;
        for _ in 0..256 {
            let (out_l, out_r) = dvf.process(0.5, -0.5, -0.6);
            if (out_l - 0.5).abs() > 1e-6 || (out_r - (-0.5)).abs() > 1e-6 {
                changed = true;
            }
        }

        assert!(changed, "DVF should modify signal at 0.3 m");
    }

    #[test]
    fn test_reset_clears_state() {
        let mut dvf = DvfFilter::new();
        dvf.update_params(0.25, 48000.0);

        // Pump signal through to fill delay elements.
        for _ in 0..512 {
            dvf.process(1.0, -1.0, 0.0);
        }

        // State must be non-zero before reset.
        assert!(
            dvf.x1_l.abs() > 1e-15 || dvf.y1_l.abs() > 1e-15,
            "Expected non-zero state before reset"
        );

        dvf.reset();

        assert!(dvf.x1_l.abs() < 1e-30, "x1_l not cleared");
        assert!(dvf.x2_l.abs() < 1e-30, "x2_l not cleared");
        assert!(dvf.y1_l.abs() < 1e-30, "y1_l not cleared");
        assert!(dvf.y2_l.abs() < 1e-30, "y2_l not cleared");
        assert!(dvf.x1_r.abs() < 1e-30, "x1_r not cleared");
        assert!(dvf.x2_r.abs() < 1e-30, "x2_r not cleared");
        assert!(dvf.y1_r.abs() < 1e-30, "y1_r not cleared");
        assert!(dvf.y2_r.abs() < 1e-30, "y2_r not cleared");
    }
}
