//! aureonoise - Burst position modulation
//! Port of beta7_tools/src/burst/engine.cpp
//! Long grains -> center pull (ASMR proximity)
//! Short grains -> edge push (bilateral separation)

use crate::stoch::Hawkes;
use crate::math::clamp01;

/// Result of burst position modulation
pub struct BurstResult {
    /// Modulated pan position [-1, 1]
    pub pan: f64,
    /// Amplitude scaling factor (>= 1.0 during bursts)
    pub amp_scale: f64,
}

/// Burst engine: modulates grain pan/amplitude based on Hawkes intensity.
///
/// When the Hawkes process fires (lambda > base), grains are repositioned:
/// - Long grains pulled toward center (intimate ASMR proximity effect)
/// - Short grains pushed toward edges (bilateral separation for EMDR)
/// - Amplitude boosted proportionally to burst weight and duration
pub struct BurstEngine {
    /// Master enable flag (synced from Params.burst)
    pub enabled: bool,
    /// Floor threshold for burst activation (default 0.35)
    pub floor: f64,
    /// Phi mix coefficient for position blending (default 0.6)
    pub phi_mix: f64,
}

impl BurstEngine {
    pub fn new() -> Self {
        Self {
            enabled: false,
            floor: 0.35,
            phi_mix: 0.6,
        }
    }

    /// Compute burst weight from Hawkes intensity.
    /// Port of C++ compute_weight: (lambda - base) / (base * 4.0), clamped [0, 1].
    /// Returns 0.0 when disabled or when Hawkes is at baseline.
    #[inline]
    pub fn compute_weight(&self, hawkes: &Hawkes) -> f64 {
        if !self.enabled {
            return 0.0;
        }
        let base = hawkes.base;
        if base <= 1.0e-12 {
            return 0.0;
        }
        clamp01((hawkes.intensity() - base) / (base * 4.0))
    }

    /// Apply position modulation to a grain's pan and amplitude.
    ///
    /// Port of beta7_tools/src/burst/engine.cpp::apply_position
    ///
    /// - `hawkes`: current Hawkes process state
    /// - `pan`: raw pan position [-1, 1]
    /// - `dur_norm`: normalized grain duration [0, 1] (long = 1)
    /// - `gap_norm`: normalized gap between grains [0, 1] (large gap = 1)
    pub fn apply_position(
        &self,
        hawkes: &Hawkes,
        pan: f64,
        dur_norm: f64,
        gap_norm: f64,
    ) -> BurstResult {
        let clamped_pan = pan.clamp(-1.0, 1.0);
        let w = self.compute_weight(hawkes);

        if w < 1.0e-6 {
            return BurstResult {
                pan: clamped_pan,
                amp_scale: 1.0,
            };
        }

        let dur = clamp01(dur_norm);
        let gap = clamp01(gap_norm);
        let sign = if clamped_pan >= 0.0 { 1.0 } else { -1.0 };
        let mut mag = clamped_pan.abs();

        // Long grains -> center pull (intimate, ASMR proximity)
        let center_pull = w * (0.35 + 0.45 * dur);
        mag = mag * (1.0 - center_pull) + center_pull * 0.12;

        // Short grains -> edge push (bilateral separation)
        let edge_push = w * (0.35 + 0.30 * (1.0 - dur));
        mag = (mag + edge_push * (1.0 - gap)).clamp(0.0, 1.0);

        // Amplitude boost: louder during bursts, proportional to duration
        let amp_scale = 1.0 + w * (0.6 + 0.4 * dur);

        BurstResult {
            pan: (sign * mag).clamp(-1.0, 1.0),
            amp_scale,
        }
    }
}

impl Default for BurstEngine {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    fn make_hot_hawkes() -> Hawkes {
        // Pump Hawkes well above baseline so compute_weight > 0
        let mut h = Hawkes::new(None, None);
        let mut rng = Rng::new(42);
        for _ in 0..200 {
            h.tick(0.001, &mut rng);
        }
        // Force lambda high to guarantee weight
        h.lambda = h.base * 3.0;
        h
    }

    #[test]
    fn test_weight_zero_when_disabled() {
        let be = BurstEngine::new(); // enabled = false
        let h = make_hot_hawkes();
        assert_eq!(be.compute_weight(&h), 0.0);
    }

    #[test]
    fn test_weight_zero_at_baseline() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let mut h = Hawkes::new(None, None);
        h.lambda = h.base; // exactly at baseline
        assert!(be.compute_weight(&h) < 1.0e-6);
    }

    #[test]
    fn test_weight_positive_above_baseline() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let h = make_hot_hawkes();
        let w = be.compute_weight(&h);
        assert!(w > 0.0, "weight should be positive, got {}", w);
        assert!(w <= 1.0, "weight should be <= 1.0, got {}", w);
    }

    #[test]
    fn test_passthrough_when_disabled() {
        let be = BurstEngine::new(); // disabled
        let h = make_hot_hawkes();
        let r = be.apply_position(&h, 0.7, 0.5, 0.5);
        assert!((r.pan - 0.7).abs() < 1e-10);
        assert!((r.amp_scale - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_long_grain_center_pull() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let h = make_hot_hawkes();
        // Long grain (dur_norm = 1.0) should pull toward center
        let r = be.apply_position(&h, 0.8, 1.0, 0.5);
        // Pan magnitude should be reduced from 0.8
        assert!(r.pan.abs() < 0.8, "center pull failed: pan={}", r.pan);
        assert!(r.pan > 0.0, "sign should be preserved");
    }

    #[test]
    fn test_short_grain_edge_push() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let h = make_hot_hawkes();
        // Short grain (dur_norm = 0.0), small gap (gap_norm = 0.0)
        let r = be.apply_position(&h, 0.3, 0.0, 0.0);
        // Pan magnitude should increase from edge push
        assert!(r.pan.abs() > 0.3, "edge push failed: pan={}", r.pan);
    }

    #[test]
    fn test_amp_scale_boosted_during_burst() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let h = make_hot_hawkes();
        let r = be.apply_position(&h, 0.5, 0.5, 0.5);
        assert!(r.amp_scale > 1.0, "amp should be boosted, got {}", r.amp_scale);
    }

    #[test]
    fn test_pan_bounded() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let h = make_hot_hawkes();
        // Extreme inputs
        for pan in [-1.0, -0.5, 0.0, 0.5, 1.0] {
            for dur in [0.0, 0.5, 1.0] {
                for gap in [0.0, 0.5, 1.0] {
                    let r = be.apply_position(&h, pan, dur, gap);
                    assert!(r.pan >= -1.0 && r.pan <= 1.0,
                        "pan out of bounds: {} for ({}, {}, {})", r.pan, pan, dur, gap);
                    assert!(r.amp_scale >= 1.0,
                        "amp_scale < 1: {} for ({}, {}, {})", r.amp_scale, pan, dur, gap);
                }
            }
        }
    }

    #[test]
    fn test_sign_preserved() {
        let mut be = BurstEngine::new();
        be.enabled = true;
        let h = make_hot_hawkes();
        let r_pos = be.apply_position(&h, 0.5, 0.5, 0.5);
        let r_neg = be.apply_position(&h, -0.5, 0.5, 0.5);
        assert!(r_pos.pan >= 0.0, "positive sign lost");
        assert!(r_neg.pan <= 0.0, "negative sign lost");
    }
}
