//! Stochastic Resonance — adaptive noise-to-pattern calibration
//!
//! Collins, Chow, Imhoff (1995). Phys Rev E. 52(4):
//! In nonlinear threshold systems, noise at -15 to -20 dB below signal
//! threshold optimizes detection of subthreshold patterns.
//!
//! The dialogue system's coherence metric proxies the brain's pattern
//! detection quality. SR adapts the noise floor gain to maintain optimal
//! conditions for the Collins effect:
//!
//! - Low coherence (weak patterns) → increase noise (SR enhancement)
//! - Optimal coherence → unity gain (SR is working)
//! - High coherence (clear patterns) → reduce noise (signal already clear)
//!
//! Gain range ±2.5 dB spans the Collins -15 to -20 dB window.

use crate::math::clamp;

/// Gain limits: 10^(∓2.5/20) ≈ [0.75, 1.33]
const GAIN_MIN: f64 = 0.75;
const GAIN_MAX: f64 = 1.33;

/// Smoothing tau (seconds) — prevents audible pumping
const GAIN_TAU: f64 = 0.5;

/// Normalized coherence at which SR is optimal
const COH_OPTIMAL: f64 = 0.5;

/// Collins window half-width in dB
const COLLINS_HALF_DB: f64 = 2.5;

/// Stochastic Resonance adaptive noise gain.
///
/// Modulates the noise floor level per-block based on dialogue coherence.
/// The gain is applied to noise samples before they enter the ring buffer.
///
/// The adaptation avoids measuring signal power directly because in this
/// granular synth, noise IS the signal material — measuring output power
/// would create a positive feedback loop that collapses to silence.
/// Instead, coherence (a spatial/temporal metric from the dialogue system)
/// serves as the signal-strength proxy, breaking the feedback loop.
pub struct StochasticResonance {
    /// Target offset in dB (Collins center: -17.5)
    target_db: f64,
    /// Current adaptive noise gain
    gain: f64,
    /// Active flag
    active: bool,
}

impl StochasticResonance {
    pub fn new() -> Self {
        Self {
            target_db: -17.5,
            gain: 1.0,
            active: false,
        }
    }

    /// Set active state. When deactivated, gain resets to 1.0.
    pub fn set_active(&mut self, on: bool) {
        self.active = on;
        if !on { self.gain = 1.0; }
    }

    /// Set Collins target dB (clamped to [-25, -10]).
    pub fn set_target_db(&mut self, db: f64) {
        self.target_db = clamp(db, -25.0, -10.0);
    }

    /// Get target dB (for validation/metadata).
    pub fn target_db(&self) -> f64 {
        self.target_db
    }

    #[inline]
    pub fn is_active(&self) -> bool { self.active }

    /// Current noise gain multiplier. Apply to noise samples before ring write.
    #[inline]
    pub fn noise_gain(&self) -> f64 { self.gain }

    /// Adapt noise gain from dialogue coherence (once per block).
    ///
    /// Coherence in [0.6, 1.8] is normalized to [0, 1].
    /// Gain moves ±2.5 dB around unity based on deviation from optimal.
    ///
    /// - `coherence`: raw dialogue coherence value
    /// - `block_len`: samples in current block
    /// - `sr`: sample rate in Hz
    pub fn adapt(&mut self, coherence: f64, block_len: usize, sr: f64) {
        if !self.active || block_len == 0 { return; }

        let alpha = clamp(block_len as f64 / (GAIN_TAU * sr), 0.0005, 0.1);

        // Coherence [0.6, 1.8] → [0, 1]
        let coh_norm = clamp((coherence - 0.6) / 1.2, 0.0, 1.0);

        // dB offset: +2.5 at low coherence (boost noise), -2.5 at high (reduce)
        let db_offset = 2.0 * COLLINS_HALF_DB * (COH_OPTIMAL - coh_norm);
        let target_gain = (10.0_f64).powf(db_offset / 20.0);
        let clamped = clamp(target_gain, GAIN_MIN, GAIN_MAX);

        // Smooth EMA
        self.gain = (1.0 - alpha) * self.gain + alpha * clamped;
    }

    /// Reset gain to unity.
    pub fn reset(&mut self) {
        self.gain = 1.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inactive_unity_gain() {
        let sr_mod = StochasticResonance::new();
        assert!(!sr_mod.is_active());
        assert!((sr_mod.noise_gain() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_optimal_coherence_near_unity() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_active(true);
        // Optimal coherence: 0.5 normalized = 1.2 raw
        for _ in 0..200 {
            sr_mod.adapt(1.2, 512, 44100.0);
        }
        assert!((sr_mod.noise_gain() - 1.0).abs() < 0.05,
            "At optimal coherence, gain should be ~1.0, got {}", sr_mod.noise_gain());
    }

    #[test]
    fn test_low_coherence_boosts() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_active(true);
        // Low coherence = 0.6 (minimum, normalized 0)
        for _ in 0..500 {
            sr_mod.adapt(0.6, 512, 44100.0);
        }
        assert!(sr_mod.noise_gain() > 1.1,
            "At low coherence, gain should boost, got {}", sr_mod.noise_gain());
        assert!(sr_mod.noise_gain() <= GAIN_MAX,
            "Gain should be clamped, got {}", sr_mod.noise_gain());
    }

    #[test]
    fn test_high_coherence_attenuates() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_active(true);
        // High coherence = 1.8 (maximum, normalized 1)
        for _ in 0..500 {
            sr_mod.adapt(1.8, 512, 44100.0);
        }
        assert!(sr_mod.noise_gain() < 0.9,
            "At high coherence, gain should attenuate, got {}", sr_mod.noise_gain());
        assert!(sr_mod.noise_gain() >= GAIN_MIN,
            "Gain should be clamped, got {}", sr_mod.noise_gain());
    }

    #[test]
    fn test_gain_clamped() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_active(true);
        for _ in 0..10000 {
            sr_mod.adapt(0.6, 512, 44100.0);
        }
        assert!(sr_mod.noise_gain() <= GAIN_MAX + 1e-9);
        assert!(sr_mod.noise_gain() >= GAIN_MIN - 1e-9);
    }

    #[test]
    fn test_deactivate_resets_gain() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_active(true);
        for _ in 0..200 {
            sr_mod.adapt(0.6, 512, 44100.0);
        }
        assert!(sr_mod.noise_gain() > 1.0);
        sr_mod.set_active(false);
        assert!((sr_mod.noise_gain() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_smooth_adaptation() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_active(true);
        // Single step from unity toward low-coherence boost
        sr_mod.adapt(0.6, 512, 44100.0);
        // Should move toward GAIN_MAX but not jump there instantly
        assert!(sr_mod.noise_gain() > 1.0);
        assert!(sr_mod.noise_gain() < GAIN_MAX,
            "Single step should not reach max, got {}", sr_mod.noise_gain());
    }

    #[test]
    fn test_target_db_clamped() {
        let mut sr_mod = StochasticResonance::new();
        sr_mod.set_target_db(-30.0);
        assert!((sr_mod.target_db() - (-25.0)).abs() < 1e-9);
        sr_mod.set_target_db(-5.0);
        assert!((sr_mod.target_db() - (-10.0)).abs() < 1e-9);
    }
}
