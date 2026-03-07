//! aureonoise - Modal resonator
//! 2-pole resonant filter bank for physically-informed material simulation.
//! Ported from Mirror7 JUCE: modules/core/modal_engine.hpp

use crate::constants::{PI, TWO_PI};
use crate::math::{clamp, clamp01};

/// ln(1000) — used in T60-to-pole-radius conversion (60 dB = factor 1000)
const LN_1000: f64 = 6.907755278982137;

/// Material preset selector
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ModalPreset {
    Off = 0,
    Wood = 1,
    Metal = 2,
    Glass = 3,
}

impl ModalPreset {
    pub fn from_i32(v: i32) -> Self {
        match v {
            1 => Self::Wood,
            2 => Self::Metal,
            3 => Self::Glass,
            _ => Self::Off,
        }
    }
}

/// Descriptor for one resonant partial (freq in Hz, decay in seconds, gain linear).
/// Used only for preset tables; the engine stores ratio-based `ModeState`.
struct ModeDesc {
    freq: f64,
    decay: f64,
    gain: f64,
}

/// Internal state for a single 2-pole resonant mode.
/// Filter: y[n] = b0 * x[n] + a1 * y[n-1] + a2 * y[n-2]
/// where a1 = 2*r*cos(w), a2 = -r^2, b0 = (1 - r)
#[derive(Clone)]
struct ModeState {
    ratio: f64, // freq relative to base_hz
    decay: f64, // raw T60 from preset (scaled by decay_norm at coefficient build)
    gain: f64,  // linear output gain
    // Coefficients
    a1: f64,
    a2: f64,
    b0: f64,
    // Per-channel state (L)
    y1_l: f64,
    y2_l: f64,
    // Per-channel state (R) — may use detuned coefficients for mirror
    a1_r: f64,
    a2_r: f64,
    b0_r: f64,
    y1_r: f64,
    y2_r: f64,
}

impl ModeState {
    fn new(ratio: f64, decay: f64, gain: f64) -> Self {
        Self {
            ratio,
            decay,
            gain,
            a1: 0.0,
            a2: 0.0,
            b0: 0.0,
            y1_l: 0.0,
            y2_l: 0.0,
            a1_r: 0.0,
            a2_r: 0.0,
            b0_r: 0.0,
            y1_r: 0.0,
            y2_r: 0.0,
        }
    }
}

/// Preset mode tables (ported from C++ `preset_modes`).
/// Values: (freq_hz, decay_seconds, gain_linear).
fn preset_table(preset: ModalPreset) -> Vec<ModeDesc> {
    match preset {
        ModalPreset::Wood => vec![
            ModeDesc { freq: 205.0, decay: 1.60, gain: 1.00 },
            ModeDesc { freq: 312.0, decay: 1.45, gain: 0.88 },
            ModeDesc { freq: 421.0, decay: 1.32, gain: 0.78 },
            ModeDesc { freq: 538.0, decay: 1.18, gain: 0.70 },
            ModeDesc { freq: 660.0, decay: 1.00, gain: 0.60 },
            ModeDesc { freq: 795.0, decay: 0.88, gain: 0.52 },
            ModeDesc { freq: 932.0, decay: 0.80, gain: 0.44 },
            ModeDesc { freq: 1105.0, decay: 0.72, gain: 0.36 },
        ],
        ModalPreset::Metal => vec![
            ModeDesc { freq: 185.0, decay: 2.40, gain: 0.95 },
            ModeDesc { freq: 260.0, decay: 2.20, gain: 0.90 },
            ModeDesc { freq: 340.0, decay: 2.00, gain: 0.85 },
            ModeDesc { freq: 420.0, decay: 1.85, gain: 0.80 },
            ModeDesc { freq: 520.0, decay: 1.70, gain: 0.74 },
            ModeDesc { freq: 630.0, decay: 1.55, gain: 0.68 },
            ModeDesc { freq: 760.0, decay: 1.40, gain: 0.62 },
            ModeDesc { freq: 905.0, decay: 1.25, gain: 0.56 },
            ModeDesc { freq: 1085.0, decay: 1.15, gain: 0.48 },
            ModeDesc { freq: 1265.0, decay: 1.05, gain: 0.42 },
        ],
        ModalPreset::Glass => vec![
            ModeDesc { freq: 480.0, decay: 1.80, gain: 0.80 },
            ModeDesc { freq: 720.0, decay: 1.60, gain: 0.72 },
            ModeDesc { freq: 960.0, decay: 1.44, gain: 0.64 },
            ModeDesc { freq: 1200.0, decay: 1.28, gain: 0.58 },
            ModeDesc { freq: 1440.0, decay: 1.14, gain: 0.52 },
            ModeDesc { freq: 1680.0, decay: 1.02, gain: 0.46 },
            ModeDesc { freq: 1920.0, decay: 0.92, gain: 0.40 },
            ModeDesc { freq: 2160.0, decay: 0.84, gain: 0.35 },
            ModeDesc { freq: 2400.0, decay: 0.76, gain: 0.30 },
        ],
        ModalPreset::Off => vec![],
    }
}

/// Modal resonator bank with stereo decorrelation and feedback.
///
/// Signal path: input + feedback * prev_output -> resonator bank -> output.
/// Right channel modes are detuned by `mirror * 0.02 * freq` for spatial width.
pub struct ModalEngine {
    modes: Vec<ModeState>,
    preset: ModalPreset,
    base_hz: f64,
    decay_scale: f64, // 0..1, maps to C++ decayNorm_
    feedback: f64,     // 0..1
    mirror: f64,       // 0..1, stereo detuning amount
    mix: f64,          // 0..1, wet/dry
    active: bool,
    sr: f64,
    // Feedback state
    prev_out_l: f64,
    prev_out_r: f64,
    // Contralateral spatial mirror
    mirror_pan: f64,       // target mirror position (negated burst centroid)
    mirror_intensity: f64, // burst weight controlling mirror strength
    pub contralateral: f64,    // user param: overall contralateral strength
}

impl ModalEngine {
    /// Create a new modal engine at the given sample rate (defaults to Wood preset, inactive).
    pub fn new(sr: f64) -> Self {
        let sr = if sr > 0.0 { sr } else { 44100.0 };
        Self {
            modes: Vec::new(),
            preset: ModalPreset::Off,
            base_hz: 220.0,
            decay_scale: 0.6,
            feedback: 0.0,
            mirror: 0.0,
            mix: 1.0,
            active: false,
            sr,
            prev_out_l: 0.0,
            prev_out_r: 0.0,
            mirror_pan: 0.0,
            mirror_intensity: 0.0,
            contralateral: 0.0,
        }
    }

    /// Load a material preset, rebuilding the mode bank.
    pub fn set_preset(&mut self, preset: ModalPreset) {
        if self.preset == preset {
            return;
        }
        self.preset = preset;
        self.build_preset_modes();
    }

    /// Update sample rate, recalculating all coefficients.
    pub fn set_sr(&mut self, sr: f64) {
        let sr = if sr > 0.0 { sr } else { 44100.0 };
        if (sr - self.sr).abs() > 1.0 {
            self.sr = sr;
            self.rebuild_coefficients();
        }
    }

    /// Set wet/dry mix (0 = fully dry, 1 = fully wet).
    pub fn set_mix(&mut self, mix: f64) {
        self.mix = clamp01(mix);
    }

    /// Set global decay multiplier (0..1). Maps to C++ decayNorm_.
    pub fn set_decay_scale(&mut self, scale: f64) {
        let scale = clamp01(scale);
        if (scale - self.decay_scale).abs() > 1e-9 {
            self.decay_scale = scale;
            self.rebuild_coefficients();
        }
    }

    /// Set feedback amount (0..1). Output is fed back into the resonator input.
    pub fn set_feedback(&mut self, fb: f64) {
        self.feedback = clamp01(fb);
    }

    /// Set stereo decorrelation (0..1).
    /// Right channel modes are detuned by up to 2% of their frequency.
    pub fn set_mirror(&mut self, mirror: f64) {
        let mirror = clamp01(mirror);
        if (mirror - self.mirror).abs() > 1e-9 {
            self.mirror = mirror;
            self.rebuild_coefficients();
        }
    }

    /// Set base frequency in Hz (20..20000). Modes are ratios of this.
    pub fn set_base_hz(&mut self, hz: f64) {
        let hz = clamp(hz, 20.0, 20000.0);
        if (hz - self.base_hz).abs() > 1e-6 {
            self.base_hz = hz;
            self.rebuild_coefficients();
        }
    }

    /// Enable or disable modal processing.
    pub fn set_active(&mut self, active: bool) {
        self.active = active;
    }

    /// Set contralateral mirror state (called once per block, not per sample).
    /// centroid: burst cluster center of mass (-1..1)
    /// intensity: burst weight controlling mirror strength
    /// amount: user param (0..1) overall contralateral strength
    pub fn set_contralateral(&mut self, centroid: f64, intensity: f64, amount: f64) {
        self.mirror_pan = -centroid;  // mirror = opposite hemisphere
        self.mirror_intensity = intensity;
        self.contralateral = amount;
    }

    /// Number of active resonant modes.
    pub fn mode_count(&self) -> usize {
        self.modes.len()
    }

    /// Process one stereo sample pair. Returns (left, right).
    /// When inactive or empty, passes input through unchanged.
    #[inline]
    pub fn process(&mut self, input_l: f64, input_r: f64) -> (f64, f64) {
        if !self.active || self.modes.is_empty() {
            return (input_l, input_r);
        }

        // Mix input with feedback
        let in_l = input_l + self.feedback * self.prev_out_l;
        let in_r = input_r + self.feedback * self.prev_out_r;

        let mut sum_l = 0.0;
        let mut sum_r = 0.0;

        for m in self.modes.iter_mut() {
            // Left channel
            let y_l = m.b0 * in_l + m.a1 * m.y1_l + m.a2 * m.y2_l;
            m.y2_l = m.y1_l;
            m.y1_l = y_l;

            // Right channel (possibly detuned coefficients)
            let y_r = m.b0_r * in_r + m.a1_r * m.y1_r + m.a2_r * m.y2_r;
            m.y2_r = m.y1_r;
            m.y1_r = y_r;

            sum_l += m.gain * y_l;
            sum_r += m.gain * y_r;
        }

        self.prev_out_l = sum_l;
        self.prev_out_r = sum_r;

        // Contralateral spatial mirror
        let mirror_amt = self.contralateral * self.mirror_intensity;
        if mirror_amt > 1e-6 {
            // mirror_pan: -1 = full left, +1 = full right
            // Route modal output toward the mirror position
            let mp = clamp(self.mirror_pan, -1.0, 1.0);
            // Equal-power pan for the mirrored portion
            let angle = (mp + 1.0) * 0.25 * PI;  // 0 to PI/2
            let mirror_l = angle.cos();
            let mirror_r = angle.sin();

            // Blend: (1-amt)*original + amt*mirrored
            let mono_wet = (sum_l + sum_r) * 0.5;
            let final_l = sum_l * (1.0 - mirror_amt) + mono_wet * mirror_l * mirror_amt;
            let final_r = sum_r * (1.0 - mirror_amt) + mono_wet * mirror_r * mirror_amt;

            let out_l = input_l * (1.0 - self.mix) + final_l * self.mix;
            let out_r = input_r * (1.0 - self.mix) + final_r * self.mix;
            return (out_l, out_r);
        }

        // Normal (no mirror) path
        let out_l = (1.0 - self.mix) * input_l + self.mix * sum_l;
        let out_r = (1.0 - self.mix) * input_r + self.mix * sum_r;

        (out_l, out_r)
    }

    /// Reset all resonator state (delay lines and feedback memory).
    pub fn reset(&mut self) {
        for m in self.modes.iter_mut() {
            m.y1_l = 0.0;
            m.y2_l = 0.0;
            m.y1_r = 0.0;
            m.y2_r = 0.0;
        }
        self.prev_out_l = 0.0;
        self.prev_out_r = 0.0;
        self.mirror_pan = 0.0;
        self.mirror_intensity = 0.0;
    }

    // -- internal -------------------------------------------------------

    /// Build mode bank from the current preset table.
    fn build_preset_modes(&mut self) {
        self.modes.clear();
        let table = preset_table(self.preset);
        if table.is_empty() {
            return;
        }
        let f0 = table[0].freq.max(1.0);
        self.modes.reserve(table.len());
        for desc in &table {
            let ratio = if f0 > 0.0 { desc.freq / f0 } else { 1.0 };
            let ratio = if ratio.is_finite() && ratio > 0.0 { ratio } else { 1.0 };
            self.modes.push(ModeState::new(ratio, desc.decay, desc.gain));
        }
        self.rebuild_coefficients();
    }

    /// Recalculate filter coefficients for all modes.
    /// Left channel uses base_hz * ratio; right channel detunes by mirror amount.
    fn rebuild_coefficients(&mut self) {
        for m in self.modes.iter_mut() {
            let target_l = self.base_hz * m.ratio;
            let freq_l = clamp(target_l, 40.0, 0.45 * self.sr);
            let t60 = (m.decay * (0.4 + 0.6 * self.decay_scale)).max(0.05);
            let r = (-LN_1000 / (t60 * self.sr)).exp();

            // Left channel coefficients
            let omega_l = TWO_PI * freq_l / self.sr;
            m.b0 = 1.0 - r;
            m.a1 = 2.0 * r * omega_l.cos();
            m.a2 = -(r * r);

            // Right channel: detune by mirror * 0.02 * freq
            let target_r = target_l * (1.0 + self.mirror * 0.02);
            let freq_r = clamp(target_r, 40.0, 0.45 * self.sr);
            let omega_r = TWO_PI * freq_r / self.sr;
            m.b0_r = 1.0 - r;
            m.a1_r = 2.0 * r * omega_r.cos();
            m.a2_r = -(r * r);

            // Clear resonator state on coefficient change (matches C++ behaviour)
            m.y1_l = 0.0;
            m.y2_l = 0.0;
            m.y1_r = 0.0;
            m.y2_r = 0.0;
        }
        self.prev_out_l = 0.0;
        self.prev_out_r = 0.0;
    }
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preset_mode_counts() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Wood);
        assert_eq!(eng.mode_count(), 8, "Wood should have 8 modes");

        eng.set_preset(ModalPreset::Metal);
        assert_eq!(eng.mode_count(), 10, "Metal should have 10 modes");

        eng.set_preset(ModalPreset::Glass);
        assert_eq!(eng.mode_count(), 9, "Glass should have 9 modes");

        eng.set_preset(ModalPreset::Off);
        assert_eq!(eng.mode_count(), 0, "Off should have 0 modes");
    }

    #[test]
    fn impulse_decays_to_zero() {
        // For each preset, feed a single impulse and verify all modes
        // decay below threshold within max_decay * sr * 3 samples.
        for preset in [ModalPreset::Wood, ModalPreset::Metal, ModalPreset::Glass] {
            let sr = 44100.0;
            let mut eng = ModalEngine::new(sr);
            eng.set_preset(preset);
            eng.set_active(true);
            eng.set_mix(1.0);
            eng.set_decay_scale(0.6);

            // Find the longest decay in this preset
            let max_decay: f64 = eng
                .modes
                .iter()
                .map(|m| m.decay * (0.4 + 0.6 * eng.decay_scale))
                .fold(0.0, f64::max);

            // Impulse
            let (_, _) = eng.process(1.0, 1.0);

            // Run for 3x the longest T60
            let n = (max_decay * sr * 3.0) as usize;
            let mut last_l = 0.0;
            let mut last_r = 0.0;
            for _ in 0..n {
                let (l, r) = eng.process(0.0, 0.0);
                last_l = l;
                last_r = r;
            }

            assert!(
                last_l.abs() < 0.001 && last_r.abs() < 0.001,
                "{:?}: signal did not decay (L={}, R={})",
                preset,
                last_l,
                last_r,
            );
        }
    }

    #[test]
    fn no_nan_or_inf_under_noise() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Metal);
        eng.set_active(true);
        eng.set_mix(0.7);
        eng.set_feedback(0.3);
        eng.set_mirror(0.5);
        eng.set_decay_scale(0.8);

        // Simple deterministic pseudo-noise (not dependent on rng module)
        let mut state: u64 = 0xDEAD_BEEF;
        for _ in 0..44100 {
            // xorshift64
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let v = (state as f64) / (u64::MAX as f64) * 2.0 - 1.0;
            let (l, r) = eng.process(v, -v);
            assert!(l.is_finite(), "NaN/Inf in left channel");
            assert!(r.is_finite(), "NaN/Inf in right channel");
        }
    }

    #[test]
    fn mirror_zero_identical_lr() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Glass);
        eng.set_active(true);
        eng.set_mix(1.0);
        eng.set_mirror(0.0);

        // Feed identical L/R; output must be identical
        let (_, _) = eng.process(1.0, 1.0);
        for _ in 0..256 {
            let (l, r) = eng.process(0.0, 0.0);
            assert!(
                (l - r).abs() < 1e-15,
                "mirror=0 should produce identical L/R, got L={} R={}",
                l,
                r,
            );
        }
    }

    #[test]
    fn mirror_nonzero_differs_lr() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Glass);
        eng.set_active(true);
        eng.set_mix(1.0);
        eng.set_mirror(0.8);

        // Impulse with identical L/R
        let (_, _) = eng.process(1.0, 1.0);

        // After some samples, L and R should diverge due to detuning
        let mut diff_found = false;
        for _ in 0..512 {
            let (l, r) = eng.process(0.0, 0.0);
            if (l - r).abs() > 1e-6 {
                diff_found = true;
                break;
            }
        }
        assert!(diff_found, "mirror>0 should produce different L/R");
    }

    #[test]
    fn passthrough_when_inactive() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Wood);
        eng.set_active(false);

        let (l, r) = eng.process(0.42, -0.37);
        assert!((l - 0.42).abs() < 1e-15);
        assert!((r - -0.37).abs() < 1e-15);
    }

    #[test]
    fn feedback_increases_sustain() {
        let sr = 44100.0;
        let n = 4000_usize;

        // Without feedback
        let mut eng_nofb = ModalEngine::new(sr);
        eng_nofb.set_preset(ModalPreset::Wood);
        eng_nofb.set_active(true);
        eng_nofb.set_mix(1.0);
        eng_nofb.set_feedback(0.0);
        let (_, _) = eng_nofb.process(1.0, 1.0);
        let mut energy_nofb = 0.0;
        for _ in 0..n {
            let (l, _) = eng_nofb.process(0.0, 0.0);
            energy_nofb += l * l;
        }

        // With feedback
        let mut eng_fb = ModalEngine::new(sr);
        eng_fb.set_preset(ModalPreset::Wood);
        eng_fb.set_active(true);
        eng_fb.set_mix(1.0);
        eng_fb.set_feedback(0.4);
        let (_, _) = eng_fb.process(1.0, 1.0);
        let mut energy_fb = 0.0;
        for _ in 0..n {
            let (l, _) = eng_fb.process(0.0, 0.0);
            energy_fb += l * l;
        }

        assert!(
            energy_fb > energy_nofb,
            "feedback should increase sustain energy ({} vs {})",
            energy_fb,
            energy_nofb,
        );
    }

    #[test]
    fn mix_zero_is_dry() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Metal);
        eng.set_active(true);
        eng.set_mix(0.0);

        let (l, r) = eng.process(0.5, -0.3);
        assert!((l - 0.5).abs() < 1e-12);
        assert!((r - -0.3).abs() < 1e-12);
    }

    #[test]
    fn reset_clears_state() {
        let mut eng = ModalEngine::new(44100.0);
        eng.set_preset(ModalPreset::Metal);
        eng.set_active(true);
        eng.set_mix(1.0);

        // Feed impulse and let it ring
        let (_, _) = eng.process(1.0, 1.0);
        for _ in 0..256 {
            eng.process(0.0, 0.0);
        }

        // Output should be non-zero before reset
        let (l_before, _) = eng.process(0.0, 0.0);
        assert!(l_before.abs() > 1e-10);

        eng.reset();

        // After reset, silence in -> silence out
        let (l_after, r_after) = eng.process(0.0, 0.0);
        assert!(l_after.abs() < 1e-15);
        assert!(r_after.abs() < 1e-15);
    }
}
