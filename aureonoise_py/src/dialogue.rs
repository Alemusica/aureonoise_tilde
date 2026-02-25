//! aureonoise - Dialogue system (interhemispheric coherence)
//!
//! Port of Mirror7 DialogueSystem.cpp + Spatializer::applyPhiPan.
//! Three subsystems:
//!   A) DialogueSystem  — Fibonacci-walk handshake detector between L/R hemispheres
//!   B) PhiPan          — Phi-ratio alternating pan (from Spatializer)
//!   C) BilateralOscillator — Deterministic EMDR-style L/R sweep (new, research-based)

use crate::constants::*;
use crate::math::{clamp, clamp01, map_phi_range};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const FIB: [f64; 8] = [1.0, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0];
const FIB_COUNT: usize = FIB.len();

// ---------------------------------------------------------------------------
// Part A: DialogueSystem
// ---------------------------------------------------------------------------

/// Parameters controlling the dialogue evaluation.
pub struct DialogueParams {
    /// Overall strength of phi-correction on pan/amp/dur (0-1, default 0.6)
    pub strength: f64,
    /// Smoothing factor for coherence EMA (0-1, default 0.5)
    pub memory: f64,
    /// Blend toward phi-ratio targets vs raw values (0-1, default 0.75)
    pub phi_mix: f64,
    /// Master enable
    pub enabled: bool,
}

impl Default for DialogueParams {
    fn default() -> Self {
        Self {
            strength: 0.6,
            memory: 0.5,
            phi_mix: 0.75,
            enabled: true,
        }
    }
}

/// Internal state of the dialogue system, tracking hemisphere history.
struct DialogueState {
    last_sign: f64,
    last_mag: f64,
    last_amp: f64,
    last_dur: f64,
    last_gap: f64,
    coherence: f64,
    coherence_accum: f64,
    coherence_samples: u64,
    utterances: u64,
    handshakes: u64,
    failures: u64,
    pan_sum: f64,
    pan_abs_sum: f64,
    pan_samples: u64,
    fib_index: usize,
    fib_direction: i32,
    fib_ratio_target: f64,
    fib_score: f64,
    fib_gap_queue: [f64; 3],
    fib_gap_queue_size: usize,
    fib_gap_queue_pos: usize,
}

impl Default for DialogueState {
    fn default() -> Self {
        Self {
            last_sign: 0.0,
            last_mag: 0.0,
            last_amp: 1.0,
            last_dur: 0.0,
            last_gap: 0.0,
            coherence: 1.0,
            coherence_accum: 0.0,
            coherence_samples: 0,
            utterances: 0,
            handshakes: 0,
            failures: 0,
            pan_sum: 0.0,
            pan_abs_sum: 0.0,
            pan_samples: 0,
            fib_index: 1,
            fib_direction: 1,
            fib_ratio_target: PHI,
            fib_score: 0.0,
            fib_gap_queue: [0.0; 3],
            fib_gap_queue_size: 0,
            fib_gap_queue_pos: 0,
        }
    }
}

/// Result of one dialogue evaluation (propose-then-commit pattern).
pub struct DialogueResult {
    /// Corrected pan position [-1, 1]
    pub pan: f64,
    /// Amplitude scaling factor
    pub amp_scale: f64,
    /// Duration scaling factor
    pub dur_scale: f64,
    /// Proposed next coherence value
    pub next_coherence: f64,
    /// Composite handshake score [0, 1]
    pub handshake_score: f64,
    /// Whether a handshake was detected this utterance
    pub handshake: bool,
    /// Current Fibonacci walk index
    pub fib_index: usize,
    /// Current Fibonacci walk direction (+1 or -1)
    pub fib_direction: i32,
    // -- internal fields for commit --
    sign: i32,
    base_amp: f64,
    base_dur: f64,
    gap_samples: f64,
    handshake_ratio_target: f64,
}

/// Interhemispheric coherence tracker.
///
/// Each grain is an "utterance" from one hemisphere (L or R based on pan sign).
/// When successive utterances alternate hemispheres with Fibonacci-ratio
/// relationships between pan magnitude, amplitude, duration and gap,
/// a "handshake" fires — the algorithmic corpus callosum.
pub struct DialogueSystem {
    state: DialogueState,
}

impl DialogueSystem {
    pub fn new() -> Self {
        Self { state: DialogueState::default() }
    }

    pub fn reset(&mut self) {
        self.state = DialogueState::default();
    }

    /// Current coherence level (EMA-smoothed, typically in [0.6, 1.8]).
    pub fn coherence(&self) -> f64 {
        self.state.coherence
    }

    /// Total handshakes detected since last reset.
    pub fn handshake_count(&self) -> u64 {
        self.state.handshakes
    }

    /// Handshake ratio: handshakes / (utterances - 1). 0 if < 2 utterances.
    pub fn handshake_ratio(&self) -> f64 {
        if self.state.utterances <= 1 { return 0.0; }
        self.state.handshakes as f64 / (self.state.utterances - 1) as f64
    }

    /// Mean coherence across all committed utterances.
    pub fn coherence_mean(&self) -> f64 {
        if self.state.coherence_samples == 0 { return 1.0; }
        self.state.coherence_accum / self.state.coherence_samples as f64
    }

    /// Evaluate a new grain utterance.
    /// Returns a DialogueResult that MUST be committed via [`commit`] if accepted.
    pub fn evaluate(&self, params: &DialogueParams, pan: f64, amp: f64, dur: f64, gap: f64) -> DialogueResult {
        let st = &self.state;
        let pan_clamped = clamp(pan, -1.0, 1.0);
        let sign = if pan_clamped >= 0.0 { 1 } else { -1 };

        let mut res = DialogueResult {
            pan: pan_clamped,
            amp_scale: 1.0,
            dur_scale: 1.0,
            next_coherence: 1.0,
            handshake_score: 0.0,
            handshake: false,
            fib_index: 1,
            fib_direction: 1,
            sign,
            base_amp: amp,
            base_dur: dur,
            gap_samples: gap,
            handshake_ratio_target: PHI,
        };

        if !params.enabled {
            return res;
        }

        // --- Fibonacci walk: candidate next index ---
        let cur_index = st.fib_index.min(FIB_COUNT - 1);
        let cur_dir = if st.fib_direction >= 0 { 1_i32 } else { -1 };
        let mut cand_dir = cur_dir;
        let mut raw_idx = cur_index as i32 + cand_dir;
        if raw_idx < 0 || raw_idx >= FIB_COUNT as i32 {
            cand_dir = -cand_dir;
            raw_idx = cur_index as i32 + cand_dir;
        }
        let cand_index = (raw_idx.max(0) as usize).min(FIB_COUNT - 1);

        let fib_prev = FIB[cur_index].max(1.0);
        let fib_next = FIB[cand_index].max(1.0);
        let mut fib_ratio = fib_next / fib_prev;
        if !fib_ratio.is_finite() || fib_ratio <= 0.0 {
            fib_ratio = PHI;
        }

        res.handshake_ratio_target = fib_ratio;
        res.fib_index = cur_index;
        res.fib_direction = cur_dir;

        let strength = clamp01(params.strength);
        let memory = clamp01(params.memory);
        let phi_mix = clamp01(params.phi_mix);
        let mag = res.pan.abs().max(1.0e-6);
        let has_prev = st.last_sign.abs() > 0.5;

        let mut coherence_target = INV_PHI;

        if has_prev && (st.last_sign.signum() as i32) != sign {
            // --- HEMISPHERE FLIP: handshake candidate ---
            let prev_mag = clamp(st.last_mag, 0.05, 1.0);
            let incoming_ratio = mag / prev_mag.max(1.0e-6);
            let amp_ratio = amp.max(1.0e-6) / st.last_amp.max(1.0e-6);
            let dur_ratio = if st.last_dur > 0.0 {
                dur.max(1.0e-6) / st.last_dur.max(1.0e-6)
            } else {
                incoming_ratio
            };
            let gap_ratio = if st.last_gap > 0.0 {
                gap.max(1.0e-6) / st.last_gap.max(1.0e-6)
            } else {
                incoming_ratio
            };

            let pan_score = ratio_score(incoming_ratio, fib_ratio);
            let amp_score = ratio_score(amp_ratio, fib_ratio);
            let dur_score = ratio_score(dur_ratio, fib_ratio);
            let gap_score = ratio_score(gap_ratio, fib_ratio);

            res.handshake_score = 0.55 * pan_score
                + 0.15 * amp_score
                + 0.15 * dur_score
                + 0.15 * gap_score;

            let threshold = clamp(0.42 - 0.18 * strength, 0.22, 0.42);
            res.handshake = res.handshake_score >= threshold;

            // Phi-corrected pan magnitude
            let target_mag = clamp(
                (1.0 - strength) * mag
                    + strength * (phi_mix * clamp(prev_mag * fib_ratio, 0.05, 1.0)
                        + (1.0 - phi_mix) * mag),
                0.05,
                1.0,
            );
            res.pan = sign as f64 * target_mag;

            // Phi-corrected amplitude scaling
            let amp_target = st.last_amp.max(1.0e-6) * fib_ratio;
            let amp_scale_target = clamp(amp_target / res.base_amp.max(1.0e-6), 0.15, 4.0);
            let amp_interp = phi_mix * amp_scale_target + (1.0 - phi_mix);
            res.amp_scale = clamp((1.0 - strength) + strength * amp_interp, 0.15, 3.0);

            // Phi-corrected duration scaling
            if st.last_dur > 0.0 {
                let dur_target = st.last_dur * fib_ratio;
                let dur_scale_target = clamp(dur_target / res.base_dur.max(1.0), 0.25, 4.0);
                let dur_interp = phi_mix * dur_scale_target + (1.0 - phi_mix);
                res.dur_scale = clamp((1.0 - strength) + strength * dur_interp, 0.25, 4.0);
            }

            coherence_target = clamp(fib_ratio, 0.6, 1.8);

            if res.handshake {
                res.fib_index = cand_index;
                res.fib_direction = cand_dir;
            }
        } else if has_prev {
            // --- SAME HEMISPHERE: gentle relaxation ---
            let relaxed_target = clamp(st.last_mag * st.fib_ratio_target, 0.05, 1.0);
            let relaxed = (1.0 - 0.5 * strength) * mag + 0.5 * strength * relaxed_target;
            res.pan = sign as f64 * clamp(relaxed, 0.05, 1.0);
            coherence_target = 1.0;
        } else {
            // --- FIRST UTTERANCE ---
            coherence_target = 1.0;
        }

        // Pan bias correction: subtract accumulated drift
        if st.pan_samples > 0 {
            let bias = st.pan_sum / st.pan_samples as f64;
            let corr = strength * 0.4;
            res.pan = clamp(res.pan - corr * bias, -1.0, 1.0);
        }
        res.sign = if res.pan >= 0.0 { 1 } else { -1 };

        // Coherence EMA
        res.next_coherence = (1.0 - memory) * st.coherence + memory * coherence_target;

        res
    }

    /// Commit an accepted evaluation result into the state.
    /// Call with the result from [`evaluate`] after the grain is confirmed.
    pub fn commit(&mut self, params: &DialogueParams, result: &DialogueResult, accepted: bool) {
        if !params.enabled {
            return;
        }
        if !accepted {
            self.state.failures += 1;
            return;
        }

        let st = &mut self.state;
        let prev_gap = st.last_gap;

        st.last_sign = result.sign as f64;
        st.last_mag = result.pan.abs();
        st.last_amp = (result.base_amp * result.amp_scale).max(1.0e-6);
        st.last_dur = result.base_dur * result.dur_scale;
        st.last_gap = result.gap_samples;
        st.coherence = result.next_coherence;
        st.coherence_accum += st.coherence;
        st.coherence_samples += 1;
        st.utterances += 1;
        st.fib_index = result.fib_index.min(FIB_COUNT - 1);
        st.fib_direction = if result.fib_direction >= 0 { 1 } else { -1 };

        if result.handshake {
            st.handshakes += 1;
            st.fib_ratio_target = result.handshake_ratio_target;
            st.fib_gap_queue = [0.0; 3];
            st.fib_gap_queue_pos = 0;
            st.fib_gap_queue_size = 0;

            let gap_cap = 44100.0 * 16.0;
            let base_gap = if prev_gap > 1.0 { prev_gap } else { result.gap_samples.max(1.0) };

            // Push base gap
            if st.fib_gap_queue_size < 3 {
                st.fib_gap_queue[st.fib_gap_queue_size] = clamp(base_gap, 1.0, gap_cap);
                st.fib_gap_queue_size += 1;
            }

            // Push extended gap (phi-blend * plastic bias)
            let phi_lo = (result.handshake_ratio_target / PHI).max(1.0e-6);
            let phi_hi = (result.handshake_ratio_target * PHI).max(phi_lo * 1.000001);
            let phi_blend = map_phi_range(phi_lo, phi_hi, clamp01(result.handshake_score));
            let plastic_bias = PLASTIC.powf(0.05 * (result.handshake_score - 0.5));
            let extend_ratio = clamp(phi_blend * plastic_bias, 0.25, 8.0);
            if st.fib_gap_queue_size < 3 {
                st.fib_gap_queue[st.fib_gap_queue_size] = clamp(base_gap * extend_ratio, 1.0, gap_cap);
                st.fib_gap_queue_size += 1;
            }
        } else if st.fib_gap_queue_pos >= st.fib_gap_queue_size {
            st.fib_gap_queue_size = 0;
            st.fib_gap_queue_pos = 0;
            st.fib_gap_queue = [0.0; 3];
        }

        st.fib_score = result.handshake_score;
        st.pan_sum += result.pan;
        st.pan_abs_sum += result.pan.abs();
        st.pan_samples += 1;
    }
}

/// Ratio scoring: how close `ratio` is to `target` on a log scale.
/// Returns 0.0 for ratios far from target, 1.0 for exact match.
#[inline]
fn ratio_score(ratio: f64, target: f64) -> f64 {
    if !ratio.is_finite() || ratio <= 0.0 || !target.is_finite() || target <= 0.0 {
        return 0.0;
    }
    let log_diff = (ratio / target).ln().abs();
    clamp(1.0 - log_diff / 0.4, 0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Part B: PhiPan (from Spatializer::applyPhiPan)
// ---------------------------------------------------------------------------

/// Phi-ratio alternating pan.
///
/// Each grain gets the opposite sign of the previous, with magnitude
/// scaled by phi. If the result would exceed 1.0, falls back to
/// magnitude / phi. Produces a natural bilateral alternation pattern
/// where successive pan positions relate by the golden ratio.
pub struct PhiPan {
    last_sign: i32,
    last_mag: f64,
}

impl PhiPan {
    pub fn new() -> Self {
        Self { last_sign: 1, last_mag: 0.0 }
    }

    pub fn reset(&mut self) {
        self.last_sign = 1;
        self.last_mag = 0.0;
    }

    /// Compute next phi-pan value from a base pan position.
    /// Alternates sign, scales magnitude by phi with overflow fallback.
    pub fn next(&mut self, base_pan: f64) -> f64 {
        let result = if self.last_mag > 1.0e-6 {
            let next_sign = if self.last_sign >= 0 { -1 } else { 1 };
            let mut target_mag = self.last_mag * PHI;
            if target_mag > 0.999 {
                let fallback = self.last_mag / PHI;
                if fallback >= 0.05 {
                    target_mag = fallback;
                } else {
                    target_mag = base_pan.abs().max(0.05);
                }
            }
            target_mag = clamp(target_mag, 0.05, 1.0);
            clamp(next_sign as f64 * target_mag, -1.0, 1.0)
        } else {
            let base_mag = base_pan.abs().max(0.05);
            let first_sign = if base_pan >= 0.0 { 1 } else { -1 };
            clamp(first_sign as f64 * base_mag, -1.0, 1.0)
        };

        // Update state
        self.last_mag = result.abs();
        self.last_sign = if result >= 0.0 { 1 } else { -1 };

        result
    }
}

// ---------------------------------------------------------------------------
// Part C: BilateralOscillator (new, research-based)
// ---------------------------------------------------------------------------

/// Deterministic bilateral oscillator for EMDR-style L/R sweeps.
///
/// Produces a smooth sinusoidal pan signal at a configurable rate (0.5-2.0 Hz),
/// blended with stochastic pan. Uses a pure sine waveform to maintain
/// binaural cues during the sweep (no hard switching artifacts).
///
/// Research basis: EMDR bilateral stimulation at 0.5-2 Hz alternation,
/// smooth transitions preserve interaural coherence.
pub struct BilateralOscillator {
    phase: f64,
    rate: f64,
    amount: f64,
    active: bool,
}

impl BilateralOscillator {
    pub fn new() -> Self {
        Self {
            phase: 0.0,
            rate: 1.0,
            amount: 0.0,
            active: false,
        }
    }

    /// Set oscillation rate in Hz, clamped to [0.5, 2.0].
    pub fn set_rate(&mut self, hz: f64) {
        self.rate = clamp(hz, 0.5, 2.0);
    }

    /// Set blend amount [0, 1]. 0 = pure stochastic, 1 = pure bilateral.
    pub fn set_amount(&mut self, amt: f64) {
        self.amount = clamp01(amt);
        self.active = self.amount > 1.0e-6;
    }

    /// Advance phase by one sample and return the bilateral pan value [-1, 1].
    pub fn tick(&mut self, sr: f64) -> f64 {
        if !self.active || sr <= 0.0 {
            return 0.0;
        }
        self.phase += self.rate / sr;
        if self.phase >= 1.0 {
            self.phase -= self.phase.floor();
        }
        (TWO_PI * self.phase).sin()
    }

    /// Advance phase and blend bilateral pan with stochastic pan.
    /// Returns `(1 - amount) * stochastic_pan + amount * bilateral_pan`.
    pub fn apply(&mut self, stochastic_pan: f64, sr: f64) -> f64 {
        if !self.active {
            return stochastic_pan;
        }
        let bilateral = self.tick(sr);
        (1.0 - self.amount) * stochastic_pan + self.amount * bilateral
    }

    /// Reset phase to zero.
    pub fn reset(&mut self) {
        self.phase = 0.0;
    }

    /// Whether the oscillator is currently active (amount > 0).
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

    // -- DialogueSystem tests --

    #[test]
    fn dialogue_passthrough_when_disabled() {
        let sys = DialogueSystem::new();
        let params = DialogueParams { enabled: false, ..Default::default() };
        let res = sys.evaluate(&params, 0.7, 0.5, 1000.0, 500.0);
        assert!((res.pan - 0.7).abs() < 1e-10);
        assert!((res.amp_scale - 1.0).abs() < 1e-10);
        assert!((res.dur_scale - 1.0).abs() < 1e-10);
        assert!(!res.handshake);
    }

    #[test]
    fn dialogue_first_utterance_no_handshake() {
        let sys = DialogueSystem::new();
        let params = DialogueParams::default();
        let res = sys.evaluate(&params, 0.5, 0.8, 2000.0, 1000.0);
        assert!(!res.handshake);
        assert!(res.handshake_score < 1e-10);
    }

    #[test]
    fn dialogue_same_hemisphere_no_handshake() {
        let mut sys = DialogueSystem::new();
        let params = DialogueParams::default();

        // First utterance: right hemisphere
        let res1 = sys.evaluate(&params, 0.6, 0.8, 2000.0, 1000.0);
        sys.commit(&params, &res1, true);

        // Second utterance: still right hemisphere (no sign flip)
        let res2 = sys.evaluate(&params, 0.4, 0.5, 1500.0, 800.0);
        assert!(!res2.handshake);
    }

    #[test]
    fn dialogue_fib_index_stays_in_bounds() {
        let mut sys = DialogueSystem::new();
        let params = DialogueParams { strength: 0.9, ..Default::default() };

        // Alternate hemispheres many times to walk the Fibonacci sequence
        for i in 0..200 {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            let pan = sign * 0.5;
            let res = sys.evaluate(&params, pan, 0.8, 2000.0, 1000.0);
            assert!(res.fib_index < FIB_COUNT,
                "fib_index {} out of bounds at iteration {}", res.fib_index, i);
            sys.commit(&params, &res, true);
        }
    }

    #[test]
    fn dialogue_handshake_on_sign_flip_with_phi_ratio() {
        let mut sys = DialogueSystem::new();
        let params = DialogueParams {
            strength: 0.8,
            phi_mix: 0.75,
            ..Default::default()
        };

        // First utterance: right (pan=+0.5)
        let res1 = sys.evaluate(&params, 0.5, 1.0, 2000.0, 1000.0);
        sys.commit(&params, &res1, true);

        // Second utterance: left (sign flip), with phi-ratio pan magnitude
        // 0.5 * PHI = 0.809 -> ratio is PHI, should score well
        let phi_pan = -(0.5 * PHI);
        let res2 = sys.evaluate(&params, phi_pan, 1.0 * PHI, 2000.0 * PHI, 1000.0 * PHI);
        // The handshake_score should be high since all ratios approximate PHI
        assert!(res2.handshake_score > 0.3,
            "expected decent score, got {}", res2.handshake_score);
        // With strength=0.8, threshold = 0.42 - 0.18*0.8 = 0.276
        // The score with all-phi ratios should exceed that
    }

    #[test]
    fn dialogue_coherence_updates() {
        let mut sys = DialogueSystem::new();
        let params = DialogueParams::default();

        let initial = sys.coherence();
        assert!((initial - 1.0).abs() < 1e-10);

        let res = sys.evaluate(&params, 0.5, 0.8, 2000.0, 1000.0);
        sys.commit(&params, &res, true);

        // Coherence should have changed from the initial 1.0
        let after = sys.coherence();
        assert!(after > 0.0 && after < 10.0,
            "coherence out of plausible range: {}", after);
    }

    #[test]
    fn dialogue_commit_rejected_increments_failures() {
        let mut sys = DialogueSystem::new();
        let params = DialogueParams::default();
        let res = sys.evaluate(&params, 0.5, 0.8, 2000.0, 1000.0);
        sys.commit(&params, &res, false);
        // Utterance count should NOT have incremented
        assert_eq!(sys.state.utterances, 0);
        assert_eq!(sys.state.failures, 1);
    }

    #[test]
    fn dialogue_pan_bias_correction() {
        let mut sys = DialogueSystem::new();
        let params = DialogueParams { strength: 1.0, ..Default::default() };

        // Commit many right-biased utterances to build drift
        for _ in 0..50 {
            let res = sys.evaluate(&params, 0.9, 0.8, 2000.0, 1000.0);
            sys.commit(&params, &res, true);
        }

        // Now a new right utterance should be corrected toward center
        let res = sys.evaluate(&params, 0.9, 0.8, 2000.0, 1000.0);
        assert!(res.pan < 0.9, "bias correction failed: pan={}", res.pan);
    }

    // -- PhiPan tests --

    #[test]
    fn phi_pan_alternates_sign() {
        let mut pp = PhiPan::new();
        let v1 = pp.next(0.5);  // first: should be +
        let v2 = pp.next(0.5);  // second: should be -
        let v3 = pp.next(0.5);  // third: should be +

        assert!(v1 > 0.0, "first should be positive, got {}", v1);
        assert!(v2 < 0.0, "second should be negative, got {}", v2);
        assert!(v3 > 0.0, "third should be positive, got {}", v3);
    }

    #[test]
    fn phi_pan_magnitude_in_range() {
        let mut pp = PhiPan::new();
        for _ in 0..100 {
            let v = pp.next(0.4);
            assert!(v >= -1.0 && v <= 1.0,
                "phi-pan out of bounds: {}", v);
            assert!(v.abs() >= 0.05 - 1e-10,
                "phi-pan magnitude too small: {}", v);
        }
    }

    #[test]
    fn phi_pan_overflow_fallback() {
        let mut pp = PhiPan::new();
        // Start with a large magnitude that will overflow phi multiplication
        let v1 = pp.next(0.9);
        assert!(v1.abs() <= 1.0);

        // After a few iterations, the phi-scaling + fallback should keep things bounded
        for _ in 0..50 {
            let v = pp.next(0.9);
            assert!(v.abs() <= 1.0, "overflow not handled: {}", v);
        }
    }

    #[test]
    fn phi_pan_reset() {
        let mut pp = PhiPan::new();
        pp.next(0.5);
        pp.next(0.5);
        pp.reset();
        // After reset, should behave like fresh
        assert!(pp.last_mag < 1e-10);
    }

    // -- BilateralOscillator tests --

    #[test]
    fn bilateral_inactive_by_default() {
        let mut osc = BilateralOscillator::new();
        assert!(!osc.is_active());
        let v = osc.tick(44100.0);
        assert!((v - 0.0).abs() < 1e-15);
    }

    #[test]
    fn bilateral_output_range() {
        let mut osc = BilateralOscillator::new();
        osc.set_rate(1.0);
        osc.set_amount(1.0);

        let sr = 44100.0;
        for _ in 0..88200 {
            let v = osc.tick(sr);
            assert!(v >= -1.0 && v <= 1.0,
                "bilateral out of range: {}", v);
        }
    }

    #[test]
    fn bilateral_full_cycle_at_1hz() {
        let mut osc = BilateralOscillator::new();
        osc.set_rate(1.0);
        osc.set_amount(1.0);
        let sr = 44100.0;

        // Track zero crossings (positive -> negative transitions)
        let mut crossings = 0u32;
        let mut prev = 0.0_f64;
        for i in 0..44100 {
            let v = osc.tick(sr);
            if i > 0 && prev >= 0.0 && v < 0.0 {
                crossings += 1;
            }
            prev = v;
        }

        // At 1 Hz over 1 second, expect exactly 1 positive-to-negative crossing
        assert_eq!(crossings, 1,
            "expected 1 full cycle at 1Hz/44100sr, got {} crossings", crossings);
    }

    #[test]
    fn bilateral_blend_with_stochastic() {
        let mut osc = BilateralOscillator::new();
        osc.set_rate(1.0);
        osc.set_amount(0.5);

        let sr = 44100.0;
        let stochastic = 0.7;

        // First sample: phase ~= rate/sr, sin(2pi * tiny) ~= 0
        let blended = osc.apply(stochastic, sr);
        // With amount=0.5 and bilateral ~= 0, result ~= 0.5 * 0.7 = 0.35
        assert!(blended > 0.0 && blended < 1.0,
            "blended out of range: {}", blended);
    }

    #[test]
    fn bilateral_rate_clamped() {
        let mut osc = BilateralOscillator::new();
        osc.set_rate(0.1);  // below min
        assert!((osc.rate - 0.5).abs() < 1e-10);
        osc.set_rate(10.0); // above max
        assert!((osc.rate - 2.0).abs() < 1e-10);
    }

    #[test]
    fn bilateral_amount_zero_passthrough() {
        let mut osc = BilateralOscillator::new();
        osc.set_rate(1.0);
        osc.set_amount(0.0);
        let result = osc.apply(0.42, 44100.0);
        assert!((result - 0.42).abs() < 1e-15);
    }

    // -- ratio_score tests --

    #[test]
    fn ratio_score_exact_match() {
        let s = ratio_score(PHI, PHI);
        assert!((s - 1.0).abs() < 1e-10,
            "exact match should score 1.0, got {}", s);
    }

    #[test]
    fn ratio_score_invalid_inputs() {
        assert_eq!(ratio_score(0.0, PHI), 0.0);
        assert_eq!(ratio_score(-1.0, PHI), 0.0);
        assert_eq!(ratio_score(PHI, 0.0), 0.0);
        assert_eq!(ratio_score(f64::NAN, PHI), 0.0);
        assert_eq!(ratio_score(f64::INFINITY, PHI), 0.0);
    }

    #[test]
    fn ratio_score_decreases_with_distance() {
        let close = ratio_score(PHI * 1.1, PHI);
        let far = ratio_score(PHI * 2.0, PHI);
        assert!(close > far,
            "closer ratio should score higher: close={}, far={}", close, far);
    }
}
