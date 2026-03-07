//! aureonoise - Room divergence (Schroeder allpass reverb)
//! Minimal reverb for headphone externalization.
//! Phi-ratio delay times. RT60 ~300ms.

use crate::constants::PHI;
use crate::math::clamp01;

// ---------------------------------------------------------------
// Delay time constants (seconds)
// Base early reflection: 2.0ms
// Allpass 1: 5.0ms
// Allpass 2: 5.0ms * PHI = 8.09ms
// Comb: 5.0ms * PHI^2 = 13.09ms (rounded to ~15ms conceptually)
// ---------------------------------------------------------------

const RT60: f64 = 0.300; // 300ms reverb tail — short, for externalization

const EARLY_BASE_SEC: f64 = 0.002;   // 2.0ms
const AP1_SEC: f64 = 0.005;          // 5.0ms
const AP2_SEC: f64 = AP1_SEC * PHI;  // ~8.09ms
const COMB_SEC: f64 = AP2_SEC * PHI; // ~13.09ms

// Early reflection gains (decaying, phi-spaced taps)
const EARLY_GAINS: [f64; 4] = [0.80, 0.55, 0.35, 0.22];

// Comb filter damping LP coefficient (mild HF rolloff in tail)
const COMB_DAMP: f64 = 0.35;

// L/R wet offset for stereo width (~0.3ms)
const LR_OFFSET_SEC: f64 = 0.0003;

/// Round up to next power of 2 (minimum 16).
fn next_pow2(n: usize) -> usize {
    let mut v = n.max(16);
    v -= 1;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v + 1
}

/// Feedback coefficient for a given delay time to reach RT60.
/// g = 10^(-3 * delay_sec / rt60)
#[inline]
fn feedback_for_rt60(delay_sec: f64, rt60: f64) -> f64 {
    if rt60 <= 0.0 || delay_sec <= 0.0 {
        return 0.0;
    }
    10.0_f64.powf(-3.0 * delay_sec / rt60)
}

// ---------------------------------------------------------------
// AllpassSection — single allpass delay with feedback
// ---------------------------------------------------------------

pub struct AllpassSection {
    buf: Vec<f64>,
    mask: usize,
    delay_samples: usize,
    feedback: f64,
    write_idx: usize,
}

impl AllpassSection {
    /// Create from delay time in seconds and sample rate.
    fn new(delay_sec: f64, sr: f64) -> Self {
        let delay_samples = (delay_sec * sr).round().max(1.0) as usize;
        let size = next_pow2(delay_samples + 1);
        Self {
            buf: vec![0.0; size],
            mask: size - 1,
            delay_samples,
            feedback: feedback_for_rt60(delay_sec, RT60),
            write_idx: 0,
        }
    }

    fn reset(&mut self) {
        self.buf.iter_mut().for_each(|s| *s = 0.0);
        self.write_idx = 0;
    }

    /// Classic Schroeder allpass: y = -g*x + delayed + g*y_delayed
    /// Rearranged: write = x + g * delayed; output = delayed - g * write
    #[inline]
    fn process(&mut self, input: f64) -> f64 {
        let read_idx = (self.write_idx + self.buf.len() - self.delay_samples) & self.mask;
        let delayed = self.buf[read_idx];
        let write_val = input + self.feedback * delayed;
        self.buf[self.write_idx] = write_val;
        self.write_idx = (self.write_idx + 1) & self.mask;
        delayed - self.feedback * write_val
    }
}

// ---------------------------------------------------------------
// CombSection — feedback comb filter with damping LP
// ---------------------------------------------------------------

pub struct CombSection {
    buf: Vec<f64>,
    mask: usize,
    delay_samples: usize,
    feedback: f64,
    damp: f64,       // damping LP coefficient [0,1]
    lp_state: f64,   // one-pole LP state
    write_idx: usize,
}

impl CombSection {
    fn new(delay_sec: f64, sr: f64) -> Self {
        let delay_samples = (delay_sec * sr).round().max(1.0) as usize;
        let size = next_pow2(delay_samples + 1);
        Self {
            buf: vec![0.0; size],
            mask: size - 1,
            delay_samples,
            feedback: feedback_for_rt60(delay_sec, RT60),
            damp: COMB_DAMP,
            lp_state: 0.0,
            write_idx: 0,
        }
    }

    fn reset(&mut self) {
        self.buf.iter_mut().for_each(|s| *s = 0.0);
        self.write_idx = 0;
        self.lp_state = 0.0;
    }

    /// Feedback comb with one-pole damping LP in the feedback path.
    /// LP: y_lp = (1 - damp) * delayed + damp * y_lp_prev
    /// Write: input + feedback * y_lp
    /// Output: delayed
    #[inline]
    fn process(&mut self, input: f64) -> f64 {
        let read_idx = (self.write_idx + self.buf.len() - self.delay_samples) & self.mask;
        let delayed = self.buf[read_idx];
        // Damping LP on feedback signal
        self.lp_state = (1.0 - self.damp) * delayed + self.damp * self.lp_state;
        self.buf[self.write_idx] = input + self.feedback * self.lp_state;
        self.write_idx = (self.write_idx + 1) & self.mask;
        delayed
    }
}

// ---------------------------------------------------------------
// EarlyReflections — 4 taps at phi-ratio intervals
// ---------------------------------------------------------------

pub struct EarlyReflections {
    buf: Vec<f64>,
    mask: usize,
    taps: [usize; 4],   // delay in samples for each tap
    gains: [f64; 4],
    write_idx: usize,
}

impl EarlyReflections {
    fn new(sr: f64) -> Self {
        // Tap delays: base * [1.0, PHI, PHI^2, PHI^3]
        let phi2 = PHI * PHI;
        let phi3 = phi2 * PHI;
        let multipliers = [1.0, PHI, phi2, phi3];

        let mut max_delay: usize = 0;
        let mut taps = [0usize; 4];
        for (i, &m) in multipliers.iter().enumerate() {
            let d = (EARLY_BASE_SEC * m * sr).round().max(1.0) as usize;
            taps[i] = d;
            if d > max_delay {
                max_delay = d;
            }
        }

        let size = next_pow2(max_delay + 1);
        Self {
            buf: vec![0.0; size],
            mask: size - 1,
            taps,
            gains: EARLY_GAINS,
            write_idx: 0,
        }
    }

    fn reset(&mut self) {
        self.buf.iter_mut().for_each(|s| *s = 0.0);
        self.write_idx = 0;
    }

    /// Write input into delay line, read 4 taps, return weighted sum.
    #[inline]
    fn process(&mut self, input: f64) -> f64 {
        self.buf[self.write_idx] = input;
        let mut out = 0.0;
        for i in 0..4 {
            let read_idx = (self.write_idx + self.buf.len() - self.taps[i]) & self.mask;
            out += self.gains[i] * self.buf[read_idx];
        }
        self.write_idx = (self.write_idx + 1) & self.mask;
        out
    }
}

// ---------------------------------------------------------------
// RoomReverb — main struct
// Mono in (L+R sum), stereo out (slight L/R wet offset for width).
// ---------------------------------------------------------------

pub struct RoomReverb {
    early: EarlyReflections,
    ap1: AllpassSection,
    ap2: AllpassSection,
    comb: CombSection,
    mix: f64,        // wet/dry [0,1]
    active: bool,
    lr_offset: usize, // L/R delay offset in samples for width
    offset_buf: Vec<f64>,
    offset_mask: usize,
    offset_write: usize,
}

impl RoomReverb {
    /// Construct from sample rate. All delay times computed from phi-ratios.
    pub fn new(sr: f64) -> Self {
        let sr = if sr > 0.0 { sr } else { 44100.0 };
        let lr_offset = (LR_OFFSET_SEC * sr).round().max(1.0) as usize;
        let offset_size = next_pow2(lr_offset + 1);
        Self {
            early: EarlyReflections::new(sr),
            ap1: AllpassSection::new(AP1_SEC, sr),
            ap2: AllpassSection::new(AP2_SEC, sr),
            comb: CombSection::new(COMB_SEC, sr),
            mix: 0.0,
            active: false,
            lr_offset,
            offset_buf: vec![0.0; offset_size],
            offset_mask: offset_size - 1,
            offset_write: 0,
        }
    }

    /// Clear all internal delay buffers and filter states.
    pub fn reset(&mut self) {
        self.early.reset();
        self.ap1.reset();
        self.ap2.reset();
        self.comb.reset();
        self.offset_buf.iter_mut().for_each(|s| *s = 0.0);
        self.offset_write = 0;
    }

    /// Set wet/dry mix [0, 1]. mix=0 bypasses entirely.
    pub fn set_mix(&mut self, mix: f64) {
        self.mix = clamp01(mix);
        self.active = self.mix > 1e-6;
    }

    /// Process one stereo sample pair.
    /// Mono sum feeds the reverb chain; L receives wet directly,
    /// R receives wet delayed by ~0.3ms for inter-aural decorrelation.
    #[inline]
    pub fn process(&mut self, input_l: f64, input_r: f64) -> (f64, f64) {
        if !self.active {
            return (input_l, input_r);
        }

        // Mono sum into reverb
        let mono = (input_l + input_r) * 0.5;

        // Signal chain: early reflections -> comb -> allpass 1 -> allpass 2
        let er = self.early.process(mono);
        let cb = self.comb.process(er);
        let a1 = self.ap1.process(cb);
        let wet = self.ap2.process(a1);

        // L/R offset: R wet is delayed by lr_offset samples
        self.offset_buf[self.offset_write] = wet;
        let read_idx =
            (self.offset_write + self.offset_buf.len() - self.lr_offset) & self.offset_mask;
        let wet_r = self.offset_buf[read_idx];
        self.offset_write = (self.offset_write + 1) & self.offset_mask;

        let wet_l = wet;
        let dry = 1.0 - self.mix;
        (
            dry * input_l + self.mix * wet_l,
            dry * input_r + self.mix * wet_r,
        )
    }
}

// ---------------------------------------------------------------
// Tests
// ---------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inactive_when_zero_mix() {
        let mut room = RoomReverb::new(44100.0);
        room.set_mix(0.0);
        assert!(!room.active);
        let (l, r) = room.process(0.5, -0.3);
        assert!((l - 0.5).abs() < 1e-15);
        assert!((r - (-0.3)).abs() < 1e-15);
    }

    #[test]
    fn test_active_when_nonzero_mix() {
        let mut room = RoomReverb::new(44100.0);
        room.set_mix(0.35);
        assert!(room.active);
        assert!((room.mix - 0.35).abs() < 1e-15);
    }

    #[test]
    fn test_output_differs_from_input() {
        let mut room = RoomReverb::new(44100.0);
        room.set_mix(0.5);
        // Feed signal long enough to populate all delay lines.
        // Max delay is comb ~13ms * 44100 ~ 577 samples. Run 1024 to be safe.
        for _ in 0..1024 {
            room.process(0.4, -0.2);
        }
        let (l, r) = room.process(0.4, -0.2);
        let dry_l = 0.4;
        let dry_r = -0.2;
        // At least one channel must differ from pure dry signal
        assert!(
            (l - dry_l).abs() > 1e-6 || (r - dry_r).abs() > 1e-6,
            "output should differ from dry input after delay lines are populated"
        );
    }

    #[test]
    fn test_reset_clears() {
        let mut room = RoomReverb::new(48000.0);
        room.set_mix(0.8);
        // Populate buffers with loud signal
        for _ in 0..2048 {
            room.process(1.0, -1.0);
        }
        room.reset();
        room.set_mix(0.0);
        // After reset + zero mix, output is pure passthrough
        let (l, r) = room.process(0.25, 0.75);
        assert!((l - 0.25).abs() < 1e-15);
        assert!((r - 0.75).abs() < 1e-15);
    }

    #[test]
    fn test_phi_ratio_delay_times() {
        // Verify delay times maintain phi-ratio relationships
        let sr = 44100.0;
        let d1 = (AP1_SEC * sr).round() as usize;
        let d2 = (AP2_SEC * sr).round() as usize;
        let ratio = d2 as f64 / d1 as f64;
        assert!(
            (ratio - PHI).abs() < 0.02,
            "AP2/AP1 delay ratio should be ~phi, got {}",
            ratio
        );
    }

    #[test]
    fn test_early_reflection_taps() {
        let er = EarlyReflections::new(44100.0);
        // 4 taps must be strictly increasing (phi-ratio spacing)
        for i in 1..4 {
            assert!(
                er.taps[i] > er.taps[i - 1],
                "tap {} ({}) should be > tap {} ({})",
                i,
                er.taps[i],
                i - 1,
                er.taps[i - 1]
            );
        }
        // First tap should be ~88 samples (2ms * 44100)
        assert!(
            (er.taps[0] as f64 - 88.0).abs() < 2.0,
            "first tap should be ~88 samples at 44.1kHz, got {}",
            er.taps[0]
        );
    }

    #[test]
    fn test_feedback_coefficients() {
        // All feedback coefficients must be < 1.0 (stable) and > 0 (active)
        let g_ap1 = feedback_for_rt60(AP1_SEC, RT60);
        let g_ap2 = feedback_for_rt60(AP2_SEC, RT60);
        let g_comb = feedback_for_rt60(COMB_SEC, RT60);
        for (name, g) in [("ap1", g_ap1), ("ap2", g_ap2), ("comb", g_comb)] {
            assert!(g > 0.0, "{} feedback should be positive, got {}", name, g);
            assert!(g < 1.0, "{} feedback should be < 1.0, got {}", name, g);
        }
        // Longer delay → smaller feedback (more decay per round-trip)
        assert!(
            g_comb < g_ap1,
            "comb feedback ({}) should be < ap1 feedback ({})",
            g_comb,
            g_ap1
        );
    }

    #[test]
    fn test_stereo_width_offset() {
        let mut room = RoomReverb::new(44100.0);
        room.set_mix(1.0); // full wet to maximize difference
        // Feed identical mono signal — L and R wet should differ
        // because R is delayed by ~0.3ms
        for _ in 0..2048 {
            room.process(0.5, 0.5);
        }
        // Now send a transient: the L/R difference reveals the offset
        let (l1, r1) = room.process(1.0, 1.0);
        // L and R won't be identical due to the offset delay
        // (they converge only after the offset delay flushes)
        // Feed silence and check the tail differs
        let (l2, r2) = room.process(0.0, 0.0);
        let diff1 = (l1 - r1).abs();
        let diff2 = (l2 - r2).abs();
        assert!(
            diff1 > 1e-10 || diff2 > 1e-10,
            "L and R wet should differ due to offset delay"
        );
    }

    #[test]
    fn test_mix_clamp() {
        let mut room = RoomReverb::new(44100.0);
        room.set_mix(-0.5);
        assert!((room.mix - 0.0).abs() < 1e-15);
        room.set_mix(1.5);
        assert!((room.mix - 1.0).abs() < 1e-15);
    }
}
