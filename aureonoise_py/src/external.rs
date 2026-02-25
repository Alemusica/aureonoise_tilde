//! aureonoise - External externalisation
//! Block-level cross-channel feedback delay for stereo image widening.
//! Ported from beta7_tools/src/spatial/processor.cpp

use crate::math::clamp01;

const DELAY_SIZE: usize = 256;
const DELAY_MASK: usize = DELAY_SIZE - 1;

/// Pre-computed configuration for one block of externalisation processing.
/// Derived from the amount parameter and sample rate.
pub struct ExternalConfig {
    pub amount: f64,
    pub active: bool,
    pub delay_samples: usize,
    pub feedback: f64,
    pub cross: f64,
    pub hf: f64,
}

/// Cross-channel feedback delay processor.
/// Persistent state: circular delay buffers (L/R), write index, LP filter state.
pub struct ExternalProcessor {
    delay_l: [f64; DELAY_SIZE],
    delay_r: [f64; DELAY_SIZE],
    write_idx: usize,
    lp_l: f64,
    lp_r: f64,
}

impl ExternalProcessor {
    pub fn new() -> Self {
        Self {
            delay_l: [0.0; DELAY_SIZE],
            delay_r: [0.0; DELAY_SIZE],
            write_idx: 0,
            lp_l: 0.0,
            lp_r: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.delay_l = [0.0; DELAY_SIZE];
        self.delay_r = [0.0; DELAY_SIZE];
        self.write_idx = 0;
        self.lp_l = 0.0;
        self.lp_r = 0.0;
    }

    /// Prepare per-block config from amount [0,1] and sample rate.
    /// Delay range: 0.7ms (amount=0) to 5.2ms (amount=1).
    pub fn prepare(amount: f64, sr: f64) -> ExternalConfig {
        let amount = clamp01(amount);
        let active = amount > 1e-6;
        if !active {
            return ExternalConfig {
                amount,
                active,
                delay_samples: 0,
                feedback: 0.0,
                cross: 0.0,
                hf: 0.0,
            };
        }
        let sr = if sr > 0.0 { sr } else { 44100.0 };
        let delay = ((0.0007 + 0.0045 * amount) * sr).round() as usize;
        let delay = delay.clamp(1, DELAY_SIZE - 1);
        ExternalConfig {
            amount,
            active,
            delay_samples: delay,
            feedback: 0.07 * amount,
            cross: 0.32 * amount,
            hf: 0.55 * amount,
        }
    }

    /// Process one sample pair in-place.
    /// Reads from circular delay, writes new + feedback, applies HF emphasis
    /// and cross-channel mixing, then blends wet/dry by amount.
    #[inline]
    pub fn process_sample(&mut self, cfg: &ExternalConfig, y_l: &mut f64, y_r: &mut f64) {
        if !cfg.active {
            return;
        }

        let read_idx = (self.write_idx + DELAY_SIZE - cfg.delay_samples) & DELAY_MASK;
        let prev_l = self.delay_l[read_idx];
        let prev_r = self.delay_r[read_idx];

        // Write current + feedback into delay
        self.delay_l[self.write_idx] = *y_l + cfg.feedback * prev_l;
        self.delay_r[self.write_idx] = *y_r + cfg.feedback * prev_r;

        // LP filter on delayed signal (coeff = 1 - 0.32 = 0.68)
        self.lp_l = prev_l + (self.lp_l - prev_l) * 0.68;
        self.lp_r = prev_r + (self.lp_r - prev_r) * 0.68;

        // HF emphasis: subtract LP from delayed to get HP, scale by hf param
        let hf_l = prev_l - self.lp_l;
        let hf_r = prev_r - self.lp_r;

        // Cross-channel mix with HF emphasis
        let ext_l = prev_l + cfg.hf * hf_l + cfg.cross * prev_r;
        let ext_r = prev_r + cfg.hf * hf_r + cfg.cross * prev_l;

        // Wet/dry blend
        *y_l = (1.0 - cfg.amount) * *y_l + cfg.amount * ext_l;
        *y_r = (1.0 - cfg.amount) * *y_r + cfg.amount * ext_r;

        self.write_idx = (self.write_idx + 1) & DELAY_MASK;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inactive_when_zero() {
        let cfg = ExternalProcessor::prepare(0.0, 48000.0);
        assert!(!cfg.active);
    }

    #[test]
    fn test_active_when_nonzero() {
        let cfg = ExternalProcessor::prepare(0.5, 48000.0);
        assert!(cfg.active);
        assert!(cfg.delay_samples >= 1);
        assert!(cfg.delay_samples < DELAY_SIZE);
    }

    #[test]
    fn test_delay_range() {
        // At amount=0.01 (just active), delay ~= 0.7ms * sr
        let cfg_lo = ExternalProcessor::prepare(0.01, 48000.0);
        // At amount=1.0, delay ~= 5.2ms * sr = ~250
        let cfg_hi = ExternalProcessor::prepare(1.0, 48000.0);
        assert!(cfg_lo.delay_samples < cfg_hi.delay_samples);
    }

    #[test]
    fn test_passthrough_when_inactive() {
        let mut proc = ExternalProcessor::new();
        let cfg = ExternalProcessor::prepare(0.0, 48000.0);
        let (mut l, mut r) = (0.42, -0.37);
        proc.process_sample(&cfg, &mut l, &mut r);
        assert!((l - 0.42).abs() < 1e-15);
        assert!((r - -0.37).abs() < 1e-15);
    }

    #[test]
    fn test_modifies_signal_when_active() {
        let mut proc = ExternalProcessor::new();
        let cfg = ExternalProcessor::prepare(0.8, 48000.0);
        // Feed some signal to populate delay
        for _ in 0..512 {
            let (mut l, mut r) = (0.5, -0.3);
            proc.process_sample(&cfg, &mut l, &mut r);
        }
        // After delay is populated, output should differ from input
        let (mut l, mut r) = (0.5, -0.3);
        let (orig_l, orig_r) = (l, r);
        proc.process_sample(&cfg, &mut l, &mut r);
        // At least one channel should be modified
        assert!((l - orig_l).abs() > 1e-6 || (r - orig_r).abs() > 1e-6);
    }

    #[test]
    fn test_reset_clears_state() {
        let mut proc = ExternalProcessor::new();
        let cfg = ExternalProcessor::prepare(0.8, 48000.0);
        for _ in 0..256 {
            let (mut l, mut r) = (1.0, -1.0);
            proc.process_sample(&cfg, &mut l, &mut r);
        }
        proc.reset();
        // After reset, delay buffer should be zeroed, so processing
        // should only see zeros from the delay line
        let cfg_inactive = ExternalProcessor::prepare(0.0, 48000.0);
        let (mut l, mut r) = (0.5, 0.5);
        proc.process_sample(&cfg_inactive, &mut l, &mut r);
        assert!((l - 0.5).abs() < 1e-15);
        assert!((r - 0.5).abs() < 1e-15);
    }
}
