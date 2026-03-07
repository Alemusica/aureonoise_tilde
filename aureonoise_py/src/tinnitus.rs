//! aureonoise - Tinnitus notch filter
//! 4th-order Butterworth notch (2 cascaded biquads).
//! Research: 10 Hz AM at tinnitus frequency → 19/28 patients showed suppression (p<0.0001).

use std::f64::consts::PI;
use crate::math::clamp;

/// Biquad filter state (Direct Form II Transposed)
#[derive(Clone, Debug)]
struct Biquad {
    b0: f64, b1: f64, b2: f64,
    a1: f64, a2: f64,
    z1: f64, z2: f64,
}

impl Biquad {
    fn new() -> Self {
        Self {
            b0: 1.0, b1: 0.0, b2: 0.0,
            a1: 0.0, a2: 0.0,
            z1: 0.0, z2: 0.0,
        }
    }

    /// Design notch filter at given frequency and Q
    fn design_notch(&mut self, freq_hz: f64, q: f64, sr: f64) {
        let w0 = 2.0 * PI * clamp(freq_hz, 20.0, sr * 0.45) / sr;
        let alpha = w0.sin() / (2.0 * q.max(0.1));
        let cos_w0 = w0.cos();

        let a0 = 1.0 + alpha;
        self.b0 = 1.0 / a0;
        self.b1 = -2.0 * cos_w0 / a0;
        self.b2 = 1.0 / a0;
        self.a1 = -2.0 * cos_w0 / a0;
        self.a2 = (1.0 - alpha) / a0;
    }

    fn reset(&mut self) {
        self.z1 = 0.0;
        self.z2 = 0.0;
    }

    #[inline]
    fn process(&mut self, x: f64) -> f64 {
        let y = self.b0 * x + self.z1;
        self.z1 = self.b1 * x - self.a1 * y + self.z2;
        self.z2 = self.b2 * x - self.a2 * y;
        y
    }
}

/// 4th-order Butterworth notch filter (2 cascaded biquads)
pub struct TinnitusNotch {
    stage1: Biquad,
    stage2: Biquad,
    freq_hz: f64,
    q: f64,
    sr: f64,
}

impl TinnitusNotch {
    pub fn new(sr: f64) -> Self {
        Self {
            stage1: Biquad::new(),
            stage2: Biquad::new(),
            freq_hz: 0.0,
            q: 6.0,
            sr,
        }
    }

    /// Configure notch. freq_hz=0 disables.
    /// Q controls bandwidth: higher Q = narrower notch.
    pub fn set_params(&mut self, freq_hz: f64, q: f64, sr: f64) {
        self.sr = sr;
        self.freq_hz = freq_hz;
        self.q = q.max(0.5);
        if freq_hz > 20.0 {
            // 4th-order Butterworth Q values
            let q1 = self.q * 0.5412; // 1/(2*cos(pi/8))
            let q2 = self.q * 1.3066; // 1/(2*cos(3*pi/8))
            self.stage1.design_notch(freq_hz, q1, sr);
            self.stage2.design_notch(freq_hz, q2, sr);
        }
    }

    pub fn reset(&mut self) {
        self.stage1.reset();
        self.stage2.reset();
    }

    /// Process one sample. Returns filtered sample.
    /// If freq_hz <= 20, passes through unchanged.
    #[inline]
    pub fn process(&mut self, x: f64) -> f64 {
        if self.freq_hz <= 20.0 {
            return x;
        }
        self.stage2.process(self.stage1.process(x))
    }
}

impl Default for TinnitusNotch {
    fn default() -> Self {
        Self::new(44100.0)
    }
}
