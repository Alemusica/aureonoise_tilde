//! aureonoise - Binaural beat generator
//! Two sine oscillators separated by beat frequency, mixed below noise floor.
//! Research: binaural beats must be SEPARATE from noise path (Wahbeh 2007).

use crate::math::*;
use crate::constants::*;

/// Binaural beat generator — stereo sine pair
pub struct BinauralBeat {
    phase_l: f64,
    phase_r: f64,
}

impl BinauralBeat {
    pub fn new() -> Self {
        Self {
            phase_l: 0.0,
            phase_r: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.phase_l = 0.0;
        self.phase_r = 0.0;
    }

    /// Generate one stereo sample of binaural beat.
    /// carrier_hz: center frequency (200-500 Hz typical)
    /// beat_hz: binaural beat frequency (0.5-40 Hz)
    /// level: amplitude (0-1), typically 0.05-0.15 (-20 to -16 dB)
    #[inline]
    pub fn process_sample(&mut self, sr: f64, carrier_hz: f64, beat_hz: f64, level: f64) -> (f64, f64) {
        if sr <= 0.0 { return (0.0, 0.0); }
        let carrier = clamp(carrier_hz, 20.0, sr * 0.5);
        let freq_l = carrier - beat_hz * 0.5;
        let freq_r = carrier + beat_hz * 0.5;

        let out_l = level * (TWO_PI * self.phase_l).sin();
        let out_r = level * (TWO_PI * self.phase_r).sin();

        self.phase_l += freq_l / sr;
        if self.phase_l >= 1.0 { self.phase_l %= 1.0; }
        self.phase_r += freq_r / sr;
        if self.phase_r >= 1.0 { self.phase_r %= 1.0; }

        (out_l, out_r)
    }
}

impl Default for BinauralBeat {
    fn default() -> Self {
        Self::new()
    }
}
