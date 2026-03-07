//! aureonoise - Isochronic tone generator
//! Sine carrier with amplitude modulation by Tukey-windowed pulse train.
//! Entrainment via monaural beats — no headphones required (unlike binaural).

use crate::math::*;
use crate::constants::*;

/// Isochronic tone — AM'd sine at target brainwave frequency
pub struct IsochronicTone {
    carrier_phase: f64,
    pulse_phase: f64,
}

impl IsochronicTone {
    pub fn new() -> Self {
        Self {
            carrier_phase: 0.0,
            pulse_phase: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.carrier_phase = 0.0;
        self.pulse_phase = 0.0;
    }

    /// Generate one mono sample.
    /// carrier_hz: tone frequency (150-180 Hz default)
    /// rate_hz: pulse rate (1-40 Hz, maps to brainwave bands)
    /// duty: pulse width fraction (0.3-0.7)
    /// level: amplitude (0-1)
    #[inline]
    pub fn process_sample(
        &mut self,
        sr: f64,
        carrier_hz: f64,
        rate_hz: f64,
        duty: f64,
        level: f64,
    ) -> f64 {
        let carrier = (TWO_PI * self.carrier_phase).sin();

        // Tukey-windowed pulse envelope
        let d = clamp(duty, 0.1, 0.9);
        let env = tukey_pulse(self.pulse_phase, d);

        // Advance phases
        self.carrier_phase += carrier_hz / sr;
        if self.carrier_phase >= 1.0 { self.carrier_phase -= 1.0; }
        self.pulse_phase += rate_hz / sr;
        if self.pulse_phase >= 1.0 { self.pulse_phase -= 1.0; }

        level * carrier * env
    }
}

/// Tukey-windowed pulse: smooth on/off transitions within duty cycle.
/// phase in [0,1), duty in (0,1).
/// Returns envelope in [0,1].
#[inline]
fn tukey_pulse(phase: f64, duty: f64) -> f64 {
    if phase >= duty {
        return 0.0;
    }
    // Tukey parameter: fraction of the ON period used for fade-in/out
    let alpha = 0.3;
    let fade = alpha * duty * 0.5;
    if fade < 1e-9 {
        return 1.0;
    }
    if phase < fade {
        // Fade in (raised cosine)
        0.5 * (1.0 - (PI * phase / fade).cos())
    } else if phase > duty - fade {
        // Fade out
        0.5 * (1.0 + (PI * (phase - (duty - fade)) / fade).cos())
    } else {
        1.0
    }
}

impl Default for IsochronicTone {
    fn default() -> Self {
        Self::new()
    }
}
