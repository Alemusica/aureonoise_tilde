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

/// Seizure risk range: 8-25 Hz auditory driving
/// (analogous to photic driving in photosensitive epilepsy).
/// Ref: Harding & Harding 1999, Fisher 2005 (Epilepsia).
const SEIZURE_RISK_LO: f64 = 8.0;
const SEIZURE_RISK_HI: f64 = 25.0;

/// Safety attenuation applied when rate is in seizure risk range.
/// -6 dB = multiply by 0.5. Hardware-level safety net.
const SEIZURE_ATTEN: f64 = 0.5;

/// Body resonance ranges (Hz):
/// - 5-8 Hz: thoracic cavity mechanical resonance
/// - 18-20 Hz: ocular globe / retinal resonance (~19 Hz peak)
/// Ref: von Gierke 1971 (AGARD), Gavreau 1966.
const BODY_RESO_LO_A: f64 = 5.0;
const BODY_RESO_HI_A: f64 = 8.0;
const BODY_RESO_LO_B: f64 = 18.0;
const BODY_RESO_HI_B: f64 = 20.0;

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

    /// Returns true when the given rate falls in the 8-25 Hz seizure risk range.
    /// Auditory driving at these frequencies can trigger seizures in susceptible
    /// individuals (analogous to photic driving).
    /// Ref: Harding & Harding 1999, Fisher 2005.
    #[inline]
    pub fn is_seizure_risk_range(rate_hz: f64) -> bool {
        rate_hz >= SEIZURE_RISK_LO && rate_hz <= SEIZURE_RISK_HI
    }

    /// Returns true when the given rate falls in body resonance ranges:
    /// 5-8 Hz (thoracic cavity) or 18-20 Hz (ocular globe).
    /// Concentrated acoustic energy at these frequencies can cause
    /// mechanical discomfort or harm.
    /// Ref: von Gierke 1971, Gavreau 1966.
    #[inline]
    pub fn is_body_resonance_range(rate_hz: f64) -> bool {
        (rate_hz >= BODY_RESO_LO_A && rate_hz <= BODY_RESO_HI_A)
            || (rate_hz >= BODY_RESO_LO_B && rate_hz <= BODY_RESO_HI_B)
    }

    /// Generate one mono sample.
    /// carrier_hz: tone frequency (150-180 Hz default)
    /// rate_hz: pulse rate (1-40 Hz, maps to brainwave bands)
    /// duty: pulse width fraction (0.3-0.7)
    /// level: amplitude (0-1)
    ///
    /// Safety: when rate_hz is in the seizure risk range (8-25 Hz),
    /// output amplitude is automatically attenuated by 6 dB (×0.5)
    /// as a hardware-level safety net independent of GUI warnings.
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

        // Advance phases (modulo for safety at extreme frequencies)
        if sr <= 0.0 { return 0.0; }
        self.carrier_phase += clamp(carrier_hz, 20.0, sr * 0.5) / sr;
        if self.carrier_phase >= 1.0 { self.carrier_phase %= 1.0; }
        self.pulse_phase += rate_hz / sr;
        if self.pulse_phase >= 1.0 { self.pulse_phase %= 1.0; }

        // Seizure risk attenuation: -6 dB when rate is in 8-25 Hz range
        let safety = if Self::is_seizure_risk_range(rate_hz) { SEIZURE_ATTEN } else { 1.0 };

        level * carrier * env * safety
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

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seizure_risk_boundaries() {
        // Below range
        assert!(!IsochronicTone::is_seizure_risk_range(7.9));
        // At boundaries (inclusive)
        assert!(IsochronicTone::is_seizure_risk_range(8.0));
        assert!(IsochronicTone::is_seizure_risk_range(25.0));
        // Inside range
        assert!(IsochronicTone::is_seizure_risk_range(15.0));
        // Above range
        assert!(!IsochronicTone::is_seizure_risk_range(25.1));
        // 40 Hz (MIT GENUS gamma) should be outside seizure risk
        assert!(!IsochronicTone::is_seizure_risk_range(40.0));
    }

    #[test]
    fn body_resonance_boundaries() {
        // Below first range
        assert!(!IsochronicTone::is_body_resonance_range(4.9));
        // First range: 5-8 Hz (thoracic)
        assert!(IsochronicTone::is_body_resonance_range(5.0));
        assert!(IsochronicTone::is_body_resonance_range(6.5));
        assert!(IsochronicTone::is_body_resonance_range(8.0));
        // Gap between ranges
        assert!(!IsochronicTone::is_body_resonance_range(10.0));
        assert!(!IsochronicTone::is_body_resonance_range(17.9));
        // Second range: 18-20 Hz (ocular)
        assert!(IsochronicTone::is_body_resonance_range(18.0));
        assert!(IsochronicTone::is_body_resonance_range(19.0));
        assert!(IsochronicTone::is_body_resonance_range(20.0));
        // Above second range
        assert!(!IsochronicTone::is_body_resonance_range(20.1));
    }

    #[test]
    fn seizure_range_attenuates_output() {
        let sr = 44100.0;
        let carrier_hz = 165.0;
        let duty = 0.5;
        let level = 1.0;

        // Collect peak amplitude at 15 Hz (seizure risk range)
        let mut iso_risk = IsochronicTone::new();
        let mut peak_risk: f64 = 0.0;
        for _ in 0..4410 {
            let s = iso_risk.process_sample(sr, carrier_hz, 15.0, duty, level);
            if s.abs() > peak_risk { peak_risk = s.abs(); }
        }

        // Collect peak amplitude at 40 Hz (outside seizure risk)
        let mut iso_safe = IsochronicTone::new();
        let mut peak_safe: f64 = 0.0;
        for _ in 0..4410 {
            let s = iso_safe.process_sample(sr, carrier_hz, 40.0, duty, level);
            if s.abs() > peak_safe { peak_safe = s.abs(); }
        }

        // Risk-range output should be ~half of safe-range output (-6 dB)
        assert!(
            peak_risk > 0.0,
            "seizure-risk signal should be non-zero"
        );
        assert!(
            peak_safe > 0.0,
            "safe signal should be non-zero"
        );
        let ratio = peak_risk / peak_safe;
        assert!(
            (ratio - 0.5).abs() < 0.05,
            "seizure-risk attenuation should be ~0.5 (-6 dB), got ratio {}",
            ratio
        );
    }

    #[test]
    fn no_attenuation_at_40hz() {
        // 40 Hz is the MIT GENUS gamma frequency — must NOT be attenuated
        let mut iso = IsochronicTone::new();
        let sr = 44100.0;
        let mut peak: f64 = 0.0;
        for _ in 0..4410 {
            let s = iso.process_sample(sr, 165.0, 40.0, 0.5, 1.0);
            if s.abs() > peak { peak = s.abs(); }
        }
        // With level=1.0, duty=0.5, carrier peak = 1.0, envelope peak = 1.0
        // so peak output should approach 1.0 (no attenuation)
        assert!(
            peak > 0.9,
            "40 Hz signal should NOT be attenuated, peak = {}",
            peak
        );
    }
}
