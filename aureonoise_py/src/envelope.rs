//! aureonoise - Envelope generator
//! Adaptive ADSR with φ-based dynamics

use pyo3::prelude::*;
use crate::math::*;

/// Envelope shape parameters
#[pyclass]
#[derive(Clone, Debug)]
pub struct EnvelopeShape {
    /// End of attack phase (0..1)
    pub attack_end: f64,
    /// End of decay phase (0..1)
    pub decay_end: f64,
    /// Start of release phase (0..1)
    pub release_start: f64,
    /// Sustain level (0..1)
    pub sustain_level: f64,
}

#[pymethods]
impl EnvelopeShape {
    #[new]
    #[pyo3(signature = (attack_end = 0.12, decay_end = 0.32, release_start = 0.82, sustain_level = 0.55))]
    pub fn new(attack_end: f64, decay_end: f64, release_start: f64, sustain_level: f64) -> Self {
        Self {
            attack_end: clamp(attack_end, 0.001, 0.75),
            decay_end: clamp(decay_end, attack_end + 0.001, 0.96),
            release_start: clamp(release_start, decay_end + 0.001, 0.999),
            sustain_level: clamp01(sustain_level),
        }
    }
}

impl Default for EnvelopeShape {
    fn default() -> Self {
        Self::new(0.12, 0.32, 0.82, 0.55)
    }
}

/// Envelope state and generator
#[pyclass]
#[derive(Clone, Debug)]
pub struct Envelope {
    /// Temperature for dynamics
    pub temperature: f64,
    // Tracking state for adaptive behavior
    last_gap_samples: f64,
    last_dur_samples: f64,
    last_center_distance: f64,
    last_io_ratio: f64,
}

#[pymethods]
impl Envelope {
    #[new]
    pub fn new() -> Self {
        Self {
            temperature: 0.45,
            last_gap_samples: 0.0,
            last_dur_samples: 0.0,
            last_center_distance: 0.0,
            last_io_ratio: 1.0,
        }
    }

    /// Create an adaptive envelope shape based on context
    pub fn make_shape(
        &mut self,
        attack_ratio: f64,
        decay_ratio: f64,
        sustain_level: f64,
        release_ratio: f64,
        gap_samples: f64,
        dur_samples: f64,
        center_distance: f64,
    ) -> EnvelopeShape {
        let sustain_level = clamp01(sustain_level);
        let attack_ratio = clamp(attack_ratio, 0.01, 4.0);
        let decay_ratio = clamp(decay_ratio, 0.01, 4.0);
        let release_ratio = clamp(release_ratio, 0.01, 4.0);
        let gap_samples = gap_samples.max(0.0);
        let dur_samples = dur_samples.max(1.0);
        let center_distance = clamp01(center_distance);

        // Update tracking
        self.last_gap_samples = gap_samples;
        self.last_dur_samples = dur_samples;
        self.last_center_distance = center_distance;

        let io_ratio = gap_samples / dur_samples;
        self.last_io_ratio = io_ratio;

        // Dynamic envelope shaping
        let io_weight = clamp(0.15 + 0.7 * (io_ratio / (io_ratio + 1.0)), 0.15, 0.95);
        let center_boost = 0.55 + 0.45 * (1.0 - center_distance);
        let dynamic_span = clamp(io_weight * center_boost, 0.15, 0.95);

        let total = attack_ratio + decay_ratio + release_ratio;
        let safe_total = if total <= 1.0e-9 { 1.0 } else { total };
        
        let mut attack_norm = (attack_ratio / safe_total) * dynamic_span;
        let mut decay_norm = (decay_ratio / safe_total) * dynamic_span;
        let mut release_norm = (release_ratio / safe_total) * dynamic_span;

        let sustain_norm = 1.0 - (attack_norm + decay_norm + release_norm);
        if sustain_norm < 0.05 {
            let deficit = 0.05 - sustain_norm;
            let sum_adr = attack_norm + decay_norm + release_norm;
            if sum_adr > 1.0e-6 {
                let scale = 1.0 - deficit / sum_adr;
                attack_norm *= scale;
                decay_norm *= scale;
                release_norm *= scale;
            }
        }

        let attack_end = clamp(attack_norm, 1.0e-4, 0.75);
        let decay_end = clamp(attack_end + decay_norm.max(1.0e-4), attack_end + 1.0e-4, 0.96);
        let rel_start = 1.0 - clamp(release_norm, 1.0e-4, 0.95);
        let release_start = clamp(rel_start.max(decay_end), decay_end + 1.0e-4, 0.999);

        EnvelopeShape {
            attack_end,
            decay_end,
            release_start,
            sustain_level,
        }
    }

    /// Evaluate envelope at phase [0, 1]
    #[staticmethod]
    pub fn eval(phase: f64, shape: &EnvelopeShape) -> f64 {
        let phase = clamp01(phase);
        
        if phase <= shape.attack_end {
            // Attack: linear ramp up
            phase / shape.attack_end.max(1.0e-6)
        } else if phase <= shape.decay_end {
            // Decay: linear ramp to sustain
            let t = (phase - shape.attack_end) / (shape.decay_end - shape.attack_end).max(1.0e-6);
            1.0 + (shape.sustain_level - 1.0) * t
        } else if phase < shape.release_start {
            // Sustain: constant level
            shape.sustain_level
        } else {
            // Release: linear ramp down
            let t = (phase - shape.release_start) / (1.0 - shape.release_start).max(1.0e-6);
            shape.sustain_level * (1.0 - t)
        }
    }
}

impl Default for Envelope {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_envelope_shape() {
        let shape = EnvelopeShape::default();
        
        // Start at 0
        let v = Envelope::eval(0.0, &shape);
        assert!(v.abs() < 0.01);
        
        // Peak at attack end
        let v = Envelope::eval(shape.attack_end, &shape);
        assert!((v - 1.0).abs() < 0.01);
        
        // Sustain in middle
        let v = Envelope::eval((shape.decay_end + shape.release_start) / 2.0, &shape);
        assert!((v - shape.sustain_level).abs() < 0.01);
        
        // End at 0
        let v = Envelope::eval(1.0, &shape);
        assert!(v.abs() < 0.01);
    }
}
