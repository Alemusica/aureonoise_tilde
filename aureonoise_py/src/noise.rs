//! aureonoise - Colored noise generators
//! White, Pink, and Brown noise with φ-based filtering

use pyo3::prelude::*;
use crate::rng::Rng;
use crate::math::*;

/// Noise color type
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NoiseColor {
    White = 0,
    Pink = 1,
    Brown = 2,
}

/// Grain effect kind
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GrainKind {
    Burst = 0,
    VhsDrop = 1,
    Stutter = 2,
    Aliaser = 3,
}

impl GrainKind {
    /// Choose grain kind based on glitch mix and random value
    pub fn choose(mix: f64, u: f64) -> Self {
        let m = clamp01(mix);
        if u < 0.25 * m {
            GrainKind::VhsDrop
        } else if u < 0.60 * m {
            GrainKind::Stutter
        } else if u < 1.00 * m {
            GrainKind::Aliaser
        } else {
            GrainKind::Burst
        }
    }
}

/// Colored noise state
#[pyclass]
#[derive(Clone, Debug)]
pub struct NoiseColorState {
    /// Current noise color
    pub color: NoiseColor,
    /// Amount of coloring (0 = white, 1 = full color)
    pub amount: f64,
    // Filter states
    z1: f64,
    z2: f64,
    z3: f64,
}

#[pymethods]
impl NoiseColorState {
    #[new]
    #[pyo3(signature = (color = NoiseColor::Pink, amount = 0.65))]
    pub fn new(color: NoiseColor, amount: f64) -> Self {
        Self {
            color,
            amount: clamp01(amount),
            z1: 0.0,
            z2: 0.0,
            z3: 0.0,
        }
    }

    /// Reset filter states
    pub fn reset(&mut self) {
        self.z1 = 0.0;
        self.z2 = 0.0;
        self.z3 = 0.0;
    }

    /// Set noise color
    pub fn set_color(&mut self, color: NoiseColor) {
        self.color = color;
    }

    /// Set coloring amount
    pub fn set_amount(&mut self, amount: f64) {
        self.amount = clamp01(amount);
    }

    /// Process one sample of noise
    pub fn process(&mut self, rng: &mut Rng) -> f64 {
        let w = rng.uni_pm1();
        let amt = self.amount;
        
        if self.color == NoiseColor::White || amt <= 1.0e-6 {
            return w;
        }

        match self.color {
            NoiseColor::White => w,
            NoiseColor::Pink => {
                // Simple pink noise approximation
                self.z1 = (1.0 - 0.02 * amt) * self.z1 + (0.02 * amt) * w;
                (1.0 - amt) * w + amt * self.z1
            }
            NoiseColor::Brown => {
                // Brown noise with resonant character
                let slow_pole = map_phi_range(0.88, 0.995, amt);
                let slower_pole = map_phi_range(0.70, 0.985, amt);
                
                self.z1 = slow_pole * self.z1 + (1.0 - slow_pole) * w;
                self.z2 = slower_pole * self.z2 + (1.0 - slower_pole) * self.z1;
                self.z3 = (0.5 + 0.5 * (1.0 - amt)) * self.z3 + (0.5 * amt) * self.z2;
                
                let brown = self.z2 + 0.35 * self.z3;
                let shaped = soft_tanh(brown * (1.0 + 1.4 * amt));
                (1.0 - amt) * w + amt * shaped
            }
        }
    }
}

impl Default for NoiseColorState {
    fn default() -> Self {
        Self::new(NoiseColor::Pink, 0.65)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_noise_range() {
        let mut rng = Rng::new(42);
        let mut noise = NoiseColorState::new(NoiseColor::Pink, 0.8);
        
        for _ in 0..10000 {
            let v = noise.process(&mut rng);
            // Output should be roughly bounded
            assert!(v.abs() < 5.0, "Noise value {} out of expected range", v);
        }
    }

    #[test]
    fn test_grain_kind_choice() {
        // With mix = 0, always Burst
        assert_eq!(GrainKind::choose(0.0, 0.5), GrainKind::Burst);
        
        // With mix = 1, distribution varies
        assert_eq!(GrainKind::choose(1.0, 0.1), GrainKind::VhsDrop);
        assert_eq!(GrainKind::choose(1.0, 0.4), GrainKind::Stutter);
        assert_eq!(GrainKind::choose(1.0, 0.8), GrainKind::Aliaser);
    }
}
