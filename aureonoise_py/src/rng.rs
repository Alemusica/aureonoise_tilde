//! aureonoise - Fast xorshift64 RNG
//! Optimized for audio DSP - not cryptographically secure

use pyo3::prelude::*;

/// Fast xorshift64 random number generator
#[pyclass]
#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
}

#[pymethods]
impl Rng {
    /// Create a new RNG with the given seed
    #[new]
    #[pyo3(signature = (seed = 0x9E3779B97F4A7C15))]
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x2545F4914F6CDD1D } else { seed },
        }
    }

    /// Reseed the RNG
    pub fn seed(&mut self, value: u64) {
        self.state = if value == 0 { 0x2545F4914F6CDD1D } else { value };
    }

    /// Generate next u64 value
    pub fn next_u64(&mut self) -> u64 {
        self.state ^= self.state >> 12;
        self.state ^= self.state << 25;
        self.state ^= self.state >> 27;
        self.state.wrapping_mul(2685821657736338717)
    }

    /// Generate uniform random in [0, 1)
    pub fn uni01(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9007199254740992.0)
    }

    /// Generate uniform random in [-1, 1)
    pub fn uni_pm1(&mut self) -> f64 {
        2.0 * self.uni01() - 1.0
    }

    /// Generate Gaussian-ish random using Box-Muller approximation
    /// (fast but not exact - good enough for audio)
    pub fn gauss_approx(&mut self) -> f64 {
        // Sum of 3 uniforms approximates Gaussian (central limit theorem)
        let u1 = self.uni_pm1();
        let u2 = self.uni_pm1();
        let u3 = self.uni_pm1();
        (u1 + u2 + u3) / 1.732050808 // normalize by sqrt(3)
    }
}

impl Default for Rng {
    fn default() -> Self {
        Self::new(0x9E3779B97F4A7C15)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rng_deterministic() {
        let mut rng1 = Rng::new(12345);
        let mut rng2 = Rng::new(12345);
        
        for _ in 0..100 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }

    #[test]
    fn test_uni01_range() {
        let mut rng = Rng::new(42);
        for _ in 0..10000 {
            let v = rng.uni01();
            assert!(v >= 0.0 && v < 1.0);
        }
    }

    #[test]
    fn test_uni_pm1_range() {
        let mut rng = Rng::new(42);
        for _ in 0..10000 {
            let v = rng.uni_pm1();
            assert!(v >= -1.0 && v < 1.0);
        }
    }
}
