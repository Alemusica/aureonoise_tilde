//! aureonoise - Weyl sequences for low-discrepancy quasi-random numbers
//! Based on φ (golden ratio) and related irrationals

use pyo3::prelude::*;
use crate::constants::*;

/// Weyl sequence generator
/// Produces low-discrepancy sequences ideal for avoiding periodicity
#[pyclass]
#[derive(Clone, Debug)]
pub struct Weyl {
    /// Current state in [0, 1)
    x: f64,
    /// Step size (irrational number for best coverage)
    step: f64,
}

#[pymethods]
impl Weyl {
    /// Create a new Weyl sequence with given initial value and step
    #[new]
    #[pyo3(signature = (x = 0.5, step = None))]
    pub fn new(x: f64, step: Option<f64>) -> Self {
        Self {
            x: x.fract().abs(),
            step: step.unwrap_or(INV_PHI),
        }
    }

    /// Create a φ-based Weyl sequence (1/φ step)
    #[staticmethod]
    pub fn phi(x: f64) -> Self {
        Self::new(x, Some(INV_PHI))
    }

    /// Create a φ²-based Weyl sequence (1/φ² step)
    #[staticmethod]
    pub fn phi_sq(x: f64) -> Self {
        Self::new(x, Some(INV_PHI_SQ))
    }

    /// Create a φ³-based Weyl sequence (1/φ³ step)
    #[staticmethod]
    pub fn phi_cu(x: f64) -> Self {
        Self::new(x, Some(INV_PHI_CU))
    }

    /// Create a √2-based Weyl sequence (1/√2 step)
    #[staticmethod]
    pub fn sqrt2(x: f64) -> Self {
        Self::new(x, Some(INV_SQRT2))
    }

    /// Create a plastic constant-based Weyl sequence (1/ρ step)
    #[staticmethod]
    pub fn plastic(x: f64) -> Self {
        Self::new(x, Some(INV_PLASTIC))
    }

    /// Set the step size
    pub fn set_step(&mut self, step: f64) {
        self.step = step;
    }

    /// Get current value without advancing
    pub fn peek(&self) -> f64 {
        self.x
    }

    /// Advance and return new value in [0, 1)
    pub fn next(&mut self) -> f64 {
        self.x += self.step;
        self.x -= self.x.floor();
        self.x
    }

    /// Reset to initial value
    pub fn reset(&mut self, x: f64) {
        self.x = x.fract().abs();
    }
}

impl Default for Weyl {
    fn default() -> Self {
        Self::new(0.5, Some(INV_PHI))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weyl_range() {
        let mut w = Weyl::phi(0.123);
        for _ in 0..10000 {
            let v = w.next();
            assert!(v >= 0.0 && v < 1.0, "Value {} out of range", v);
        }
    }

    #[test]
    fn test_weyl_coverage() {
        // Weyl sequences should cover [0,1) fairly uniformly
        let mut w = Weyl::phi(0.0);
        let n = 1000;
        let mut bins = [0u32; 10];
        
        for _ in 0..n {
            let v = w.next();
            let bin = (v * 10.0).floor() as usize;
            if bin < 10 {
                bins[bin] += 1;
            }
        }
        
        // Each bin should have roughly n/10 = 100 samples
        for (i, &count) in bins.iter().enumerate() {
            assert!(count > 50 && count < 150, 
                "Bin {} has {} samples, expected ~100", i, count);
        }
    }
}
