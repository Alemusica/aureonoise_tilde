//! aureonoise - Stochastic processes
//! Ornstein-Uhlenbeck, Lattice dynamics, and Hawkes processes

use pyo3::prelude::*;
use crate::rng::Rng;
use crate::math::clamp01;
use crate::constants::*;

/// Ornstein-Uhlenbeck process for smooth stochastic modulation
#[pyclass]
#[derive(Clone, Debug)]
pub struct OrnsteinUhlenbeck {
    /// Current value
    pub y: f64,
    /// Time constant (relaxation rate)
    pub tau: f64,
    /// Volatility (noise amplitude)
    pub sigma: f64,
}

#[pymethods]
impl OrnsteinUhlenbeck {
    #[new]
    #[pyo3(signature = (tau = 0.5, sigma = 0.3))]
    pub fn new(tau: f64, sigma: f64) -> Self {
        Self {
            y: 0.0,
            tau: tau.max(1.0e-6),
            sigma,
        }
    }

    /// Step the process forward by dt seconds
    pub fn step(&mut self, dt: f64, mu: f64, rng: &mut Rng) -> f64 {
        let a = (-dt / self.tau.max(1.0e-6)).exp();
        let s = self.sigma * (1.0 - a * a).max(0.0).sqrt();
        self.y = a * self.y + (1.0 - a) * mu + s * rng.uni_pm1();
        self.y
    }

    /// Reset to initial state
    pub fn reset(&mut self) {
        self.y = 0.0;
    }

    /// Get current value
    pub fn value(&self) -> f64 {
        self.y
    }
}

impl Default for OrnsteinUhlenbeck {
    fn default() -> Self {
        Self::new(0.5, 0.3)
    }
}

/// 3D Lattice for complex spatial modulation
#[pyclass]
#[derive(Clone)]
pub struct Lattice {
    /// Grid dimensions
    pub x_dim: usize,
    pub y_dim: usize,
    pub z_dim: usize,
    /// Coupling strength (φ-corrected default: 1/φ³)
    pub eps: f64,
    /// Nonlinearity strength
    pub gamma: f64,
    /// Noise amplitude
    pub sigma: f64,
    /// Current state
    data: Vec<f64>,
    /// Temporary buffer for updates
    tmp: Vec<f64>,
}

#[pymethods]
impl Lattice {
    #[new]
    #[pyo3(signature = (x = 8, y = 8, z = 4))]
    pub fn new(x: usize, y: usize, z: usize) -> Self {
        let x = x.max(2);
        let y = y.max(2);
        let z = z.max(1);
        let n = x * y * z;
        Self {
            x_dim: x,
            y_dim: y,
            z_dim: z,
            eps: INV_PHI_CU, // φ-corrected: was 0.18
            gamma: PHI,      // φ-corrected: was 1.4
            sigma: 0.06,
            data: vec![0.0; n],
            tmp: vec![0.0; n],
        }
    }

    /// Resize the lattice
    pub fn resize(&mut self, x: usize, y: usize, z: usize) {
        let x = x.max(2);
        let y = y.max(2);
        let z = z.max(1);
        let n = x * y * z;
        self.x_dim = x;
        self.y_dim = y;
        self.z_dim = z;
        self.data = vec![0.0; n];
        self.tmp = vec![0.0; n];
    }

    /// Step the lattice dynamics
    pub fn step(&mut self, rng: &mut Rng) {
        let x_dim = self.x_dim;
        let y_dim = self.y_dim;
        let z_dim = self.z_dim;
        
        for k in 0..z_dim {
            for j in 0..y_dim {
                for i in 0..x_dim {
                    let p = self.idx(i, j, k);
                    
                    // 6-neighbor sum with periodic boundaries
                    let neighbors = 
                        self.act(self.data[self.idx(i.wrapping_add(1) % x_dim, j, k)]) +
                        self.act(self.data[self.idx(i.wrapping_sub(1) % x_dim, j, k)]) +
                        self.act(self.data[self.idx(i, j.wrapping_add(1) % y_dim, k)]) +
                        self.act(self.data[self.idx(i, j.wrapping_sub(1) % y_dim, k)]) +
                        self.act(self.data[self.idx(i, j, k.wrapping_add(1) % z_dim)]) +
                        self.act(self.data[self.idx(i, j, k.wrapping_sub(1) % z_dim)]);
                    
                    self.tmp[p] = (1.0 - self.eps) * self.act(self.data[p]) 
                        + (self.eps / 6.0) * neighbors 
                        + self.sigma * rng.uni_pm1();
                }
            }
        }
        
        std::mem::swap(&mut self.data, &mut self.tmp);
    }

    /// Probe the lattice at a position u in [0, 1]
    pub fn probe(&self, u: f64) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        let n = self.data.len();
        let p = (clamp01(u) * n as f64).floor() as usize % n;
        self.data[p]
    }

    /// Reset lattice to zero
    pub fn reset(&mut self) {
        self.data.fill(0.0);
        self.tmp.fill(0.0);
    }
}

impl Lattice {
    #[inline]
    fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        let i = i % self.x_dim;
        let j = j % self.y_dim;
        let k = k % self.z_dim;
        (k * self.y_dim + j) * self.x_dim + i
    }

    #[inline]
    fn act(&self, v: f64) -> f64 {
        (self.gamma * v).tanh()
    }
}

impl Default for Lattice {
    fn default() -> Self {
        Self::new(8, 8, 4)
    }
}

/// Hawkes process for burst clustering
#[pyclass]
#[derive(Clone, Debug)]
pub struct Hawkes {
    /// Current intensity
    pub lambda: f64,
    /// Base intensity
    pub base: f64,
    /// Decay rate
    pub beta: f64,
}

#[pymethods]
impl Hawkes {
    #[new]
    #[pyo3(signature = (base = None, beta = None))]
    pub fn new(base: Option<f64>, beta: Option<f64>) -> Self {
        Self {
            lambda: 0.0,
            // φ-corrected defaults
            base: base.unwrap_or(4.0 * PHI),      // was 4.0
            beta: beta.unwrap_or(30.0 / PHI),     // was 30.0, slower decay
        }
    }

    /// Tick the process, returns true if an event occurred
    pub fn tick(&mut self, dt: f64, rng: &mut Rng) -> bool {
        self.lambda = self.base + (self.lambda - self.base) * (-self.beta * dt).exp();
        let p = 1.0 - (-self.lambda * dt).exp();
        let event = rng.uni01() < p;
        if event {
            // φ-corrected boost (was 0.7)
            self.lambda += self.base * HAWKES_BOOST;
        }
        event
    }

    /// Reset to base intensity
    pub fn reset(&mut self) {
        self.lambda = self.base;
    }

    /// Get current intensity
    pub fn intensity(&self) -> f64 {
        self.lambda
    }
}

impl Default for Hawkes {
    fn default() -> Self {
        Self::new(None, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ou_mean_reversion() {
        let mut ou = OrnsteinUhlenbeck::new(0.1, 0.1);
        let mut rng = Rng::new(42);
        
        // Push far from mean
        ou.y = 10.0;
        
        // Should revert toward 0
        for _ in 0..1000 {
            ou.step(0.01, 0.0, &mut rng);
        }
        
        assert!(ou.y.abs() < 2.0, "OU should revert to mean, got {}", ou.y);
    }

    #[test]
    fn test_lattice_bounded() {
        let mut lat = Lattice::new(4, 4, 2);
        let mut rng = Rng::new(42);
        
        for _ in 0..1000 {
            lat.step(&mut rng);
        }
        
        // Values should be bounded by tanh
        for &v in &lat.data {
            assert!(v.abs() < 2.0, "Lattice value {} should be bounded", v);
        }
    }

    #[test]
    fn test_hawkes_bursts() {
        let mut hawkes = Hawkes::new(None, None);
        let mut rng = Rng::new(42);
        
        let mut events = 0;
        for _ in 0..10000 {
            if hawkes.tick(0.001, &mut rng) {
                events += 1;
            }
        }
        
        // Should have some events but not too many
        assert!(events > 10 && events < 1000, "Got {} events", events);
    }
}
