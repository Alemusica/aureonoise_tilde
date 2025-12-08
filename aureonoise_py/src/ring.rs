//! aureonoise - Ring buffer with fractional delay reading
//! Lagrange interpolation for smooth ITD

use pyo3::prelude::*;
use crate::constants::*;

/// Ring buffer for delay-based grain reading
#[pyclass]
#[derive(Clone)]
pub struct RingBuffer {
    data: Vec<f64>,
    write_index: usize,
}

#[pymethods]
impl RingBuffer {
    #[new]
    pub fn new() -> Self {
        Self {
            data: vec![0.0; RING_SIZE],
            write_index: 0,
        }
    }

    /// Clear the buffer
    pub fn clear(&mut self) {
        self.data.fill(0.0);
        self.write_index = 0;
    }

    /// Write a sample to the buffer
    #[inline]
    pub fn write(&mut self, sample: f64) {
        self.data[self.write_index] = sample;
        self.write_index = (self.write_index + 1) & RING_MASK;
    }

    /// Get current write index
    pub fn get_write_index(&self) -> usize {
        self.write_index
    }

    /// Read with linear interpolation at fractional position
    pub fn read_linear(&self, pos: f64) -> f64 {
        let i = pos.floor() as i64;
        let f = pos - i as f64;
        
        let i0 = (i as usize) & RING_MASK;
        let i1 = (i as usize + 1) & RING_MASK;
        
        self.data[i0] * (1.0 - f) + self.data[i1] * f
    }
}

impl Default for RingBuffer {
    fn default() -> Self {
        Self::new()
    }
}

impl RingBuffer {
    /// 4-point Lagrange interpolation for high-quality fractional delay
    #[inline]
    pub fn lagrange3(&self, pos: f64) -> f64 {
        let i = pos.floor() as i64;
        let f = pos - i as f64;
        
        let i_1 = ((i - 1) as usize) & RING_MASK;
        let i0 = (i as usize) & RING_MASK;
        let i1 = ((i + 1) as usize) & RING_MASK;
        let i2 = ((i + 2) as usize) & RING_MASK;
        
        let x_1 = self.data[i_1];
        let x0 = self.data[i0];
        let x1 = self.data[i1];
        let x2 = self.data[i2];
        
        let f1 = f - 1.0;
        let f2 = f - 2.0;
        
        let c_1 = -f * f1 * f2 / 6.0;
        let c0 = (f + 1.0) * f1 * f2 / 2.0;
        let c1 = -f * (f + 1.0) * f2 / 2.0;
        let c2 = f * (f + 1.0) * f1 / 6.0;
        
        c_1 * x_1 + c0 * x0 + c1 * x1 + c2 * x2
    }

    /// Read stereo with ITD (interaural time difference)
    /// Returns (left, right) samples
    #[inline]
    pub fn read_stereo_itd(&self, base: usize, itd: f64) -> (f64, f64) {
        let ad = itd.abs();
        let b = ad.ceil() + 2.0;
        
        let dl = b + if itd > 0.0 { itd } else { 0.0 };
        let dr = b + if itd < 0.0 { -itd } else { 0.0 };
        
        let pos_l = base as f64 - dl;
        let pos_r = base as f64 - dr;
        
        (self.lagrange3(pos_l), self.lagrange3(pos_r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ring_write_read() {
        let mut ring = RingBuffer::new();
        
        // Write some samples
        for i in 0..100 {
            ring.write(i as f64);
        }
        
        // Read back should work
        let (l, r) = ring.read_stereo_itd(ring.get_write_index().wrapping_sub(1), 0.0);
        assert!((l - r).abs() < 0.001); // Zero ITD = same sample
    }

    #[test]
    fn test_lagrange_interpolation() {
        let mut ring = RingBuffer::new();
        
        // Write a known pattern
        for i in 0..RING_SIZE {
            ring.write((i as f64 / RING_SIZE as f64).sin());
        }
        
        // Interpolation should be smooth
        let v1 = ring.lagrange3(100.0);
        let v2 = ring.lagrange3(100.5);
        let v3 = ring.lagrange3(101.0);
        
        // Values should be monotonic or at least continuous
        assert!((v2 - v1).abs() < 0.1);
        assert!((v3 - v2).abs() < 0.1);
    }
}
