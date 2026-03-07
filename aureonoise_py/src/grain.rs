//! aureonoise - Individual grain processor

use pyo3::prelude::*;
use crate::constants::*;
use crate::noise::GrainKind;
use crate::envelope::EnvelopeShape;

/// A single audio grain
#[pyclass]
#[derive(Clone)]
pub struct Grain {
    /// Is this grain active?
    pub on: bool,
    /// Current age in samples
    pub age: u32,
    /// Total duration in samples
    pub dur: u32,
    /// Amplitude
    pub amp: f64,
    /// Pan position [-1, 1]
    pub pan: f64,
    /// Left channel gain (equal power)
    pub pan_l: f64,
    /// Right channel gain (equal power)
    pub pan_r: f64,
    /// ITD in samples (fractional)
    pub itd: f64,
    /// ILD left gain
    pub gain_l: f64,
    /// ILD right gain
    pub gain_r: f64,
    /// Crossfeed amount
    pub crossfeed: f64,
    /// Focus (binaural localization sharpness)
    pub focus: f64,
    /// IPD allpass coefficient
    pub ipd_coeff: f64,
    /// IPD filter state L
    pub ipd_z_l: f64,
    /// IPD filter state R
    pub ipd_z_r: f64,
    /// Head shadow filter coefficient
    pub shadow_a: f64,
    /// Head shadow filter state L
    pub shadow_z_l: f64,
    /// Head shadow filter state R
    pub shadow_z_r: f64,
    /// Apply shadow to left channel
    pub shadow_left: bool,
    /// Apply shadow to right channel
    pub shadow_right: bool,
    /// Sample rate hold count
    pub sr_hold_n: i32,
    /// Sample rate hold counter
    pub sr_hold_cnt: i32,
    /// Held sample L
    pub held_l: f64,
    /// Held sample R
    pub held_r: f64,
    /// Quantization levels (0 = disabled)
    pub q_levels: i32,
    /// Grain effect kind
    #[pyo3(get, set)]
    pub kind: GrainKind,
    /// Envelope shape
    pub env: EnvelopeShape,
    /// Ring buffer read offset for inter-grain decorrelation.
    /// Each grain reads from a different region of the ring buffer,
    /// producing uncorrelated noise content across concurrent grains.
    pub ring_offset: usize,
}

#[pymethods]
impl Grain {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset grain to inactive state
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Check if grain is finished
    pub fn is_finished(&self) -> bool {
        !self.on || self.age >= self.dur
    }

    /// Get current phase [0, 1]
    pub fn phase(&self) -> f64 {
        if self.dur == 0 {
            1.0
        } else {
            self.age as f64 / self.dur as f64
        }
    }
}

impl Default for Grain {
    fn default() -> Self {
        Self {
            on: false,
            age: 0,
            dur: 0,
            amp: 0.0,
            pan: 0.0,
            pan_l: 1.0,
            pan_r: 1.0,
            itd: 0.0,
            gain_l: 1.0,
            gain_r: 1.0,
            crossfeed: 0.0,
            focus: 0.0,
            ipd_coeff: 0.0,
            ipd_z_l: 0.0,
            ipd_z_r: 0.0,
            shadow_a: 0.0,
            shadow_z_l: 0.0,
            shadow_z_r: 0.0,
            shadow_left: false,
            shadow_right: false,
            sr_hold_n: 1,
            sr_hold_cnt: 1,
            held_l: 0.0,
            held_r: 0.0,
            q_levels: 0,
            kind: GrainKind::Burst,
            env: EnvelopeShape::default(),
            ring_offset: 0,
        }
    }
}

/// Collection of grains
#[pyclass]
pub struct GrainPool {
    grains: Vec<Grain>,
}

#[pymethods]
impl GrainPool {
    #[new]
    pub fn new() -> Self {
        Self {
            grains: (0..MAX_GRAINS).map(|_| Grain::default()).collect(),
        }
    }

    /// Find a free grain slot, returns index or -1 if none available
    pub fn find_free(&self) -> i32 {
        for (i, g) in self.grains.iter().enumerate() {
            if !g.on {
                return i as i32;
            }
        }
        -1
    }

    /// Get number of active grains
    pub fn active_count(&self) -> usize {
        self.grains.iter().filter(|g| g.on).count()
    }

    /// Reset all grains
    pub fn reset_all(&mut self) {
        for g in &mut self.grains {
            g.reset();
        }
    }
}

impl Default for GrainPool {
    fn default() -> Self {
        Self::new()
    }
}

impl GrainPool {
    /// Get mutable reference to grain by index
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Grain> {
        self.grains.get_mut(index)
    }

    /// Get reference to grain by index
    pub fn get(&self, index: usize) -> Option<&Grain> {
        self.grains.get(index)
    }

    /// Iterate over all grains mutably
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Grain> {
        self.grains.iter_mut()
    }

    /// Iterate over all grains
    pub fn iter(&self) -> impl Iterator<Item = &Grain> {
        self.grains.iter()
    }
}
