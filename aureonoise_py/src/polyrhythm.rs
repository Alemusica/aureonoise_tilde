//! aureonoise - Polyrhythm clock
//! Two pulse trains at ratio p:q per hemisphere.
//! Coincidence = structural handshake opportunity.

use crate::math::clamp;

/// Polyrhythm clock — two pulse streams at p:q ratio
pub struct PolyrhythmClock {
    phase_p: f64,
    phase_q: f64,
    /// Left hemisphere ratio numerator
    pub p: u32,
    /// Right hemisphere ratio denominator
    pub q: u32,
    /// Base rate in Hz (both streams scaled by this)
    pub base_rate: f64,
    /// Active flag
    pub active: bool,
    // Output state
    p_pulse: bool,
    q_pulse: bool,
    coincidence: bool,
    /// Tolerance window for coincidence detection (fraction of shorter period)
    tolerance: f64,
}

impl PolyrhythmClock {
    pub fn new() -> Self {
        Self {
            phase_p: 0.0,
            phase_q: 0.0,
            p: 3,
            q: 2,
            base_rate: 0.5,
            active: false,
            p_pulse: false,
            q_pulse: false,
            coincidence: false,
            tolerance: 0.08,
        }
    }

    pub fn reset(&mut self) {
        self.phase_p = 0.0;
        self.phase_q = 0.0;
        self.p_pulse = false;
        self.q_pulse = false;
        self.coincidence = false;
    }

    /// Advance clock by elapsed_samples. Returns (p_pulse, q_pulse, coincidence).
    pub fn tick(&mut self, sr: f64, elapsed_samples: f64) -> (bool, bool, bool) {
        if !self.active || sr <= 0.0 || self.base_rate <= 0.0 {
            return (false, false, false);
        }

        let dt = elapsed_samples / sr;
        let p_rate = self.p as f64 * self.base_rate;
        let q_rate = self.q as f64 * self.base_rate;

        let old_p = self.phase_p;
        let old_q = self.phase_q;

        self.phase_p += p_rate * dt;
        self.phase_q += q_rate * dt;

        // Detect wraps (pulse on phase crossing integer boundary)
        self.p_pulse = self.phase_p.floor() > old_p.floor();
        self.q_pulse = self.phase_q.floor() > old_q.floor();

        // Coincidence: both pulses within tolerance window
        self.coincidence = if self.p_pulse && self.q_pulse {
            true
        } else if self.p_pulse || self.q_pulse {
            // Check if the other phase is near an integer boundary
            let tol = self.tolerance;
            if self.p_pulse {
                let q_frac = self.phase_q.fract();
                q_frac < tol || q_frac > (1.0 - tol)
            } else {
                let p_frac = self.phase_p.fract();
                p_frac < tol || p_frac > (1.0 - tol)
            }
        } else {
            false
        };

        // Wrap phases to prevent accumulation
        if self.phase_p >= 1000.0 { self.phase_p -= self.phase_p.floor(); }
        if self.phase_q >= 1000.0 { self.phase_q -= self.phase_q.floor(); }

        (self.p_pulse, self.q_pulse, self.coincidence)
    }

    /// Get pan modulation from polyrhythm:
    /// p_pulse pushes left (-1), q_pulse pushes right (+1),
    /// coincidence → center (0). Returns pan offset in [-1, 1].
    pub fn pan_offset(&self, amount: f64) -> f64 {
        let amt = clamp(amount, 0.0, 1.0);
        if self.coincidence {
            0.0 // convergence point
        } else if self.p_pulse {
            -amt
        } else if self.q_pulse {
            amt
        } else {
            0.0
        }
    }

    /// Get pan modulation with asymmetric phi-lattice distribution.
    /// p_pulse distributes grains in left hemisphere using phi position,
    /// q_pulse distributes in right hemisphere.
    /// `phi_pos`: position from PhiLattice [0,1] for within-hemisphere spread.
    pub fn pan_offset_asymmetric(&self, amount: f64, phi_pos: f64) -> f64 {
        let amt = clamp(amount, 0.0, 1.0);
        if self.coincidence {
            0.0
        } else if self.p_pulse {
            -amt * (0.3 + 0.7 * clamp(phi_pos, 0.0, 1.0))
        } else if self.q_pulse {
            amt * (0.3 + 0.7 * clamp(phi_pos, 0.0, 1.0))
        } else {
            0.0
        }
    }

    /// Was there a coincidence on the last tick?
    pub fn had_coincidence(&self) -> bool {
        self.coincidence
    }
}

impl Default for PolyrhythmClock {
    fn default() -> Self {
        Self::new()
    }
}
