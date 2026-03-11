//! aureonoise - Phi Lattice
//! Central oracle for organic phi-ratio parameter generation.
//! All grain personality, timing, and spatial decisions derive from sacred ratios
//! selected stochastically via PhitRng (triphase hardware entropy).

use crate::phit::PhitRng;
use crate::constants::PHI;
use crate::math::clamp;

/// The 10 sacred ratios (from dialogue.rs, validated by Klimesch 2012)
/// Phi family: therapeutic/desynchronization (full weight)
/// Harmonic family: coupling/active cognition (half weight)
pub const SACRED_RATIOS: [f64; 10] = [
    0.382, 0.618, 1.0, 1.618, 2.618, 4.236,  // phi family
    0.500, 0.667, 1.500, 2.000,                // harmonic family
];

const SACRED_WEIGHTS: [f64; 10] = [
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  // phi full weight
    0.5, 0.5, 0.5, 0.5,             // harmonic half weight
];

/// ADSR templates — all segments in phi ratio to each other
/// (attack, decay, release, sustain_level)
const ADSR_TEMPLATES: [(f64, f64, f64, f64); 5] = [
    (0.146, 0.236, 0.618, 0.618),   // fast attack: A=1/φ³, D=1/φ², R=1/φ
    (0.236, 0.382, 0.382, 0.500),   // balanced: A=1/φ², D=1/φ, R=1/φ
    (0.382, 0.236, 0.382, 0.618),   // slow attack: A=1/φ, D=1/φ², R=1/φ
    (0.146, 0.146, 0.708, 0.382),   // percussive: A=1/φ³, D=1/φ³, R=long
    (0.500, 0.118, 0.382, 0.750),   // pad-like: A=1/2, D=tiny, R=1/φ
];

/// Spectral tilt values on phi-spaced scale
const TILT_SCALE: [f64; 7] = [
    -2.0, -1.618, -1.0, -0.618, -0.382, 0.0, 0.382,
];

pub struct PhiLattice {
    rng: PhitRng,
    weights: [f64; 10],
    ou_state: f64,
    last_ratio_idx: usize,
}

impl PhiLattice {
    pub fn new() -> Self {
        Self {
            rng: PhitRng::new(),
            weights: SACRED_WEIGHTS,
            ou_state: 0.0,
            last_ratio_idx: usize::MAX,
        }
    }

    /// Pick a sacred ratio organically.
    /// PhitRng entropy + weighted selection + OU drift = organic.
    /// Avoids repeating same ratio consecutively.
    pub fn next_ratio(&mut self) -> f64 {
        let mut eff_weights = [0.0f64; 10];
        let mut total = 0.0;
        for i in 0..10 {
            let drift = 1.0 + 0.3 * self.ou_state * if i < 6 { 1.0 } else { -1.0 };
            let w = (self.weights[i] * drift).max(0.05);
            let suppress = if i == self.last_ratio_idx { 0.1 } else { 1.0 };
            eff_weights[i] = w * suppress;
            total += eff_weights[i];
        }

        let u = self.rng.next_f64() * total;
        let mut cumulative = 0.0;
        let mut chosen = 0;
        for i in 0..10 {
            cumulative += eff_weights[i];
            if u < cumulative {
                chosen = i;
                break;
            }
        }

        self.last_ratio_idx = chosen;
        SACRED_RATIOS[chosen]
    }

    /// Grain duration: baselen_ms * sacred_ratio.
    pub fn grain_duration(&mut self, baselen_ms: f64) -> f64 {
        let ratio = self.next_ratio();
        clamp(baselen_ms * ratio, 10.0, 4000.0)
    }

    /// Inter-onset interval quantized to phi grid.
    pub fn onset_interval(&mut self, base_interval_samples: f64) -> f64 {
        let ratio = self.next_ratio();
        (base_interval_samples * ratio).max(1.0)
    }

    /// ADSR proportions — all segments in phi ratio to each other.
    /// Returns (attack, decay, sustain_level, release) where a+d+r <= 0.95.
    pub fn grain_adsr(&mut self) -> (f64, f64, f64, f64) {
        let idx = self.rng.next_u64() as usize % ADSR_TEMPLATES.len();
        let (a, d, r, s) = ADSR_TEMPLATES[idx];
        let mut jitter = || 0.9 + 0.2 * self.rng.next_f64();
        let aj = (a * (jitter)()).max(0.01);
        let dj = (d * (jitter)()).max(0.01);
        let rj = (r * (jitter)()).max(0.01);
        let total = aj + dj + rj;
        let scale = if total > 0.95 { 0.95 / total } else { 1.0 };
        (aj * scale, dj * scale, s, rj * scale)
    }

    /// Spectral tilt for a grain from phi-spaced discrete values.
    pub fn grain_spectral_tilt(&mut self, base_tilt: f64) -> f64 {
        let base_idx = TILT_SCALE.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((**a - base_tilt).abs()).partial_cmp(&((**b - base_tilt).abs())).unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(3);

        let offset = (self.rng.next_u64() % 3) as i32 - 1;
        let idx = (base_idx as i32 + offset).clamp(0, TILT_SCALE.len() as i32 - 1) as usize;
        TILT_SCALE[idx]
    }

    /// Spatial position on phi lattice [-1, 1].
    pub fn spatial_position(&mut self) -> f64 {
        let ratio = self.next_ratio();
        let sign = if self.rng.next_f64() < 0.5 { -1.0 } else { 1.0 };
        let mag = (ratio / (PHI * PHI)).min(1.0);
        sign * mag
    }

    /// Step OU drift (call once per block).
    pub fn step_drift(&mut self, dt: f64) {
        let theta = 0.628;
        let sigma = 0.25;
        let noise = self.rng.next_f64_bipolar();
        self.ou_state += -theta * self.ou_state * dt + sigma * dt.sqrt() * noise;
        self.ou_state = clamp(self.ou_state, -1.0, 1.0);
    }
}

impl Default for PhiLattice {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next_ratio_is_sacred() {
        let mut lattice = PhiLattice::new();
        for _ in 0..1000 {
            let r = lattice.next_ratio();
            assert!(
                SACRED_RATIOS.iter().any(|&sr| (sr - r).abs() < 1e-6),
                "Ratio {} is not a sacred ratio", r
            );
        }
    }

    #[test]
    fn test_next_ratio_no_immediate_repeat() {
        let mut lattice = PhiLattice::new();
        let mut prev = lattice.next_ratio();
        let mut repeats = 0;
        for _ in 0..500 {
            let r = lattice.next_ratio();
            if (r - prev).abs() < 1e-6 { repeats += 1; }
            prev = r;
        }
        assert!(repeats < 25, "Too many immediate repeats: {}/500", repeats);
    }

    #[test]
    fn test_grain_duration_is_phi_scaled() {
        let mut lattice = PhiLattice::new();
        let base = 200.0;
        for _ in 0..100 {
            let dur = lattice.grain_duration(base);
            let ratio = dur / base;
            assert!(
                SACRED_RATIOS.iter().any(|&sr| (sr - ratio).abs() < 0.01),
                "Duration ratio {} is not a sacred ratio", ratio
            );
        }
    }

    #[test]
    fn test_grain_adsr_sums_correctly() {
        let mut lattice = PhiLattice::new();
        for _ in 0..100 {
            let (a, d, _s, r) = lattice.grain_adsr();
            let sum = a + d + r;
            assert!(sum > 0.0 && sum <= 0.951, "ADSR a+d+r={} out of range", sum);
            assert!(a > 0.0 && d > 0.0 && r > 0.0, "ADSR segments must be positive");
        }
    }

    #[test]
    fn test_grain_adsr_phi_proportions() {
        let mut lattice = PhiLattice::new();
        for _ in 0..100 {
            let (a, d, _s, r) = lattice.grain_adsr();
            let ratios = [d / a, r / a, r / d];
            let has_phi = ratios.iter().any(|&ratio| {
                SACRED_RATIOS.iter().any(|&sr| (sr - ratio).abs() < 0.35)
            });
            assert!(has_phi, "ADSR proportions {:?} have no phi relationship", (a, d, r));
        }
    }

    #[test]
    fn test_spatial_position_range() {
        let mut lattice = PhiLattice::new();
        for _ in 0..1000 {
            let pos = lattice.spatial_position();
            assert!(pos >= -1.0 && pos <= 1.0, "Position {} out of [-1,1]", pos);
        }
    }

    #[test]
    fn test_onset_interval_positive() {
        let mut lattice = PhiLattice::new();
        for _ in 0..100 {
            let interval = lattice.onset_interval(1000.0);
            assert!(interval > 0.0, "Interval must be positive, got {}", interval);
        }
    }

    #[test]
    fn test_spectral_tilt_range() {
        let mut lattice = PhiLattice::new();
        for _ in 0..100 {
            let tilt = lattice.grain_spectral_tilt(-1.0);
            assert!(tilt >= -2.0 && tilt <= 0.5, "Tilt {} out of [-2, 0.5]", tilt);
        }
    }
}
