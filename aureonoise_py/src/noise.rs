//! aureonoise - Colored noise generators
//! White, Pink, Brown noise with φ-based filtering + Aureo/Quantum/Velvet modes

use pyo3::prelude::*;
use crate::rng::Rng;
use crate::weyl::Weyl;
use crate::math::*;
use crate::constants::*;

// ─── Enums ───────────────────────────────────────────────────────────────────

/// Noise color type (classic White/Pink/Brown)
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NoiseColor {
    White = 0,
    Pink = 1,
    Brown = 2,
}

/// Noise generation mode (classic + harmonic/sparse modes)
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NoiseMode {
    White = 0,
    Pink = 1,
    Brown = 2,
    Aureo = 3,
    Quantum = 4,
    Velvet = 5,
}

/// Grain effect kind
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GrainKind {
    Burst = 0,
    Stutter = 2,
    Aliaser = 3,
}

impl GrainKind {
    /// Choose grain kind based on glitch mix and random value
    pub fn choose(mix: f64, u: f64) -> Self {
        let m = clamp01(mix);
        if u < 0.40 * m {
            GrainKind::Stutter
        } else if u < 1.00 * m {
            GrainKind::Aliaser
        } else {
            GrainKind::Burst
        }
    }
}

// ─── Pink filter ─────────────────────────────────────────────────────────────

/// Paul Kellet's 6-stage pink noise filter
/// Produces 1/f spectral slope (±0.5dB) across 20Hz-20kHz
#[derive(Clone, Debug)]
pub struct PinkFilter {
    b: [f64; 7],
}

impl PinkFilter {
    pub fn new() -> Self {
        Self { b: [0.0; 7] }
    }

    pub fn reset(&mut self) {
        self.b = [0.0; 7];
    }

    /// Paul Kellet's refined method: 1/f ±0.5dB over 20Hz-20kHz
    pub fn process(&mut self, white: f64) -> f64 {
        self.b[0] = 0.99886 * self.b[0] + white * 0.0555179;
        self.b[1] = 0.99332 * self.b[1] + white * 0.0750759;
        self.b[2] = 0.96900 * self.b[2] + white * 0.1538520;
        self.b[3] = 0.86650 * self.b[3] + white * 0.3104856;
        self.b[4] = 0.55000 * self.b[4] + white * 0.5329522;
        self.b[5] = -0.7616 * self.b[5] - white * 0.0168980;
        let pink = self.b[0] + self.b[1] + self.b[2] + self.b[3]
                 + self.b[4] + self.b[5] + self.b[6] + white * 0.5362;
        self.b[6] = white * 0.115926;
        pink * 0.11  // normalize
    }
}

// ─── Classic colored noise (White/Pink/Brown) ────────────────────────────────

/// Colored noise state
#[pyclass]
#[derive(Clone, Debug)]
pub struct NoiseColorState {
    /// Current noise color
    pub color: NoiseColor,
    /// Amount of coloring (0 = white, 1 = full color)
    pub amount: f64,
    // Pink filter (Kellet 6-stage)
    pink: PinkFilter,
    // Brown: single-pole leaky integrator state
    z1: f64,
}

#[pymethods]
impl NoiseColorState {
    #[new]
    #[pyo3(signature = (color = NoiseColor::Pink, amount = 0.65))]
    pub fn new(color: NoiseColor, amount: f64) -> Self {
        Self {
            color,
            amount: clamp01(amount),
            pink: PinkFilter::new(),
            z1: 0.0,
        }
    }

    /// Reset filter states
    pub fn reset(&mut self) {
        self.pink.reset();
        self.z1 = 0.0;
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
                // Paul Kellet 6-stage: true 1/f ±0.5dB over 20Hz-20kHz
                (1.0 - amt) * w + amt * self.pink.process(w)
            }
            NoiseColor::Brown => {
                // Brown noise: single-pole leaky integrator (~21 Hz cutoff at 44.1kHz).
                // Pure integrator output (no white mix) gives PSD ∝ 1/f² = -6 dB/oct.
                // Normalization 0.04 keeps peaks within soft_tanh linear range.
                self.z1 = 0.997 * self.z1 + w;
                let brown = self.z1 * 0.04;
                (1.0 - amt) * w + amt * brown
            }
        }
    }
}

impl Default for NoiseColorState {
    fn default() -> Self {
        Self::new(NoiseColor::Pink, 0.65)
    }
}

// ─── Continuous spectral tilt ────────────────────────────────────────────────

/// Continuous spectral slope filter.
/// noise_slope: 0.0=White, -1.0=Pink, -2.0=Brown (continuous interpolation).
/// Extends to +0.5 (brightened white) via gentle highpass.
///
/// Pink path: Paul Kellet 6-stage IIR → PSD ∝ 1/f (beta=1, -3 dB/oct).
/// Brown path: single-pole leaky integrator (fc ~21 Hz at 44.1kHz),
/// pure output (no white mix) → PSD ∝ 1/f² (beta=2, -6 dB/oct).
/// Normalization 0.04 keeps peaks in soft_tanh linear range.
#[derive(Clone, Debug)]
pub struct SpectralTilt {
    pink: PinkFilter,
    /// Single-pole leaky integrator state for brown path
    z1: f64,
}

impl SpectralTilt {
    pub fn new() -> Self {
        Self {
            pink: PinkFilter::new(),
            z1: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.pink.reset();
        self.z1 = 0.0;
    }

    /// Process white noise sample through continuous spectral tilt.
    /// slope in [-2.0, +0.5]. Always runs both pink and brown to keep
    /// filter states warm, then crossfades.
    #[inline]
    pub fn process(&mut self, white: f64, slope: f64) -> f64 {
        let pink = self.pink.process(white);

        // Brown: single-pole leaky integrator (fc ~21 Hz at 44.1kHz).
        // Pure output gives PSD ∝ 1/f² (beta=2, -6 dB/oct).
        self.z1 = 0.997 * self.z1 + white;
        let brown = self.z1 * 0.04;

        let s = clamp(slope, -2.0, 0.5);

        if s >= 0.0 {
            // [0, +0.5]: white with gentle high-shelf boost
            // Simple highpass emphasis: white + s * (white - pink)
            white + s * (white - pink)
        } else if s >= -1.0 {
            // [0, -1]: white → pink
            let t = -s;
            (1.0 - t) * white + t * pink
        } else {
            // [-1, -2]: pink → brown
            let t = -s - 1.0;
            (1.0 - t) * pink + t * brown
        }
    }
}

impl Default for SpectralTilt {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Aureo harmonic state ────────────────────────────────────────────────────

/// Maximum partials for Aureo/Quantum harmonic stacks
pub const MAX_PARTIALS: usize = 32;

/// Size of the golden-relation wash pool
const RELATION_POOL: usize = 64;

/// Rational frequency ratios for the harmonic stack
const RATIONAL_SET: [f64; 10] = [
    1.0, 1.5, 1.25, 4.0 / 3.0, 5.0 / 3.0,
    1.75, 1.2, 1.4, 2.25, 15.0 / 8.0,
];

/// Irrational frequency ratios (phi, pi, sqrt-based)
const IRRATIO_SET: [f64; 5] = [
    SQRT2,                      // sqrt(2)
    PHI,                        // phi
    PI * 0.5,                   // pi/2
    1.2720196495140689,         // sqrt(phi)
    1.7724538509055159,         // sqrt(pi)
];

/// Single partial in the Aureo harmonic stack
#[derive(Clone, Debug)]
struct AureoPartial {
    ratio: f64,
    phase_phi: f64,
    phase_pi: f64,
    amp: f64,
}

impl Default for AureoPartial {
    fn default() -> Self {
        Self { ratio: 1.0, phase_phi: 0.0, phase_pi: 0.0, amp: 0.0 }
    }
}

/// Aureo harmonic noise state
/// Ported from C++ NoiseAureoState: harmonic stack with phi/pi ratios,
/// Planck-curve decay, prime-number stride, and golden-relation wash pool.
#[derive(Clone, Debug)]
pub struct NoiseAureoState {
    pub base_freq: f64,
    pub phi_pi_mix: f64,
    pub planck_decay: f64,
    pub prime_stride: f64,
    pub active_partials: i32,
    pub velvet_on: bool,
    pub velvet_amount: f64,

    // Velvet-only smoothing
    velvet_a: f64,
    velvet_z: f64,

    // Partials
    partials: Vec<AureoPartial>,

    // Golden-relation wash pool
    rel_pool: [f64; RELATION_POOL],
    rel_mean: f64,
    rel_ready: bool,
    rel_cursor: usize,
    rel_stride: usize,
    rel_weyl: Weyl,
}

impl NoiseAureoState {
    pub fn new() -> Self {
        let mut s = Self {
            base_freq: 220.0,
            phi_pi_mix: 0.5,
            planck_decay: 0.3,
            prime_stride: 1.0,
            active_partials: 12,
            velvet_on: true,
            velvet_amount: 0.10,
            velvet_a: 0.985,
            velvet_z: 0.0,
            partials: Vec::with_capacity(MAX_PARTIALS),
            rel_pool: [1.0; RELATION_POOL],
            rel_mean: 1.0,
            rel_ready: false,
            rel_cursor: 0,
            rel_stride: 5,
            rel_weyl: Weyl::new(0.5, Some(1.0 / (RELATION_POOL as f64 * PHI))),
        };
        s.setup_ratios();
        s
    }

    /// Initialize partial frequency ratios (rational + irrational interleaving)
    fn setup_ratios(&mut self) {
        self.partials.clear();
        let r_count = RATIONAL_SET.len();
        let i_count = IRRATIO_SET.len();
        for idx in 0..MAX_PARTIALS {
            let mut ratio = RATIONAL_SET[idx % r_count];
            if (idx / r_count) % 2 == 1 {
                ratio *= IRRATIO_SET[(idx / r_count) % i_count];
            }
            self.partials.push(AureoPartial {
                ratio,
                phase_phi: 0.0,
                phase_pi: 0.0,
                amp: 0.0,
            });
        }
    }

    /// Build the 64-layer golden-relation wash pool
    fn rebuild_relations(&mut self, rng: &mut Rng) {
        let jitter = 0.001 * rng.uni_pm1();
        let mut sum = 0.0;
        for i in 0..RELATION_POOL {
            let u = ((i as f64 + 0.5) * INV_PHI + jitter).fract().abs();
            let v = ((i as f64 + 0.5) * INV_PLASTIC + 0.5 * jitter).fract().abs();
            let th1 = TWO_PI * u;
            let th2 = TWO_PI * v;
            let mut w = 1.0
                + 0.06 * th1.sin()
                + 0.04 * th2.cos()
                + 0.02 * (th1 + th2).sin();
            w *= PHI.powf(0.03 * th2.sin());
            let pmod = (i % 7) as f64 - 3.0;
            w *= 1.0 + 0.01 * pmod;
            self.rel_pool[i] = w;
            sum += w;
        }
        // Normalize to mean=1, clamp to [0.9, 1.1]
        self.rel_mean = if sum > 1.0e-9 { sum / RELATION_POOL as f64 } else { 1.0 };
        let mut new_sum = 0.0;
        for i in 0..RELATION_POOL {
            let w = clamp(self.rel_pool[i] / self.rel_mean, 0.90, 1.10);
            self.rel_pool[i] = w;
            new_sum += w;
        }
        self.rel_mean = new_sum / RELATION_POOL as f64;
        if self.rel_mean <= 1.0e-9 { self.rel_mean = 1.0; }

        self.rel_cursor = ((clamp01(rng.uni01()) * (RELATION_POOL - 1) as f64).floor()) as usize;
        static PRIMES: [usize; 5] = [3, 5, 7, 11, 13];
        self.rel_stride = PRIMES[(self.rel_cursor + 2) % 5];
        self.rel_weyl.reset(clamp01(rng.uni01()));
        self.rel_weyl.set_step(1.0 / (RELATION_POOL as f64 * PHI));
        self.rel_ready = true;
    }

    /// Read from the wash pool with triangular smoothing
    #[inline]
    fn pool_smooth(&self, idx: usize) -> f64 {
        let n = RELATION_POOL;
        let w0 = self.rel_pool[(idx + n - 1) % n];
        let w1 = self.rel_pool[idx % n];
        let w2 = self.rel_pool[(idx + 1) % n];
        0.25 * w0 + 0.50 * w1 + 0.25 * w2
    }

    pub fn reset(&mut self) {
        self.setup_ratios();
        self.rel_ready = false;
        self.rel_mean = 1.0;
        self.rel_cursor = 0;
        self.rel_stride = 5;
        self.velvet_z = 0.0;
    }

    /// Generate one sample (Aureo/Quantum shared path)
    pub fn process(&mut self, rng: &mut Rng, sr: f64) -> f64 {
        if !self.rel_ready {
            self.rebuild_relations(rng);
        }

        // Advance the relation-pool cursor via Weyl sequence
        let u = self.rel_weyl.next();
        let next_cursor = (u * RELATION_POOL as f64).floor() as usize % RELATION_POOL;
        if next_cursor != self.rel_cursor {
            self.rel_cursor = next_cursor;
        }

        let base_comp = if self.rel_mean > 1.0e-12 { 1.0 / self.rel_mean } else { 1.0 };
        let n = (self.active_partials as usize).min(self.partials.len());
        let mut real_acc = 0.0;
        let mut imag_acc = 0.0;

        // Pre-compute pool weights to avoid borrow conflict in the loop
        let base_freq = self.base_freq;
        let prime_stride = self.prime_stride;
        let planck_decay = self.planck_decay;
        let rel_cursor = self.rel_cursor;
        let rel_stride = self.rel_stride;

        for i in 0..n {
            // Pool lookup outside mutable borrow
            let ridx = rel_cursor + i * rel_stride;
            let rw = self.pool_smooth(ridx);

            let p = &mut self.partials[i];
            let freq = base_freq * base_comp * p.ratio * rw;
            let omega = TWO_PI * freq / sr;
            p.phase_phi += omega * prime_stride;
            p.phase_pi += omega * (1.0 + prime_stride * 0.23);
            if p.phase_phi > TWO_PI { p.phase_phi = p.phase_phi % TWO_PI; }
            if p.phase_pi > TWO_PI { p.phase_pi = p.phase_pi % TWO_PI; }

            // Planck-curve amplitude decay
            let target = (-planck_decay * i as f64 / n.max(1) as f64).exp();
            p.amp += (target - p.amp) * 0.02;

            real_acc += p.amp * p.phase_phi.cos();
            imag_acc += p.amp * p.phase_pi.cos();
        }

        let mix = clamp01(self.phi_pi_mix);
        let mut out = (1.0 - mix) * real_acc + mix * imag_acc;

        if self.velvet_on && self.velvet_amount > 1.0e-6 {
            out += self.velvet_amount * rng.uni_pm1();
        }
        out
    }

    /// Velvet-only mode: lightly smoothed sparse noise
    pub fn process_velvet_only(&mut self, rng: &mut Rng) -> f64 {
        let w = rng.uni_pm1();
        let y = (1.0 - self.velvet_a) * w + self.velvet_a * self.velvet_z;
        self.velvet_z = y;
        self.velvet_amount * y
    }
}

// ─── NoiseGen: unified dispatch ──────────────────────────────────────────────

/// Unified noise generator that dispatches across all 6 modes.
/// Classic modes (White/Pink/Brown) delegate to NoiseColorState.
/// Aureo/Quantum use the harmonic stack. Velvet uses sparse impulses.
#[pyclass]
#[derive(Clone, Debug)]
pub struct NoiseGen {
    pub mode: NoiseMode,
    pub classic: NoiseColorState,
    pub aureo: NoiseAureoState,

    // Aureo parameters (public for configuration from lib.rs)
    pub aureo_decay: f64,
    pub aureo_stride: f64,
    pub aureo_harmonics: i32,

    // Quantum parameters
    pub quantum_detail: f64,
    pub quantum_base: f64,

    // Velvet parameters
    pub velvet_density: f64,
}

#[pymethods]
impl NoiseGen {
    #[new]
    #[pyo3(signature = (mode = NoiseMode::Pink))]
    pub fn new(mode: NoiseMode) -> Self {
        Self {
            mode,
            classic: NoiseColorState::default(),
            aureo: NoiseAureoState::new(),
            aureo_decay: 0.3,
            aureo_stride: 1.0,
            aureo_harmonics: 12,
            quantum_detail: 0.7,
            quantum_base: 220.0,
            velvet_density: 2000.0,
        }
    }

    /// Set the active noise mode
    pub fn set_mode(&mut self, mode: NoiseMode) {
        self.mode = mode;
        // Map classic modes to NoiseColor for the classic state
        match mode {
            NoiseMode::White => self.classic.set_color(NoiseColor::White),
            NoiseMode::Pink  => self.classic.set_color(NoiseColor::Pink),
            NoiseMode::Brown => self.classic.set_color(NoiseColor::Brown),
            _ => {}
        }
    }

    /// Set Aureo decay (Planck-curve steepness)
    pub fn set_aureo_decay(&mut self, decay: f64) {
        self.aureo_decay = clamp01(decay);
        self.aureo.planck_decay = self.aureo_decay;
    }

    /// Set Aureo prime stride
    pub fn set_aureo_stride(&mut self, stride: f64) {
        self.aureo_stride = stride.max(0.1);
        self.aureo.prime_stride = self.aureo_stride;
    }

    /// Set Aureo harmonic count (1-32)
    pub fn set_aureo_harmonics(&mut self, n: i32) {
        self.aureo_harmonics = n.clamp(1, MAX_PARTIALS as i32);
        self.aureo.active_partials = self.aureo_harmonics;
    }

    /// Set Quantum detail (0-1, maps to partial count)
    pub fn set_quantum_detail(&mut self, detail: f64) {
        self.quantum_detail = clamp01(detail);
    }

    /// Set Quantum base frequency (Hz)
    pub fn set_quantum_base(&mut self, base: f64) {
        self.quantum_base = clamp(base, 20.0, 2000.0);
    }

    /// Set Velvet impulse density (impulses/sec)
    pub fn set_velvet_density(&mut self, density: f64) {
        self.velvet_density = clamp(density, 1.0, 96000.0);
    }

    /// Reset all internal state
    pub fn reset(&mut self) {
        self.classic.reset();
        self.aureo.reset();
    }

    /// Generate one sample at the given sample rate
    pub fn next_sample(&mut self, rng: &mut Rng, sr: f64) -> f64 {
        match self.mode {
            NoiseMode::White | NoiseMode::Pink | NoiseMode::Brown => {
                self.classic.process(rng)
            }
            NoiseMode::Aureo => {
                self.configure_aureo(sr);
                self.aureo.process(rng, sr)
            }
            NoiseMode::Quantum => {
                self.configure_quantum();
                // Quantum always includes velvet residual
                let prev_velvet = self.aureo.velvet_on;
                self.aureo.velvet_on = true;
                let out = self.aureo.process(rng, sr);
                self.aureo.velvet_on = prev_velvet;
                out
            }
            NoiseMode::Velvet => {
                self.process_velvet(rng, sr)
            }
        }
    }
}

impl NoiseGen {
    /// Velvet noise with PhitRng — hardware entropy for unpredictable impulse timing
    pub fn process_velvet_phit(&self, rng: &mut crate::phit::PhitRng, sr: f64) -> f64 {
        let sr = if sr > 0.0 { sr } else { DEFAULT_SR };
        let prob = self.velvet_density / sr;
        let u = rng.next_f64();
        if u < prob {
            if rng.next_f64() < 0.5 { 1.0 } else { -1.0 }
        } else {
            0.0
        }
    }

    /// Push Aureo parameters into the aureo state before processing
    fn configure_aureo(&mut self, sr: f64) {
        self.aureo.planck_decay = clamp01(self.aureo_decay);
        self.aureo.prime_stride = self.aureo_stride.max(0.1);
        self.aureo.active_partials = self.aureo_harmonics.clamp(1, MAX_PARTIALS as i32);
        // Default base freq: sr * 0.25 * 1/sqrt(2) when not overridden
        self.aureo.base_freq = (sr * 0.25 * INV_SQRT2).max(10.0);
        self.aureo.velvet_on = true;
        self.aureo.velvet_amount = 0.10;
    }

    /// Push Quantum parameters into the aureo state before processing
    fn configure_quantum(&mut self) {
        let detail = clamp01(self.quantum_detail);
        let min_p = 6;
        let max_p = MAX_PARTIALS as i32;
        let target = (min_p as f64 + detail * (max_p - min_p) as f64).round() as i32;
        self.aureo.active_partials = target.clamp(4, max_p);
        self.aureo.planck_decay = map_phi_range(0.12, 0.52, 1.0 - detail);
        self.aureo.base_freq = clamp(self.quantum_base, 20.0, 2000.0);
        self.aureo.velvet_amount = 0.10;
    }

    /// Velvet noise: sparse random impulses
    fn process_velvet(&self, rng: &mut Rng, sr: f64) -> f64 {
        let sr = if sr > 0.0 { sr } else { DEFAULT_SR };
        let prob = self.velvet_density / sr;
        let u = rng.uni01();
        if u < prob {
            // Impulse fires: random sign
            if rng.uni01() < 0.5 { 1.0 } else { -1.0 }
        } else {
            0.0
        }
    }
}

impl Default for NoiseGen {
    fn default() -> Self {
        Self::new(NoiseMode::Pink)
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_noise_range() {
        let mut rng = Rng::new(42);
        let mut noise = NoiseColorState::new(NoiseColor::Pink, 0.8);

        for _ in 0..10000 {
            let v = noise.process(&mut rng);
            assert!(v.abs() < 5.0, "Noise value {} out of expected range", v);
        }
    }

    #[test]
    fn test_grain_kind_choice() {
        assert_eq!(GrainKind::choose(0.0, 0.5), GrainKind::Burst);
        assert_eq!(GrainKind::choose(1.0, 0.1), GrainKind::Stutter);
        assert_eq!(GrainKind::choose(1.0, 0.5), GrainKind::Aliaser);
        assert_eq!(GrainKind::choose(1.0, 0.8), GrainKind::Aliaser);
    }

    /// All 6 modes produce non-zero, non-NaN output
    #[test]
    fn test_all_modes_nonzero_valid() {
        let sr = 44100.0;
        let modes = [
            NoiseMode::White,
            NoiseMode::Pink,
            NoiseMode::Brown,
            NoiseMode::Aureo,
            NoiseMode::Quantum,
            NoiseMode::Velvet,
        ];
        for &mode in &modes {
            let mut rng = Rng::new(12345);
            let mut gen = NoiseGen::new(mode);
            gen.set_aureo_harmonics(12);
            gen.set_quantum_detail(0.7);
            gen.set_quantum_base(220.0);
            gen.set_velvet_density(2000.0);
            let mut any_nonzero = false;
            for _ in 0..44100 {
                let v = gen.next_sample(&mut rng, sr);
                assert!(!v.is_nan(), "Mode {:?} produced NaN", mode);
                assert!(!v.is_infinite(), "Mode {:?} produced Inf", mode);
                if v.abs() > 1.0e-12 {
                    any_nonzero = true;
                }
            }
            assert!(any_nonzero, "Mode {:?} produced only zeros over 1 second", mode);
        }
    }

    /// Aureo: output amplitude stays bounded with 32 partials
    #[test]
    fn test_aureo_bounded_32_partials() {
        let sr = 44100.0;
        let mut rng = Rng::new(42);
        let mut gen = NoiseGen::new(NoiseMode::Aureo);
        gen.set_aureo_harmonics(32);
        gen.set_aureo_decay(0.3);
        gen.set_aureo_stride(1.0);

        let mut max_abs = 0.0_f64;
        for _ in 0..44100 {
            let v = gen.next_sample(&mut rng, sr);
            max_abs = max_abs.max(v.abs());
        }
        // 32 partials with Planck decay should not blow up.
        // Each partial has amp <= 1.0 and decays exponentially,
        // so the sum should stay well under 32.
        assert!(
            max_abs < 20.0,
            "Aureo 32-partial peak {} exceeds safety bound", max_abs
        );
    }

    /// Velvet: impulse count should be roughly density over 1 second
    #[test]
    fn test_velvet_density() {
        let sr = 44100.0;
        let density = 2000.0;
        let mut rng = Rng::new(42);
        let mut gen = NoiseGen::new(NoiseMode::Velvet);
        gen.set_velvet_density(density);

        let n_samples = sr as usize;  // 1 second
        let mut impulse_count = 0usize;
        for _ in 0..n_samples {
            let v = gen.next_sample(&mut rng, sr);
            if v.abs() > 0.5 {
                impulse_count += 1;
            }
        }
        // Expect roughly 2000 impulses, allow +/- 20%
        let lo = (density * 0.80) as usize;
        let hi = (density * 1.20) as usize;
        assert!(
            impulse_count >= lo && impulse_count <= hi,
            "Velvet impulse count {} outside expected range [{}, {}]",
            impulse_count, lo, hi
        );
    }
}
