# Phi Lattice Vision — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Transform Aureonoise into a universal phi lattice — grains with personality, organic phi timing, triphase entropy, spatial phi correlation, and VHS removal.

**Architecture:** New `PhiLattice` struct as parameter oracle (`src/phi_lattice.rs`), wired into existing `Engine` via `spawn_grain()` and `schedule_gap_samples()`. `PhitRng` replaces `Rng` for grain decisions. Existing safety/validation preserved.

**Tech Stack:** Rust (PyO3), Python (DearPyGui), maturin build

**Build & test commands:**
```bash
unset CONDA_PREFIX && source .venv/bin/activate
maturin develop
pytest tests/ -v
python -m aureonoise.validate --strict
```

---

## Task 1: VHS Removal — Clean the Codebase

Remove all VHS-related code. User directive: "togliere VHS".

**Files:**
- Modify: `src/noise.rs:36-57` (GrainKind enum + choose())
- Modify: `src/lib.rs:103-105` (Params vhs_wow/vhs_flutter)
- Modify: `src/lib.rs:312-313` (Params defaults)
- Modify: `src/lib.rs:937-1006` (VHS LFO in process_block)
- Modify: `src/lib.rs:1063,1081,1116-1120` (spawn_grain vhs_mod param, grain processing VhsDrop)
- Modify: `src/lib.rs:1306` (spawn_grain signature)
- Modify: `python/aureonoise/presets.py` (remove vhs_wow/vhs_flutter from all presets)
- Modify: `python/aureonoise/app.py` (remove VHS sliders)
- Test: `tests/test_core.py`, full suite

**Step 1: Remove GrainKind::VhsDrop from noise.rs**

In `src/noise.rs`, remove the `VhsDrop` variant and update `choose()`:

```rust
// noise.rs — GrainKind enum (line 36-57)
pub enum GrainKind {
    Burst = 0,
    // VhsDrop = 1,  REMOVED
    Stutter = 2,
    Aliaser = 3,
}

impl GrainKind {
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
```

**Step 2: Remove VHS params from lib.rs Params struct**

Remove these fields from `Params` (lines 103-105):
```rust
// DELETE these:
pub vhs_wow: f64,
pub vhs_flutter: f64,
```

Remove their defaults (lines 312-313):
```rust
// DELETE these:
vhs_wow: 0.35,
vhs_flutter: 0.25,
```

**Step 3: Remove VHS LFO from process_block**

In `src/lib.rs` process_block:
- Remove lines 937-939 (wow_hz, flt_hz, inc_wow, inc_flt calculation)
- Remove lines 998-1006 (lfo_wow_phase/lfo_flut_phase update, wow/flt/vhs_mod calculation)
- Remove `vhs_mod` from spawn_grain call (line 1063)
- Remove `lfo_wow_phase` and `lfo_flut_phase` fields from Engine struct

**Step 4: Remove vhs_mod from spawn_grain signature**

Change `spawn_grain` signature (line 1306) from:
```rust
fn spawn_grain(&mut self, wi: usize, vhs_mod: f64, itd_scale: f64, coh_spatial_mod: f64, poly_pan_offset: f64)
```
to:
```rust
fn spawn_grain(&mut self, wi: usize, itd_scale: f64, coh_spatial_mod: f64, poly_pan_offset: f64)
```

**Step 5: Remove VhsDrop match arm from grain processing**

In process_block grain loop (lines 1116-1120), remove:
```rust
GrainKind::VhsDrop => {
    let att = 0.5 + 0.5 * (1.0 - vhs_mod.abs());
    s_l *= att;
    s_r *= att;
}
```

Also remove `vhs_mod * 0.25 * itd_scale` from the ITD read (line 1081):
```rust
// BEFORE:
let itd = grain.itd + vhs_mod * 0.25 * itd_scale;
// AFTER:
let itd = grain.itd;
```

**Step 6: Remove VHS from presets and GUI**

In `python/aureonoise/presets.py`: remove `vhs_wow` and `vhs_flutter` from `_FULL_DEFAULTS` and all preset overrides.

In `python/aureonoise/app.py`: remove VHS slider controls (search for "vhs").

**Step 7: Build and test**

```bash
maturin develop && pytest tests/ -v && python -m aureonoise.validate --strict
```

Expected: all existing tests pass (except any that explicitly test VHS — remove those too). Validation 253+ PASS.

**Step 8: Commit**

```bash
git add src/noise.rs src/lib.rs python/aureonoise/presets.py python/aureonoise/app.py
git commit -m "feat: remove VHS effects (user directive: togliere VHS)"
```

---

## Task 2: PhiLattice Core — New Module

Create `src/phi_lattice.rs` — the central phi ratio oracle.

**Files:**
- Create: `src/phi_lattice.rs`
- Modify: `src/lib.rs:1-28` (add `mod phi_lattice;` and `pub use`)
- Test: Rust unit tests in `phi_lattice.rs`

**Step 1: Write failing Rust tests**

Create `src/phi_lattice.rs` with tests first:

```rust
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
        // Allow very rare repeats (<5%) but not systematic
        assert!(repeats < 25, "Too many immediate repeats: {}/500", repeats);
    }

    #[test]
    fn test_grain_duration_is_phi_scaled() {
        let mut lattice = PhiLattice::new();
        let base = 200.0; // ms
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
    fn test_grain_adsr_sums_to_one() {
        let mut lattice = PhiLattice::new();
        for _ in 0..100 {
            let (a, d, _s, r) = lattice.grain_adsr();
            let sum = a + d + r;
            assert!(sum > 0.0 && sum <= 1.0, "ADSR a+d+r={} out of range", sum);
            assert!(a > 0.0 && d > 0.0 && r > 0.0, "ADSR segments must be positive");
        }
    }

    #[test]
    fn test_grain_adsr_phi_proportions() {
        let mut lattice = PhiLattice::new();
        for _ in 0..100 {
            let (a, d, _s, r) = lattice.grain_adsr();
            // Check that at least one pair of segments has a phi ratio relationship
            let ratios = [d / a, r / a, r / d];
            let has_phi = ratios.iter().any(|&ratio| {
                SACRED_RATIOS.iter().any(|&sr| (sr - ratio).abs() < 0.15)
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
```

**Step 2: Run tests to verify they fail**

```bash
maturin develop 2>&1 | tail -5
```

Expected: compilation error (PhiLattice not defined yet).

**Step 3: Implement PhiLattice**

```rust
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
/// Each template: (attack_ratio, decay_ratio, release_ratio) normalized
/// sustain_level is separate
const ADSR_TEMPLATES: [(f64, f64, f64, f64); 5] = [
    // (attack, decay, release, sustain_level)
    (0.146, 0.236, 0.618, 0.618),   // A=1/φ³, D=1/φ², R=1/φ, S=1/φ (fast attack)
    (0.236, 0.382, 0.382, 0.500),   // A=1/φ², D=1/φ, R=1/φ (balanced)
    (0.382, 0.236, 0.382, 0.618),   // A=1/φ, D=1/φ², R=1/φ (slow attack)
    (0.146, 0.146, 0.708, 0.382),   // A=1/φ³, D=1/φ³, R=long (percussive)
    (0.500, 0.118, 0.382, 0.750),   // A=1/2, D=tiny, R=1/φ (pad-like)
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
            last_ratio_idx: usize::MAX, // force no-repeat on first call
        }
    }

    /// Pick a sacred ratio organically.
    /// PhitRng entropy + weighted selection + OU drift = organic.
    /// Avoids repeating same ratio consecutively.
    pub fn next_ratio(&mut self) -> f64 {
        // Compute effective weights with OU drift
        let mut eff_weights = [0.0f64; 10];
        let mut total = 0.0;
        for i in 0..10 {
            let drift = 1.0 + 0.3 * self.ou_state * if i < 6 { 1.0 } else { -1.0 };
            let w = (self.weights[i] * drift).max(0.05);
            // Suppress last-used ratio
            let suppress = if i == self.last_ratio_idx { 0.1 } else { 1.0 };
            eff_weights[i] = w * suppress;
            total += eff_weights[i];
        }

        // Weighted random selection via PhitRng
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
    /// Clamps to avoid extremely short or long grains.
    pub fn grain_duration(&mut self, baselen_ms: f64) -> f64 {
        let ratio = self.next_ratio();
        clamp(baselen_ms * ratio, 10.0, 4000.0)
    }

    /// Inter-onset interval quantized to phi grid.
    /// base_interval_samples * sacred_ratio, clamped positive.
    pub fn onset_interval(&mut self, base_interval_samples: f64) -> f64 {
        let ratio = self.next_ratio();
        (base_interval_samples * ratio).max(1.0)
    }

    /// ADSR proportions — all segments in phi ratio to each other.
    /// Returns (attack_end, decay_end_offset, sustain_level, release_offset)
    /// where attack_end + decay_offset + sustain_gap + release_offset ≈ 1.0
    pub fn grain_adsr(&mut self) -> (f64, f64, f64, f64) {
        let idx = self.rng.next_u64() as usize % ADSR_TEMPLATES.len();
        let (a, d, r, s) = ADSR_TEMPLATES[idx];
        // Small PhitRng jitter (±10%) for organic variation
        let jitter = || 0.9 + 0.2 * self.rng.next_f64();
        let aj = (a * jitter()).max(0.01);
        let dj = (d * jitter()).max(0.01);
        let rj = (r * jitter()).max(0.01);
        // Renormalize so a+d+r <= 0.95 (leave 5% for sustain minimum)
        let total = aj + dj + rj;
        let scale = if total > 0.95 { 0.95 / total } else { 1.0 };
        (aj * scale, dj * scale, s, rj * scale)
    }

    /// Spectral tilt for a grain from phi-spaced discrete values.
    /// Centered around base_tilt, selects from TILT_SCALE.
    pub fn grain_spectral_tilt(&mut self, base_tilt: f64) -> f64 {
        // Find closest tilt value to base_tilt
        let base_idx = TILT_SCALE.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((**a - base_tilt).abs()).partial_cmp(&((**b - base_tilt).abs())).unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(3);

        // Offset by -1, 0, or +1 (stochastic)
        let offset = (self.rng.next_u64() % 3) as i32 - 1;
        let idx = (base_idx as i32 + offset).clamp(0, TILT_SCALE.len() as i32 - 1) as usize;
        TILT_SCALE[idx]
    }

    /// Spatial position on phi lattice [-1, 1].
    /// Positions quantized to sacred ratio divisions.
    pub fn spatial_position(&mut self) -> f64 {
        let ratio = self.next_ratio();
        let sign = if self.rng.next_f64() < 0.5 { -1.0 } else { 1.0 };
        // Map ratio to position magnitude: ratio/φ² clamped to [0, 1]
        let mag = (ratio / (PHI * PHI)).min(1.0);
        sign * mag
    }

    /// Step OU drift (call once per block). Slow evolution of weight emphasis.
    /// theta ~0.628 = 2π×0.1Hz (heart-brain coherence rhythm)
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
```

**Step 4: Register module in lib.rs**

Add after line 28 (`mod polyrhythm;`):
```rust
mod phi_lattice;
```

Add to pub use section (after line 52):
```rust
pub use phi_lattice::PhiLattice;
```

**Step 5: Build and run tests**

```bash
maturin develop && cargo test phi_lattice -p aureonoise -- --nocapture
```

Expected: 8 tests PASS.

**Step 6: Commit**

```bash
git add src/phi_lattice.rs src/lib.rs
git commit -m "feat: add PhiLattice core — sacred ratio oracle with triphase entropy"
```

---

## Task 3: Grain Personality — Extend Grain Struct

Add per-grain ADSR, spectral tilt, and tilt filter state to `Grain`.

**Files:**
- Modify: `src/grain.rs:11-79` (Grain struct fields)
- Modify: `src/grain.rs:108-143` (Default impl)
- Test: Rust compile + existing tests pass

**Step 1: Add personality fields to Grain**

In `src/grain.rs`, add after `dist_lp_z_r` (line 78):

```rust
    // Phi lattice grain personality
    /// Per-grain attack proportion [0, 1] from PhiLattice
    pub grain_attack: f64,
    /// Per-grain decay proportion [0, 1] from PhiLattice
    pub grain_decay: f64,
    /// Per-grain sustain level [0, 1] from PhiLattice
    pub grain_sustain: f64,
    /// Per-grain release proportion [0, 1] from PhiLattice
    pub grain_release: f64,
    /// Per-grain spectral tilt [-2, 0.5] from PhiLattice
    pub grain_tilt: f64,
    /// Tilt filter state (single-pole LP for per-grain coloring)
    pub tilt_z: f64,
```

**Step 2: Update Default impl**

In Default for Grain (line 108), add to the struct literal:

```rust
    grain_attack: 0.18,
    grain_decay: 0.28,
    grain_sustain: 0.55,
    grain_release: 0.30,
    grain_tilt: 0.0,
    tilt_z: 0.0,
```

**Step 3: Build and verify**

```bash
maturin develop && pytest tests/test_core.py -v
```

Expected: compiles, all core tests pass (new fields are just data, no behavior change yet).

**Step 4: Commit**

```bash
git add src/grain.rs
git commit -m "feat: add grain personality fields (ADSR, spectral tilt) to Grain struct"
```

---

## Task 4: Wire PhiLattice into Engine

Add `PhiLattice` to Engine struct, call `step_drift()` per block, populate grain personality in `spawn_grain()`.

**Files:**
- Modify: `src/lib.rs` (Engine struct, new(), process_block, spawn_grain)
- Test: `maturin develop && pytest tests/ -v`

**Step 1: Add new Params fields**

In `Params` struct (after `coherence_spatial` ~line 269):

```rust
    // Phi lattice (universal phi ratio oracle)
    #[pyo3(get, set)]
    pub phi_lattice_on: bool,
    /// Per-grain personality amount (0=all same, 1=max variation)
    #[pyo3(get, set)]
    pub phi_personality: f64,
    /// Timing snap to phi grid (0=free, 1=strict)
    #[pyo3(get, set)]
    pub phi_timing_strength: f64,
    /// Spatial snap to phi lattice (0=free, 1=strict)
    #[pyo3(get, set)]
    pub phi_spatial_strength: f64,
```

In `Params::default()` (after `coherence_spatial: false`):

```rust
    phi_lattice_on: true,
    phi_personality: 0.5,
    phi_timing_strength: 0.5,
    phi_spatial_strength: 0.5,
```

**Step 2: Add PhiLattice to Engine struct**

Find Engine struct definition. Add field:
```rust
    phi_lattice: phi_lattice::PhiLattice,
```

In `Engine::new()`, initialize:
```rust
    phi_lattice: phi_lattice::PhiLattice::new(),
```

**Step 3: Call step_drift in process_block**

In process_block, after the macro-modulation OU step (~line 921), add:

```rust
    // Phi lattice OU drift (slow evolution of ratio preferences)
    if self.params.phi_lattice_on {
        self.phi_lattice.step_drift(block_dt);
    }
```

**Step 4: Wire grain personality into spawn_grain**

In `spawn_grain()`, replace the envelope creation (lines 1554-1563):

```rust
    // Envelope — phi lattice personality or global params
    let (env_a, env_d, env_s, env_r) = if self.params.phi_lattice_on && self.params.phi_personality > 1e-6 {
        let (la, ld, ls, lr) = self.phi_lattice.grain_adsr();
        let p = clamp01(self.params.phi_personality);
        // Blend between global params and lattice personality
        (
            (1.0 - p) * self.params.env_attack + p * la,
            (1.0 - p) * self.params.env_decay + p * ld,
            (1.0 - p) * self.params.env_sustain + p * ls,
            (1.0 - p) * self.params.env_release + p * lr,
        )
    } else {
        (self.params.env_attack, self.params.env_decay, self.params.env_sustain, self.params.env_release)
    };

    // Store personality in grain
    grain.grain_attack = env_a;
    grain.grain_decay = env_d;
    grain.grain_sustain = env_s;
    grain.grain_release = env_r;

    grain.env = self.envelope.make_shape(
        env_a, env_d, env_s, env_r,
        gap_samples, len, pan.abs(),
    );
```

**Step 5: Wire per-grain spectral tilt**

In `spawn_grain()`, after envelope creation, add:

```rust
    // Per-grain spectral tilt from phi lattice
    if self.params.phi_lattice_on && self.params.phi_personality > 1e-6 {
        grain.grain_tilt = self.phi_lattice.grain_spectral_tilt(self.params.noise_slope);
    } else {
        grain.grain_tilt = self.params.noise_slope;
    }
```

In the grain processing loop (process_block, after ring buffer read ~line 1083), add per-grain tilt filter:

```rust
    // Per-grain spectral tilt (lightweight single-pole coloring)
    let tilt_diff = grain.grain_tilt - self.params.noise_slope;
    if tilt_diff.abs() > 0.1 {
        let alpha = 0.997_f64.powf(1.0 + tilt_diff.abs());
        grain.tilt_z = alpha * grain.tilt_z + (1.0 - alpha) * s_l;
        let tilt_s = if tilt_diff < 0.0 {
            grain.tilt_z  // darker: use LP output
        } else {
            s_l + 0.3 * (s_l - grain.tilt_z)  // brighter: boost HP residual
        };
        // Apply to both channels proportionally
        let ratio = if s_l.abs() > 1e-12 { tilt_s / s_l } else { 1.0 };
        s_l = tilt_s;
        s_r *= ratio;
    }
```

**Step 6: Bilateral onset clamp for grain personality**

In spawn_grain, after computing `env_a`, add safety clamp for bilateral presets:

```rust
    // Bilateral safety: clamp attack so onset < 5ms for CC stimulation
    let dur_ms = grain.dur as f64 * 1000.0 / self.sr;
    if self.params.bilateral_on && env_a * dur_ms > 5.0 {
        let env_a = 5.0 / dur_ms;
        grain.grain_attack = env_a;
        // Recompute envelope with clamped attack
        grain.env = self.envelope.make_shape(
            env_a, env_d, env_s, env_r,
            gap_samples, len, pan.abs(),
        );
    }
```

**Step 7: Build and test**

```bash
maturin develop && pytest tests/ -v && python -m aureonoise.validate --strict
```

Expected: all tests pass. Validation 253+ PASS (phi_personality defaults to 0.5, so behavior changes organically within existing tolerances).

**Step 8: Commit**

```bash
git add src/lib.rs src/grain.rs
git commit -m "feat: wire PhiLattice into Engine — per-grain ADSR personality + spectral tilt"
```

---

## Task 5: Organic Phi Timing

Replace exponential inter-arrival with phi-quantized timing. Grain duration also from phi lattice.

**Files:**
- Modify: `src/lib.rs:1366-1371` (spawn_grain duration)
- Modify: `src/lib.rs:1578-1630` (schedule_gap_samples)
- Test: existing suite + new Rust test

**Step 1: Replace grain duration with phi lattice**

In `spawn_grain()`, replace lines 1368-1371:

```rust
    // BEFORE:
    // let base = clamp(self.params.baselen_ms, ...) * 0.001 * self.sr;
    // let kexp = (2.0 * u1 - 1.0) * clamp01(self.params.len_phi);
    // let len = clamp(base * PHI.powf(kexp), ...);

    // AFTER: phi lattice organic duration
    let base_ms = clamp(self.params.baselen_ms, MIN_BASE_LENGTH_MS, 2000.0);
    let len = if self.params.phi_lattice_on && self.params.phi_timing_strength > 1e-6 {
        let phi_dur_ms = self.phi_lattice.grain_duration(base_ms);
        let free_kexp = (2.0 * u1 - 1.0) * clamp01(self.params.len_phi);
        let free_dur_ms = base_ms * PHI.powf(free_kexp);
        let t = clamp01(self.params.phi_timing_strength);
        // Blend: at 0 = free exponential, at 1 = strict phi lattice
        let dur_ms = (1.0 - t) * free_dur_ms + t * phi_dur_ms;
        clamp(dur_ms * 0.001 * self.sr, MIN_GRAIN_SAMPLES, self.sr * 4.0)
    } else {
        let base = base_ms * 0.001 * self.sr;
        let kexp = (2.0 * u1 - 1.0) * clamp01(self.params.len_phi);
        clamp(base * PHI.powf(kexp), MIN_GRAIN_SAMPLES, self.sr * 4.0)
    };
    grain.dur = len as u32;
```

**Step 2: Replace exponential gap with phi-quantized timing**

In `schedule_gap_samples()`, after the theta-gamma nesting section (~line 1604), replace the exponential inter-arrival:

```rust
    // AFTER theta-gamma nesting and rate calculation:
    if rate <= 1e-6 {
        return (self.sr * 0.25).max(1.0) as i32;
    }

    let base_interval = self.sr / rate;

    if self.params.phi_lattice_on && self.params.phi_timing_strength > 1e-6 {
        // Phi-quantized timing: interval from sacred ratio
        let phi_interval = self.phi_lattice.onset_interval(base_interval);

        // Hawkes burst modulation (unchanged)
        let hawkes_scale = if self.params.burst {
            1.0 / (1.0 + 0.3 * self.hawkes.lambda / rate.max(1e-3))
        } else {
            1.0
        };

        let t = clamp01(self.params.phi_timing_strength);
        // Blend with free exponential
        let u = self.rng.uni01().max(1.0e-12);
        let free_gap = (-u.ln() / rate.max(1e-3)) * self.sr;
        let gap = (1.0 - t) * free_gap + t * phi_interval * hawkes_scale;

        let base_rate = clamp(self.params.rate, 0.0, MAX_EVENT_RATE_HZ);
        let cap = if base_rate > 1e-6 { (30.0_f64.min(4.0 / base_rate)) * self.sr } else { 0.25 * self.sr };
        gap.min(cap).round().max(1.0) as i32
    } else {
        // Original exponential inter-arrival (unchanged)
        let t = self.sample_counter as f64 / self.sr;
        let mut lambda = rate * (1.0 + 0.2 * (TWO_PI * (t * INV_PHI)).sin());
        lambda = lambda.max(1e-3);
        if self.params.burst { lambda += 0.3 * self.hawkes.lambda; }
        let u = self.rng.uni01().max(1.0e-12);
        let mut gap_sec = -u.ln() / lambda;
        let base_rate = clamp(self.params.rate, 0.0, MAX_EVENT_RATE_HZ);
        let cap_sec = if base_rate > 1e-6 { 30.0_f64.min(4.0 / base_rate) } else { 0.25 };
        gap_sec = gap_sec.min(cap_sec);
        (gap_sec * self.sr).round().max(1.0) as i32
    }
```

**Step 3: Build and test**

```bash
maturin develop && pytest tests/ -v && python -m aureonoise.validate --strict
```

Expected: all tests pass. Overlap ratio may shift slightly — check validate output. If any preset fails temporal_continuity, adjust `phi_timing_strength` in that preset.

**Step 4: Commit**

```bash
git add src/lib.rs
git commit -m "feat: organic phi timing — grain duration and onset intervals from sacred ratios"
```

---

## Task 6: Triphase Entropy Wiring

Replace `Rng` with `PhitRng` for grain spawn decisions and velvet noise.

**Files:**
- Modify: `src/lib.rs:1316-1322` (spawn_grain quasi-random draws)
- Modify: `src/lib.rs:1042-1048` (noise generation for velvet)
- Modify: `src/noise.rs:611-623` (add process_velvet_phit method)
- Test: existing suite

**Step 1: Add process_velvet_phit to NoiseGen**

In `src/noise.rs`, add method to `impl NoiseGen` (non-pymethods block, after `process_velvet` ~line 622):

```rust
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
```

**Step 2: Wire PhitRng for velvet noise in process_block**

In process_block noise generation (~line 1042), change the velvet path:

```rust
    let mut nz = if self.params.noise_mode <= 2 {
        let white = self.rng.uni_pm1();
        self.spectral_tilt.process(white, self.params.noise_slope)
    } else if self.params.noise_mode == 5 {
        // Velvet: use PhitRng for triphase entropy
        let impulse = self.noise_gen.process_velvet_phit(&mut self.phit_rng, self.sr);
        self.spectral_tilt.process(impulse, self.params.noise_slope)
    } else {
        self.noise_gen.next_sample(&mut self.rng, self.sr)
    };
```

**Step 3: Use PhitRng for grain spawn quasi-random draws**

In `spawn_grain()`, replace lines 1320-1322:

```rust
    // BEFORE:
    // let u4 = self.rng.uni01();
    // let u5 = self.rng.uni01();
    // let u6 = self.rng.uni01();

    // AFTER: triphase entropy for grain decisions (brain can't predict)
    let u4 = if self.params.phi_lattice_on { self.phit_rng.next_f64() } else { self.rng.uni01() };
    let u5 = if self.params.phi_lattice_on { self.phit_rng.next_f64() } else { self.rng.uni01() };
    let u6 = if self.params.phi_lattice_on { self.phit_rng.next_f64() } else { self.rng.uni01() };
```

**Step 4: Build and test**

```bash
maturin develop && pytest tests/ -v && python -m aureonoise.validate --strict
```

Expected: all pass. PhitRng produces statistically equivalent output to Rng — behavior unchanged, entropy source upgraded.

**Step 5: Commit**

```bash
git add src/lib.rs src/noise.rs
git commit -m "feat: wire triphase entropy (PhitRng) into grain spawn + velvet noise"
```

---

## Task 7: Spatial Phi Lattice

Pan positions from phi lattice instead of uniform Weyl sequence.

**Files:**
- Modify: `src/lib.rs:1360-1365` (spawn_grain pan calculation)
- Test: existing suite

**Step 1: Replace pan initialization with phi lattice**

In `spawn_grain()`, replace lines 1360-1363:

```rust
    // BEFORE:
    // let mut pan = 2.0 * u3 - 1.0;
    // if self.params.lattice {
    //     pan += 0.25 * (2.0 * lat_u - 1.0) + 0.35 * oup;
    // }

    // AFTER: phi lattice spatial position
    let mut pan = if self.params.phi_lattice_on && self.params.phi_spatial_strength > 1e-6 {
        let phi_pan = self.phi_lattice.spatial_position();
        let free_pan = 2.0 * u3 - 1.0;
        let s = clamp01(self.params.phi_spatial_strength);
        (1.0 - s) * free_pan + s * phi_pan
    } else {
        2.0 * u3 - 1.0
    };
    if self.params.lattice {
        pan += 0.25 * (2.0 * lat_u - 1.0) + 0.35 * oup;
    }
```

**Step 2: Build and test**

```bash
maturin develop && pytest tests/ -v && python -m aureonoise.validate --strict
```

Expected: all pass. Spatial positions now cluster on phi lattice points.

**Step 3: Commit**

```bash
git add src/lib.rs
git commit -m "feat: spatial phi lattice — grain pan positions from sacred ratios"
```

---

## Task 8: Polyrhythm Asymmetric Grain Routing

Extend polyrhythm to support asymmetric grain distribution (3+2, 2+1, etc.).

**Files:**
- Modify: `src/polyrhythm.rs:95-109` (pan_offset method)
- Modify: `src/lib.rs:1441-1446` (polyrhythm pan in spawn_grain)
- Test: existing suite

**Step 1: Extend pan_offset for asymmetric distribution**

In `src/polyrhythm.rs`, replace `pan_offset` method (line 98):

```rust
    /// Get pan modulation with asymmetric distribution.
    /// p_pulse distributes grains in left hemisphere using phi spacing,
    /// q_pulse distributes in right hemisphere.
    /// `phi_pos` is an optional position from PhiLattice [0,1] for within-hemisphere spread.
    pub fn pan_offset_asymmetric(&self, amount: f64, phi_pos: f64) -> f64 {
        let amt = clamp(amount, 0.0, 1.0);
        if self.coincidence {
            0.0 // convergence — center
        } else if self.p_pulse {
            // Left hemisphere: spread within [-1, 0] using phi position
            -amt * (0.3 + 0.7 * clamp(phi_pos, 0.0, 1.0))
        } else if self.q_pulse {
            // Right hemisphere: spread within [0, 1] using phi position
            amt * (0.3 + 0.7 * clamp(phi_pos, 0.0, 1.0))
        } else {
            0.0
        }
    }
```

Keep the original `pan_offset` for backward compatibility.

**Step 2: Wire in spawn_grain**

In `spawn_grain()`, replace polyrhythm section (~lines 1441-1446):

```rust
    // Polyrhythm asymmetric routing
    if self.params.polyrhythm_on && poly_pan_offset.abs() > 1e-6 {
        if self.params.phi_lattice_on {
            let phi_pos = self.phi_lattice.spatial_position().abs(); // [0, 1]
            let asym_offset = self.polyrhythm.pan_offset_asymmetric(
                self.params.polyrhythm_amount, phi_pos
            );
            pan = clamp(pan + asym_offset, -1.0, 1.0);
        } else {
            pan = clamp(pan + poly_pan_offset, -1.0, 1.0);
        }
        grain.pan = pan;
    }
```

**Step 3: Build and test**

```bash
maturin develop && pytest tests/ -v
```

Expected: all pass.

**Step 4: Commit**

```bash
git add src/polyrhythm.rs src/lib.rs
git commit -m "feat: polyrhythm asymmetric grain routing — 3+2 phi-spaced distribution"
```

---

## Task 9: Preset Updates

Update presets with phi lattice params and velvet noise defaults for therapeutic presets.

**Files:**
- Modify: `python/aureonoise/presets.py` (_FULL_DEFAULTS + preset overrides)
- Test: `python -m aureonoise.validate --strict`

**Step 1: Add phi lattice params to _FULL_DEFAULTS**

In `presets.py`, add to `_FULL_DEFAULTS`:

```python
    "phi_lattice_on": True,
    "phi_personality": 0.5,
    "phi_timing_strength": 0.5,
    "phi_spatial_strength": 0.5,
```

**Step 2: Update therapeutic presets**

For each bilateral preset (EMDR, Hemispheric Bridge, CC Gentle, CC Maximum):
```python
    "phi_personality": 0.3,       # less variation for clinical precision
    "phi_timing_strength": 0.7,   # strong phi timing for therapeutic effect
    "phi_spatial_strength": 0.8,  # strong spatial phi correlation
```

For ambient presets (Sleep Pink, Theta Drift, etc.):
```python
    "phi_personality": 0.7,       # high variation for rich texture
    "phi_timing_strength": 0.5,   # moderate phi timing
    "phi_spatial_strength": 0.5,  # moderate spatial correlation
```

For entrainment presets (Gamma Focus, Delta Reset, etc.):
```python
    "phi_personality": 0.4,       # moderate variation
    "phi_timing_strength": 0.6,   # moderate-strong phi timing
    "phi_spatial_strength": 0.6,  # moderate-strong spatial correlation
    "noise_mode": 5,              # velvet noise (Alessio's preference)
    "velvet_density": 3000.0,     # dense enough for continuous texture
```

**Step 3: Run validation**

```bash
python -m aureonoise.validate --strict --verbose
```

Expected: 253+ PASS, 0 FAIL. Check that overlap ratios and onset times remain within bounds. Adjust phi_personality/phi_timing_strength if any preset violates therapeutic constraints.

**Step 4: Commit**

```bash
git add python/aureonoise/presets.py
git commit -m "feat: update presets with phi lattice params + velvet noise for entrainment"
```

---

## Task 10: Validation Updates

Add phi lattice validation checks and update existing profiles.

**Files:**
- Modify: `python/aureonoise/validate.py` (new checks, profile updates)
- Test: `python -m aureonoise.validate --strict`

**Step 1: Add phi_ratio_adherence check**

In `validate.py`, add new measurement function:

```python
def measure_phi_ratio_adherence(mono: np.ndarray, sr: float = 44100.0) -> Dict:
    """Measure how well grain intervals cluster on sacred phi ratios.
    Source: Klimesch 2012 (phi ratios in neural oscillations)."""
    SACRED = [0.382, 0.500, 0.618, 0.667, 1.0, 1.500, 1.618, 2.0, 2.618, 4.236]
    # Measure onset intervals via spectral flux
    onsets = measure_grain_onset_rate(mono, sr)
    if onsets.get("onset_count", 0) < 10:
        return {"phi_adherence": 0.0, "sample_count": 0}

    # Get intervals between consecutive onsets
    # (onset times extracted from the flux peaks)
    # Compute ratio of consecutive intervals
    # Count how many fall within 10% of a sacred ratio
    # Return adherence score [0, 1]
    return {"phi_adherence": score, "sample_count": n}
```

**Step 2: Add entropy_quality check**

```python
def measure_entropy_quality(left: np.ndarray, right: np.ndarray) -> Dict:
    """Verify PhitRng entropy in rendered audio via spectral feature analysis.
    Source: NIST SP 800-22 (randomness testing)."""
    # Compute short-time spectral variation
    # High entropy → low autocorrelation of spectral features
    # Return: autocorr_lag1 (should be < 0.3 for good entropy)
    return {"spectral_autocorr": autocorr, "entropy_sufficient": autocorr < 0.3}
```

**Step 3: Add checks to validation suite**

Add check 24 (phi_ratio_adherence) and check 25 (entropy_quality) to the main validation loop. Only run when phi_lattice_on=True in the preset.

**Step 4: Update TherapeuticProfile**

Add `phi_lattice: bool = True` field to TherapeuticProfile. Set to False for presets that opt out.

**Step 5: Run full validation**

```bash
python -m aureonoise.validate --strict --verbose
```

Expected: 260+ checks, 0 FAIL.

**Step 6: Commit**

```bash
git add python/aureonoise/validate.py
git commit -m "feat: add phi ratio adherence and entropy quality validation checks"
```

---

## Task 11: GUI Updates

Add phi lattice controls, remove VHS, update preset groups.

**Files:**
- Modify: `python/aureonoise/app.py` (new sliders, remove VHS)
- Test: launch app manually

**Step 1: Add phi lattice slider group**

In `app.py`, find where slider groups are defined. Add new group:

```python
with dpg.collapsing_header(label="Phi Lattice", default_open=True):
    dpg.add_slider_float(label="Personality", default_value=0.5,
        min_value=0.0, max_value=1.0,
        callback=lambda s, a: self.set_param("phi_personality", a))
    dpg.add_slider_float(label="Timing Strength", default_value=0.5,
        min_value=0.0, max_value=1.0,
        callback=lambda s, a: self.set_param("phi_timing_strength", a))
    dpg.add_slider_float(label="Spatial Strength", default_value=0.5,
        min_value=0.0, max_value=1.0,
        callback=lambda s, a: self.set_param("phi_spatial_strength", a))
    dpg.add_checkbox(label="Phi Lattice", default_value=True,
        callback=lambda s, a: self.set_param("phi_lattice_on", a))
```

**Step 2: Verify VHS sliders are removed** (should be done in Task 1)

Search `app.py` for any remaining "vhs" references and remove.

**Step 3: Test GUI launch**

```bash
python -m aureonoise
```

Verify: phi lattice controls visible, VHS controls gone, presets load correctly.

**Step 4: Commit**

```bash
git add python/aureonoise/app.py
git commit -m "feat: add phi lattice GUI controls, verify VHS removal"
```

---

## Task 12: Final Validation & Integration Test

Run full test suite, validate all presets, verify no regressions.

**Files:**
- Test: `tests/`, `python -m aureonoise.validate`

**Step 1: Full pytest**

```bash
maturin develop && pytest tests/ -v
```

Expected: 196+ passed, 0 failed.

**Step 2: Full validation**

```bash
python -m aureonoise.validate --strict --verbose 2>&1 | tail -20
```

Expected: 260+ PASS, 0 FAIL, warnings acceptable.

**Step 3: Rust tests**

```bash
cargo test -p aureonoise -- --nocapture 2>&1 | tail -20
```

Expected: all Rust tests pass including new phi_lattice tests.

**Step 4: Listen test**

```bash
python -m aureonoise
```

Load EMDR Bilateral, Hemispheric Bridge, Gamma Focus. Verify:
- Grains have organic variation (not mechanical repetition)
- Timing feels alive (not metronomic)
- Spatial positions cluster on phi points
- No VHS artifacts

**Step 5: Commit**

```bash
git add -A
git commit -m "feat: phi lattice vision — complete integration and validation"
```
