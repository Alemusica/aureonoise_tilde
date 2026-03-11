# Phi Lattice Vision — Design Document

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Transform Aureonoise into a universal phi lattice where everything is dialogue — grains have individual personality, timing falls on organic phi ratios, spatial geometry is phi-correlated, noise comes from triphase entropy, and VHS effects are removed.

**Architecture:** Hybrid approach (C). New `PhiLattice` struct as parameter oracle, integrated incrementally into existing modules. Existing safety/validation framework (253/256 PASS) preserved. PhitRng (triphase hardware entropy) replaces standard Rng for all stochastic grain decisions.

**Tech Stack:** Rust DSP core (PyO3), existing module structure. No new dependencies.

---

## 1. PhiLattice Core — `src/phi_lattice.rs`

### What

Central struct that generates phi-ratio parameters for the entire system. Not a replacement for existing modules — an oracle they consult.

### Design

```rust
use crate::phit::PhitRng;
use crate::constants::PHI;

/// The 10 sacred ratios (from dialogue.rs, validated by Klimesch 2012)
/// Phi family: therapeutic/desynchronization (weight 1.0)
/// Harmonic family: coupling/active cognition (weight 0.5)
const PHI_RATIOS: [f64; 10] = [
    0.382, 0.618, 1.0, 1.618, 2.618, 4.236,  // phi
    0.500, 0.667, 1.500, 2.000,                // harmonic
];

const PHI_WEIGHTS: [f64; 10] = [
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  // phi full weight
    0.5, 0.5, 0.5, 0.5,             // harmonic half weight
];

pub struct PhiLattice {
    rng: PhitRng,              // triphase entropy source
    weights: [f64; 10],        // organic weights (drift slowly via OU)
    ou_state: f64,             // OU drift for weight evolution
    ou_theta: f64,             // OU time constant (slow: ~0.1 Hz)
    last_ratio_idx: usize,     // avoid immediate repetition
}
```

### API

```rust
impl PhiLattice {
    pub fn new() -> Self;

    /// Pick a phi ratio organically.
    /// Always returns a sacred ratio, but WHICH one is stochastic.
    /// PhitRng entropy + weighted selection + OU drift = organic.
    /// Avoids repeating the same ratio consecutively.
    pub fn next_ratio(&mut self) -> f64;

    /// Grain duration: baselen_ms * sacred_ratio
    pub fn grain_duration(&mut self, baselen_ms: f64) -> f64;

    /// Inter-onset interval quantized to phi grid
    pub fn onset_interval(&mut self, base_interval_samples: f64) -> f64;

    /// ADSR proportions — all segments in phi ratio to each other.
    /// Returns (attack, decay, sustain_level, release) normalized to [0,1].
    /// Example: A=1/phi^2, D=1/phi, S=0.618, R=1.0 (relative)
    pub fn grain_adsr(&mut self) -> (f64, f64, f64, f64);

    /// Spectral tilt for a grain from phi-spaced discrete values.
    /// Returns a value in [-2.0, +0.5] chosen from a phi-spaced scale.
    pub fn grain_spectral_tilt(&mut self, base_tilt: f64) -> f64;

    /// Spatial position on phi lattice [-1, 1].
    pub fn spatial_position(&mut self) -> f64;

    /// Step OU drift (call once per block, not per sample)
    pub fn step_drift(&mut self, dt: f64);
}
```

### Ownership

- Created in `Engine::new()` as `self.phi_lattice: PhiLattice`
- Consulted in `spawn_grain()` and `schedule_gap_samples()`
- `step_drift()` called once per block in `process_block()` setup section

---

## 2. Grain Personality — `src/grain.rs` + `src/envelope.rs`

### What

Each grain born with individual ADSR, spectral color, and noise character. Currently grains only have `envelope_shape` (0=linear, 1=Hann) and shared attack/release from params.

### Changes to `Grain` struct

```rust
// NEW fields in Grain
pub grain_attack: f64,      // individual attack proportion [0,1]
pub grain_decay: f64,       // individual decay proportion [0,1]
pub grain_sustain: f64,     // individual sustain level [0,1]
pub grain_release: f64,     // individual release proportion [0,1]
pub grain_tilt: f64,        // per-grain spectral tilt [-2, +0.5]
```

### How it works

1. In `spawn_grain()`: call `phi_lattice.grain_adsr()` to get per-grain ADSR proportions
2. Pass these to `envelope.make_shape()` instead of using global `env_attack`/`env_release`
3. Call `phi_lattice.grain_spectral_tilt(params.noise_slope)` for per-grain color
4. Store `grain_tilt` in grain struct
5. In process loop: when reading from ring buffer, apply per-grain spectral tilt via a lightweight single-pole filter (not a full SpectralTilt per grain — too expensive)

### Per-grain spectral tilt (lightweight)

```rust
// In grain processing loop, after ring buffer read:
// Simple single-pole for per-grain coloring
let tilt_diff = grain.grain_tilt - base_tilt;
if tilt_diff.abs() > 0.1 {
    // Only apply if significantly different from base
    let alpha = 0.997_f64.powf(1.0 + tilt_diff.abs());
    grain.tilt_z = alpha * grain.tilt_z + (1.0 - alpha) * sample;
    sample = if tilt_diff < 0.0 { grain.tilt_z } else { sample + 0.3 * (sample - grain.tilt_z) };
}
```

Cost: 1 multiply + 1 add per grain per sample (negligible).

### Bilateral constraint

For CC/EMDR presets: `grain_attack` clamped so that `grain_attack * dur_ms < 5.0ms`. The lattice generates, the therapeutic constraint clamps.

---

## 3. Organic Phi Timing — `spawn_grain()` + `schedule_gap_samples()`

### What

Grain onsets always fall on phi ratios relative to the base rate, but which ratio is organic/stochastic. Currently: exponential inter-arrival with Hawkes modulation + theta-gamma nesting.

### Changes to `schedule_gap_samples()`

```rust
// CURRENT: gap = -ln(u) / lambda (exponential)
// NEW: gap = base_interval * phi_lattice.next_ratio()
// Still modulated by Hawkes burst intensity and theta-gamma nesting

fn schedule_gap_samples(&mut self) {
    // Check dialogue's queued Fibonacci gaps first (unchanged)
    if let Some(gap) = self.dialogue.pop_queued_gap() {
        self.samples_to_next = gap;
        return;
    }

    // Base interval from rate
    let base_interval = self.sr / effective_rate;

    // Phi-quantized interval
    let phi_interval = self.phi_lattice.onset_interval(base_interval);

    // Hawkes modulation (unchanged — burst intensity scales the interval)
    let hawkes_scale = 1.0 / (1.0 + self.hawkes.intensity() * burst_sensitivity);

    // Theta-gamma nesting quantization (unchanged)
    let gap = phi_interval * hawkes_scale;
    // ... theta-gamma snap ...

    self.samples_to_next = gap.round().max(1.0) as usize;
}
```

### Changes to `spawn_grain()` duration

```rust
// CURRENT: dur = baselen_ms * PHI^kexp, kexp = (2*u1-1) * len_phi
// NEW: dur = phi_lattice.grain_duration(baselen_ms)
// This produces baselen_ms * sacred_ratio, organic selection
let dur_ms = self.phi_lattice.grain_duration(baselen_ms);
let dur_samples = (dur_ms * self.sr / 1000.0).round().max(1.0) as usize;
```

---

## 4. Triphase Entropy Integration — `src/lib.rs`

### What

Replace standard `Rng` (xorshift64, deterministic seed) with `PhitRng` (hardware entropy from CPU clock phase) for grain personality decisions. The brain literally cannot predict the pattern because the entropy source is the CPU's own phase noise.

### Strategy

NOT a wholesale replacement. `Rng` stays for:
- OU process steps (needs stable statistics, not maximal entropy)
- Hawkes ticks (same)
- Noise generation base white noise (high throughput, doesn't need hardware entropy)

`PhitRng` used for:
- All grain spawn decisions via `PhiLattice` (which wraps PhitRng)
- Velvet noise impulse timing (sparse, so cost is low)
- Burst position modulation entropy

### Implementation

`PhiLattice` owns its `PhitRng`. Engine also keeps `self.phit_rng` for direct use in velvet noise.

In noise generation path (process_block per-sample noise section):
```rust
NoiseMode::Velvet => {
    // Use PhitRng for velvet impulse decisions
    self.noise_gen.process_velvet_phit(&mut self.phit_rng, self.sr)
}
```

New method in `NoiseGen`:
```rust
pub fn process_velvet_phit(&self, rng: &mut PhitRng, sr: f64) -> f64 {
    let prob = self.velvet_density / sr;
    let u = rng.next_f64();
    if u < prob {
        if rng.next_f64() < 0.5 { 1.0 } else { -1.0 }
    } else {
        0.0
    }
}
```

---

## 5. VHS Removal

### What

Remove `GrainKind::VhsDrop` and all VHS-related code. User directive: "togliere VHS".

### Files affected

- `src/noise.rs`: Remove `GrainKind::VhsDrop` variant. `GrainKind::choose()` redistributes probabilities to Burst/Stutter/Aliaser.
- `src/lib.rs`: Remove `vhs_lfo`-related code in process_block (wow/flutter VHS simulation), remove `vhs_mod` parameter from spawn_grain. Remove VHS aging/distortion in grain processing loop.
- `src/lib.rs` Params: Remove `glitch_vhs_amount` parameter if it exists.
- `python/aureonoise/presets.py`: Remove any VHS-related preset parameters.
- `python/aureonoise/app.py`: Remove VHS slider/controls from GUI.

### Migration

`GrainKind::choose()` becomes:
```rust
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
```

---

## 6. Spatial Phi Lattice

### What

All spatial positions geometrically correlated by phi. Currently: pan is quasi-random (Weyl) with dialogue correction and bilateral oscillator. Need: positions that fall on a phi lattice in virtual space.

### Design

`PhiLattice::spatial_position()` generates pan positions from the sacred ratio set:
```rust
pub fn spatial_position(&mut self) -> f64 {
    // Position from [-1, 1] quantized to phi lattice points
    let ratio = self.next_ratio();
    let sign = if self.rng.next_f64() < 0.5 { -1.0 } else { 1.0 };
    // Map ratio to position: 0.382 -> near center, 1.618 -> wide
    let pos = (ratio / PHI).min(1.0);
    sign * pos
}
```

### Integration in spawn_grain

```rust
// CURRENT: pan = 2*u3 - 1 (uniform random from Weyl)
// NEW: pan = phi_lattice.spatial_position()
// Then dialogue, burst, bilateral, polyrhythm still modify as before
let base_pan = self.phi_lattice.spatial_position();
```

The existing dialogue correction, hemisphere coupling, bilateral oscillator, and polyrhythm offset all still apply on top. The lattice sets the starting position, the pipeline refines it.

### Polyrhythm asymmetric distribution

Currently polyrhythm gives: p_pulse -> left, q_pulse -> right, coincidence -> center. Extend to support asymmetric grain counts:

```rust
// In spawn_grain, after polyrhythm pan_offset:
// 3:2 means 3 grains per p-cycle, 2 per q-cycle
// The phi_lattice distributes them organically within each hemisphere
if poly_clock.p_pulse {
    // 3 grains queued for left hemisphere
    // spread within left half using phi spacing
    pan = -phi_lattice.spatial_position().abs();
} else if poly_clock.q_pulse {
    // 2 grains queued for right hemisphere
    pan = phi_lattice.spatial_position().abs();
}
```

This gives the "3 da una parte, 2 dall'altra o mescolati" pattern.

---

## 7. Velvet Noise as Default Therapeutic Noise

### What

Velvet noise is Alessio's preferred noise type. Currently it's mode 5 (opt-in). For therapeutic presets, velvet should be the recommended base.

### Changes

- Update therapeutic preset defaults to use `noise_mode: 5` (Velvet) where appropriate
- Keep white/pink/brown available for specific therapeutic profiles that need colored noise
- Velvet impulse timing uses PhitRng (Section 4) for maximal entropy

### Velvet + spectral tilt

Velvet produces sparse impulses (+1/-1). To get colored velvet:
1. Velvet impulses fed through the existing SpectralTilt filter
2. This produces colored sparse noise — velvet density with spectral shape

```rust
// In process_block noise section:
NoiseMode::Velvet => {
    let impulse = self.noise_gen.process_velvet_phit(&mut self.phit_rng, self.sr);
    // Apply spectral tilt to velvet impulses
    self.spectral_tilt.process(impulse, params.noise_slope)
}
```

---

## 8. Parameters — New Params Fields

```rust
// New Params fields (PyO3-exposed)
pub phi_lattice_on: bool,       // enable phi lattice (default true)
pub phi_personality: f64,       // 0-1: how much per-grain personality (0=all same, 1=max variation)
pub phi_timing_strength: f64,   // 0-1: how strictly timing snaps to phi grid
pub phi_spatial_strength: f64,  // 0-1: how strictly positions fall on phi lattice
```

These allow continuous control: at 0, the system behaves like before (backward compatible). At 1, full phi lattice vision.

---

## 9. Validation Impact

### Existing checks preserved

All 256 existing checks remain valid. The phi lattice adds organic variation within bounds, doesn't violate therapeutic constraints.

### Therapeutic clamps

- CC onset < 5ms: `grain_attack` clamped regardless of lattice output
- EMDR rate 1.0-1.3 Hz: bilateral oscillator rate unchanged
- Spectral slope: per-grain tilt centered on preset's base slope, validation measures aggregate
- Overlap ratio: maintained by rate + duration, lattice varies duration around the phi center

### New validation checks to add

- Check 24: `phi_ratio_adherence` — verify that grain durations/gaps actually cluster on sacred ratios
- Check 25: `entropy_quality` — verify PhitRng monobit/chi-squared on rendered audio spectral features

---

## 10. Files Modified

| File | Change | Owner |
|------|--------|-------|
| `src/phi_lattice.rs` | **NEW** — PhiLattice struct | engine-specialist |
| `src/lib.rs` | Wire PhiLattice into spawn_grain, schedule_gap_samples, process_block. Remove VHS code. Wire PhitRng for velvet. | engine-specialist |
| `src/grain.rs` | Add personality fields to Grain struct | engine-specialist |
| `src/envelope.rs` | Accept per-grain ADSR from lattice | engine-specialist |
| `src/noise.rs` | Remove VhsDrop. Add `process_velvet_phit()`. | noise-specialist |
| `src/burst.rs` | Use PhitRng for burst entropy | noise-specialist |
| `src/polyrhythm.rs` | Asymmetric grain routing (3+2 pattern) | effects-specialist |
| `python/aureonoise/presets.py` | Add phi_lattice params, update noise modes, remove VHS | validation-specialist |
| `python/aureonoise/app.py` | Add phi_lattice controls, remove VHS GUI | app-specialist |
| `python/aureonoise/validate.py` | Add checks 24-25, update profiles | validation-specialist |

---

## Implementation Order

1. **VHS Removal** (clean first, no regressions)
2. **PhiLattice core** (src/phi_lattice.rs — new file, no existing code broken)
3. **Grain personality** (grain.rs + envelope.rs extension)
4. **Triphase wiring** (PhitRng into lattice + velvet)
5. **Organic phi timing** (schedule_gap_samples + spawn_grain duration)
6. **Spatial phi lattice** (spawn_grain pan)
7. **Polyrhythm asymmetric routing** (polyrhythm.rs extension)
8. **Velvet as default** (preset updates)
9. **Validation updates** (new checks 24-25)
10. **GUI updates** (new controls, VHS removal)

Each task builds on the previous. TDD: test first, implement, validate.
