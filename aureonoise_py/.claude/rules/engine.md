---
globs:
  - "src/lib.rs"
  - "src/grain.rs"
  - "src/ring.rs"
  - "src/envelope.rs"
  - "src/rng.rs"
  - "src/weyl.rs"
  - "src/math.rs"
  - "src/constants.rs"
  - "tests/test_core.py"
---

# Engine Core — Regole di contesto

## Zero allocation in process_block

Il callback audio e real-time. Qualsiasi allocation (Vec, Box, String, format!) causa glitch udibili.

Pattern corretto:
```rust
// Pre-allocare in new() o reset()
let mut buffer = [0.0f64; MAX_GRAINS];

// Usare in process_block()
for i in 0..active_count {
    buffer[i] = compute(i);
}
```

Pattern VIETATO:
```rust
// MAI in process_block
let v = Vec::new();  // allocation
let s = format!("debug {}", x);  // allocation
let b = Box::new(data);  // allocation
```

## Grain pool

`GrainPool` gestisce un array fisso di MAX_GRAINS (32) grani. Nessuna allocation dinamica.
- `spawn()`: trova il primo slot libero (grain.on == false), ritorna Option<&mut Grain>
- Se tutti i 32 slot sono occupati, lo spawn fallisce silenziosamente (non e un errore)
- `active_count()`: conta i grani con `on == true`

## Ring buffer masking

RING_SIZE = 131072 = 2^17. L'indice e mascherato con `& (RING_SIZE - 1)` per wrap-around senza divisione.

```rust
self.buffer[self.write_pos & (RING_SIZE - 1)] = sample;
```

Non cambiare RING_SIZE a un valore non-potenza-di-2.

## soft_tanh output stage

```rust
fn soft_tanh(x: f64) -> f64 {
    if x.abs() < 1.0 { x } else { x.signum() * (1.0 - 1.0 / (x.abs() * 2.0 + 1.0)) }
}
```

Applicato come `soft_tanh(OUT_DRIVE * sample)` con OUT_DRIVE = 1.2. Questo genera armoniche da saturazione — e intenzionale per warmth. Se modifichi OUT_DRIVE, riesegui `tests/test_noise_slope.py`.

## Envelope ADSR

L'envelope in envelope.rs usa segmenti lineari. Questo colora lo spettro dei grani corti (agisce come filtro bandpass). Se serve un'opzione raised-cosine/Hann per grani spettralmente neutri, il cambio va qui.

## Test

```bash
maturin develop && pytest tests/test_core.py -v
```
