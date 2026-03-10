---
globs:
  - "src/dialogue.rs"
  - "src/phi_model.rs"
  - "src/binaural.rs"
  - "src/phit.rs"
  - "tests/test_bilateral_quality.py"
  - "tests/test_phi_model.py"
---

# Spatial & Bilateral — Regole di contesto

## coherence_mean() range

`coherence_mean()` ritorna valori in [0.6, 1.8]. NON in [0, 1].

Il `coherence_target` e clamped a [0.6, 1.8] alla riga ~350 di dialogue.rs:
```rust
self.coherence_target = self.coherence_target.clamp(0.6, 1.8);
```

Questo e by design — il range rappresenta la distanza dal rapporto phi ideale. Se qualcuno segnala "coherence > 1.0 e un bug", la risposta e: no, e by design.

## EMDR bilateral rate

Range validato dalla ricerca: 1.0-1.3 Hz (Rousseau 2020, Shapiro 2001).
Il preset EMDR Bilateral usa 1.2 Hz (centro del range validato).
Il BilateralOscillator ha 4 fasi: L→center, center→R, R→center, center→L.
Full cycle = 4 fasi = 1/rate secondi.

## Corpus callosum onset

Onset < 5ms richiesto per stimolazione discreta del corpus callosum.
Calcolo: `onset_ms = env_attack * grain_dur_samples * 1000 / sample_rate`

Preset bilaterali corretti:
- EMDR Bilateral: env_attack = 0.03 → ~3.6ms onset
- Hemispheric Bridge: env_attack = 0.04 → ~4ms onset
- CC Gentle: env_attack = 0.03 → ~3.6ms onset
- CC Maximum: env_attack = 0.04 → ~4ms onset

## 10 SACRED_RATIOS

```rust
const SACRED_RATIOS: [f64; 10] = [
    PHI, INV_PHI, PHI*PHI, INV_PHI*INV_PHI,
    2.0, 0.5, PHI*PHI*PHI, INV_PHI*INV_PHI*INV_PHI,
    3.0, 1.0/3.0
];
```

Il ratio_score confronta il rapporto L/R amplitude contro ognuna di queste. Non ridurre il set.

## phi_model.rs — Stato delle funzioni

| Funzione | Stato | Param | Costo |
|----------|-------|-------|-------|
| `compute_head_result()` | CONNESSA | (sempre) | Basso (1 trig per grano) |
| `design_pinna_response()` | CONNESSA (amp-only) | `spat_pinna` 0-1 | Basso (amp_scale a spawn) |
| `compute_torso_response()` | DISCONNESSA | — | Medio (serve delay ~64 samples per grano) |
| `compute_distance_response()` | CONNESSA | `spat_distance` 0-1 | Basso (gain + LP per grano) |
| `compute_air_absorption()` | CONNESSA | `spat_distance` 0-1 | Basso (gain a spawn, ref 4kHz) |

Pinna full 5-tap FIR deferred: max delay ~12 samples, serve buffer [f64; 16] per grano per canale.
Torso deferred: delay ~40-96 samples, troppo grande per buffer per-grano.

## 0.1 Hz macro-modulation

OU `ou_macro` (tau=1.59, sigma=0.25) modula width ±15% per ritmo cuore-cervello.
Attivo solo con `feedback_on=true`. Stepped once per block (non per sample).
Applicato come moltiplicatore su `coh_spatial_mod` in process_block.

## Binaural beats — segnale separato

I binaural beats vanno aggiunti DOPO il processing dei grani nel signal chain, non miscelati nel ring buffer. Il tono sinusoidale puro e necessario per l'entrainment (Wahbeh 2007).

## Test

```bash
maturin develop && pytest tests/test_bilateral_quality.py tests/test_phi_model.py -v
```
