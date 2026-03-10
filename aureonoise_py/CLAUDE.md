# Aureonoise — Claude Guidelines

## Project Overview

Sintetizzatore di rumore granulare terapeutico. Rust DSP core (PyO3/maturin) + Python GUI (DearPyGui 2.x) + Python audio wrapper (sounddevice). Il motore produce noise colorato con sintesi granulare, spazializzazione bilaterale phi-ratio, e modulazioni terapeutiche (binaural beats, isochronic tones, tinnitus notch filtering, modal resonators).

**Stack:** Rust → `_core.cpython-312-darwin.so` via maturin → Python

**Build:** `maturin develop` (dev) / `maturin develop --release` (perf)

**Run:** `python -m aureonoise` oppure `python python/aureonoise/app.py`

**Validate:** `python -m aureonoise.validate --strict`

**Test:** `pytest tests/`

## Architecture — Key Facts

- **Rust entry**: `src/lib.rs` — `Engine` struct, `Params` (60+ parametri), `process_block()` loop per-sample, `spawn_grain()` scheduling, PyO3 exports
- **Signal chain**: noise → SR gain → tinnitus notch → ring buffer → grain spawn/process → modal resonator → binaural/isochronic → DVF near-field → externalization → room reverb → soft_tanh clip
- **Phi ratio**: φ=1.618 permea tutto — pan, timing, spatial taps, frequency ratios, head geometry
- **Python layer**: `audio.py` (sounddevice I/O), `app.py` (DearPyGui GUI), `presets.py` (20 preset), `validate.py` (validation scientifica), `analysis.py` (analisi spettrale/stereo)
- **Tests**: `pytest tests/` — 7 file di test (core, bilateral, burst, noise slope, phi model, preset validation, therapeutic)

## Boundaries — Chi Tocca Cosa

Ogni file sorgente appartiene a un solo agente. Le ownership sono dichiarate nei file agente in `.claude/agents/` e nelle rules in `.claude/rules/`. Nessun agente modifica file di un altro senza dichiararlo esplicitamente.

**Infrastruttura condivisa (read-only per tutti):**
- `src/constants.rs` — costanti globali (PHI, OUT_DRIVE, MAX_GRAINS)
- `src/math.rs` — utility matematiche (clamp, soft_tanh, woodworth_itd)
- `Cargo.toml`, `pyproject.toml` — configurazione build

## Moduli con agente specializzato

6 agenti con ownership esclusiva. Ogni modifica Rust **deve** passare `python -m aureonoise.validate --strict` prima di essere considerata completa.

### Engine Core (`engine-specialist`)
- **Scope**: `src/lib.rs`, `src/grain.rs`, `src/ring.rs`, `src/envelope.rs`, `src/rng.rs`, `src/weyl.rs`, `src/math.rs`, `src/constants.rs`
- **Tests**: `tests/test_core.py`
- **Dipendenze read-only**: tutti gli altri moduli Rust (li chiama, non li possiede)
- **Prodotto**: signal chain funzionante, zero-allocation nel callback, spawn_grain corretto
- **Responsabilita terapeutica**: real-time safety, no glitch, no allocation in process_block

### Spatial & Bilateral (`spatial-specialist`)
- **Scope**: `src/dialogue.rs`, `src/phi_model.rs`, `src/binaural.rs`, `src/phit.rs`
- **Tests**: `tests/test_bilateral_quality.py`, `tests/test_phi_model.py`
- **Dipendenze read-only**: `src/lib.rs` (capire process_block), `src/constants.rs`, `src/math.rs`
- **Prodotto**: stimolazione bilaterale corretta per corpus callosum, head model, binaural beats
- **Responsabilita terapeutica**: EMDR bilateral 1.0-1.3 Hz, onset <5ms, phi-ratio panning, ITD/ILD Woodworth

### Noise & Stochastic (`noise-specialist`)
- **Scope**: `src/noise.rs`, `src/stoch.rs`, `src/burst.rs`
- **Tests**: `tests/test_noise_slope.py`, `tests/test_burst.py`
- **Dipendenze read-only**: `src/lib.rs`, `src/grain.rs` (envelope interaction), `src/constants.rs`, `src/math.rs`, `src/rng.rs`, `src/weyl.rs`
- **Prodotto**: noise spettralmente corretto (pink -1 dB/oct, brown -2 dB/oct), Hawkes timing, OU modulation
- **Responsabilita terapeutica**: correttezza spettrale verificata, stochastic resonance calibrata

### Effects & Modulation (`effects-specialist`)
- **Scope**: `src/dvf.rs`, `src/isochronic.rs`, `src/modal.rs`, `src/polyrhythm.rs`, `src/room.rs`, `src/tinnitus.rs`, `src/external.rs`
- **Tests**: nessun test dedicato (copertura via test_therapeutic.py)
- **Dipendenze read-only**: `src/lib.rs`, `src/constants.rs`, `src/math.rs`
- **Prodotto**: DVF near-field, 40Hz gamma isochronic, modal contralateral mirror, tinnitus notch, room reverb
- **Responsabilita terapeutica**: 40Hz MIT GENUS protocol, tinnitus notch >20dB, body resonance safety

### Scientific Validation (`validation-specialist`)
- **Scope**: `python/aureonoise/validate.py`, `python/aureonoise/analysis.py`, `python/aureonoise/presets.py`, `tests/`
- **Dipendenze read-only**: tutto `src/` (verifica, non modifica), `python/aureonoise/audio.py`
- **Prodotto**: validation suite con soglie derivate dalla letteratura, presets terapeutici, test suite completa
- **Responsabilita terapeutica**: GUARDIAN — ogni modifica a qualsiasi modulo deve passare la validation suite. Le soglie sono derivate da paper peer-reviewed. Modificarle richiede citazione.

### Application & GUI (`app-specialist`)
- **Scope**: `python/aureonoise/app.py`, `python/aureonoise/audio.py`, `python/aureonoise/__init__.py`
- **Dipendenze read-only**: `python/aureonoise/presets.py`, `python/aureonoise/analysis.py`, `src/lib.rs` (PyO3 API)
- **Prodotto**: GUI DearPyGui, audio engine real-time, session management
- **Responsabilita terapeutica**: safety warnings, session timer, volume limiter, coherence visualization

## Code Conventions

- **Rust**: no allocation in `process_block()`, `#[inline]` su hot path, `f64` per DSP
- **Python**: type hints, dataclass per dati strutturati, numpy per signal processing
- **Imports**: stdlib → third-party → local
- **Nomenclatura**: snake_case Rust, snake_case Python, kebab-case per file agente
- **Commit**: inglese, presente indicativo, Co-Authored-By se da agente
- **Soglie terapeutiche**: ogni soglia in validate.py ha commento con fonte bibliografica

## Common Patterns

- **Nuovo modulo DSP Rust**: creare `src/nome.rs`, aggiungere `mod nome;` in `lib.rs`, esporre via PyO3 se serve Python-side
- **Nuovo preset**: aggiungere in `presets.py` con tutti i parametri (no bleed da defaults), aggiungere `TherapeuticProfile` in `validate.py`
- **Nuovo test**: aggiungere in `tests/test_*.py`, assicurarsi che `maturin develop` sia stato eseguito

## Research Documents

- `research/SYNTHESIS.md` — Sintesi terapeutica, 4 tier di evidenza, safety data
- `research/NEUROSCIENCE_INTEGRATION.md` — 8 aree neuroscientifiche, parameter mapping
- `research/RUSSIAN_UNCONVENTIONAL.md` — BAC, Slezin, Bekhtereva, ricerca russa non convenzionale
- `research/papers/` — Paper peer-reviewed archiviati con abstract e key findings

## Critical Issues — Priority

### P0 — Safety — RESOLVED
1. ~~No amplitude limiting~~ → Body resonance notch filters at 6.5 Hz (Q=2.5) and 19 Hz (Q=4.0)
2. ~~No epilepsy warning~~ → `isochronic_seizure_risk()` + -6dB attenuation 8-25 Hz + GUI warning
3. ~~No volume limiter~~ → Hard ceiling at 0.891 (-1 dBFS) in audio.py
4. ~~No session dose control~~ → Session timer + dose warnings in GUI

### P1 — Correctness — RESOLVED
5. ~~SpectralTilt brown~~ → Single-pole + white mix, raw ~-0.8 dB/oct
6. ~~Grain envelope linear~~ → envelope_shape param (0=linear, 1=Hann raised cosine)
7. ~~Coherence threshold~~ → Range [0.0, 2.0] in validate.py
8. ~~EMDR bilateral rate~~ → 1.2 Hz (Rousseau 2020 center of range)

### P2 — Completeness — RESOLVED
9. ~~phi_model disconnected~~ → Pinna amplitude + distance (gain+LP) + air absorption wired in spawn_grain
10. ~~0.1 Hz macro-modulation~~ → OU process (tau=1.59, sigma=0.25) modulates width ±15%
11. ~~BAC feedback loop~~ → Already existed: coherence → temperature + bilateral_rate. Enhanced with macro OU

### Remaining — Open
- phi_model.rs torso response (requires per-grain delay line ~64 samples, deferred)
- Pinna full 5-tap FIR (currently amplitude-only, FIR deferred for performance reasons)

## Dependencies & Environment

- **Rust**: edition 2021, `pyo3` 0.22
- **Python**: 3.12, `dearpygui` 2.x, `sounddevice`, `numpy`, `scipy` (optional)
- **Build**: `pip install maturin && maturin develop`
- **macOS**: PortAudio via `brew install portaudio` per sounddevice

## Learned Rules

<!-- Regole accumulate da correzioni. Formato: [DATA] DON'T: ... DO: ... Context: ... -->
- [2026-03-07] DON'T: assumere che bilateral = continuous pan. DO: usare eventi discreti con onset <5ms e gap di silenzio per stimolazione corpus callosum. Context: EMDR research richiede trasferimento callosale discreto.
- [2026-03-07] DON'T: assumere SR 44100 Hz. DO: auto-detect native SR del device audio. Context: PortAudio errore con SR mismatch.
- [2026-03-10] DON'T: aspettare coherence_mean() < 1.0. DO: aspettare coherence in [0.6, 1.8] — il range e by design in dialogue.rs. Context: coherence_target clamped a [0.6, 1.8].
- [2026-03-10] DON'T: modificare soglie terapeutiche senza citazione. DO: ogni soglia in validate.py deve avere commento con paper di riferimento. Context: "non si scherza con queste cose."
