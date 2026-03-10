---
globs:
  - "python/aureonoise/validate.py"
  - "python/aureonoise/analysis.py"
  - "python/aureonoise/presets.py"
  - "tests/test_therapeutic.py"
  - "tests/test_preset_validation.py"
  - "tests/test_bilateral_quality.py"
  - "tests/test_burst.py"
  - "tests/test_core.py"
  - "tests/test_noise_slope.py"
  - "tests/test_phi_model.py"
---

# Validation & Testing — Regole di contesto

## Soglie — citazione obbligatoria

Ogni soglia in validate.py ha (o deve avere) un commento con la fonte:
```python
EMDR_RATE = (1.0, 1.3)  # Rousseau 2020, Shapiro 2001
ITD_MAX_US = 800  # Woodworth model, 17.5cm head width
CC_MAX_ATTACK_MS = 5.0  # Corpus callosum IHTT literature
```

Se aggiungi una nuova soglia, DEVI citare il paper. Se non trovi un paper, la soglia e provvisoria e va marcata `# PROVISIONAL — needs citation`.

## _FULL_DEFAULTS — no parameter bleed

```python
_FULL_DEFAULTS = {
    "noise_mode": 0,
    "spectral_tilt": 0.0,
    "grain_rate": 20.0,
    # ... TUTTI i 60+ parametri con valori default
}

def get_preset(name):
    params = _FULL_DEFAULTS.copy()  # SEMPRE partire da defaults
    params.update(PRESET_OVERRIDES[name])  # sovrascrivere solo quelli specifici
    return params
```

MAI fare `prev_preset.update(new_values)` — causa parameter bleed.

## TherapeuticProfile — struttura

```python
@dataclass
class TherapeuticProfile:
    name: str
    category: str  # "bilateral", "entrainment", "general", "tinnitus"
    bilateral: bool
    isochronic: bool
    binaural: bool
    tinnitus_notch: bool
    expected_slope: tuple  # (min, max) dB/oct
    expected_stereo_corr: tuple  # (min, max) Pearson
    notes: str = ""
```

Ogni preset terapeutico ha un profilo. I profili determinano quali check vengono eseguiti.

## coherence range

```python
# CORRETTO
if coherence < 0.0 or coherence > 2.0:  # engine range [0.6, 1.8] with margin
    fail("coherence out of range")

# SBAGLIATO
if coherence > 1.0:  # engine produce fino a 1.8 by design
    fail("coherence > 1")
```

## Spectral slope — fattori confondenti

La slope misurata include:
1. Noise generator (noise.rs) — il contributo primario
2. Grain envelope windowing (envelope.rs) — filtraggio bandpass sui grani corti
3. soft_tanh output saturation (lib.rs) — genera armoniche, comprime picchi

Le soglie di tolleranza (es. pink ±1.0 dB/oct) tengono conto di #2 e #3. Se un altro specialista modifica envelope o output stage, le soglie potrebbero necessitare ricalibrazione.

## render_offline per test

```python
from aureonoise import Engine, Params
from aureonoise.audio import render_offline

params = Params()
# set params...
left, right = render_offline(params, duration_sec=5.0, sample_rate=44100.0)
# analyze...
```

Usare SEMPRE render_offline per test deterministici. Mai il real-time engine.

## Safety checks — stato attuale

Implementati:
- signal_present, peak_level, spectral_slope, spectral_fit_r2
- stereo_correlation, itd_us, ild_db, bilateral_symmetry, lr_alternation
- binaural_decorrelation, isochronic_am, cc_onset_time
- handshake_engagement, coherence, plv
- tinnitus_notch_depth

**NON implementati (P0 — urgente):**
- body_resonance_energy (5-8 Hz, 19 Hz)
- contraindication_screening (epilepsy risk 8-25 Hz)
- session_dose_validation (duration limits per frequency band)

## Test commands

```bash
# Validation suite completa
python -m aureonoise.validate --strict --verbose

# Singolo preset
python -m aureonoise.validate --preset "EMDR Bilateral" --verbose

# Test suite pytest
pytest tests/ -v

# Solo test specifici
pytest tests/test_noise_slope.py -v
pytest tests/test_bilateral_quality.py -v
```
