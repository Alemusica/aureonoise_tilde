# GUI Redesign — Design Document

**Date:** 2026-03-10
**Status:** Approved
**Scope:** `python/aureonoise/app.py`

## Problem

Current GUI has 14 tabs in a monolithic 1159-line file. The "Feedback" tab is a junk drawer with 6 unrelated features. `run()` is 500+ lines of inline layout. No separation between layout builders and callback logic.

## Decision

**4 tabs + collapsing headers** (Approach C with B expandability).

Primary user intent: stare bene. Preset-first, parameters accessible for exploration via expandable sections.

## Architecture

### Tab Structure

```
Tab 1: Presets (landing — always visible, no collapsing)
  - Transport: Play / Stop / Reset + status + session timer
  - Safety warnings panel
  - L/R meters
  - Preset grid (15 colored buttons, 5x3)
  - Active preset name + feature summary

Tab 2: Sound
  [▸] Noise — mode radio, conditional params (aureo/quantum/velvet)
  [▸] Grain — rate, base length, phi spread
  [▸] Envelope — ADSR
  [▸] Timbre — VHS wow/flutter, glitch mix, SR crush, bit crush
  [▸] Spectral — slope (dB/oct)
  [▸] Texture — modal resonator (on/off, preset, mix, decay, mirror, feedback,
                 contralateral) + stochastic (thermo/lattice/burst + params) +
                 polyrhythm (on/off, P, Q, rate, amount)

Tab 3: Space & Brain
  [▸] Stereo Field — width, ITD, ILD, head shadow, hemisphere coupling, IPD
  [▸] Bilateral — on/off, rate, amount
  [▸] Externalization — distance, elevation, room mix
  [▸] Entrainment — binaural (on/off, carrier, beat, level) +
                     isochronic (on/off, carrier, rate, duty, level)
                     mutual exclusion preserved
  [▸] Dialogue — on/off, strength, memory, phi mix +
                  metrics (coherence meter, handshake count/ratio, mean coherence)
  [▸] Feedback — feedback on, temp ramp, bilateral nesting, coherence spatial
  [▸] Tinnitus — notch center Hz, Q factor

Tab 4: Monitor
  [▸] Analysis — verdict, spectral (slope, centroid, spread),
                  stereo (correlation, ITD, ILD, bilateral sym),
                  level (RMS, peak), warnings
  [▸] System — audio device selector, seed, sample rate, block size, phi constants
```

### Default Expansion

- **Presets tab:** everything visible (no collapsing headers, it's the landing page)
- **Sound/Space & Brain/Monitor:** first section expanded, rest collapsed

### Code Architecture

Current: single `AureonoiseApp` class, `run()` builds everything inline.

Refactored:
```python
class AureonoiseApp:
    def run(self):
        dpg.create_context()
        self._setup_theme()
        with dpg.window(tag="main", ...):
            with dpg.tab_bar():
                self._build_presets_tab()
                self._build_sound_tab()
                self._build_space_brain_tab()
                self._build_monitor_tab()
        self._setup_viewport()
        self._start_meter_thread()
        self._main_loop()

    # ── Tab builders ──
    def _build_presets_tab(self): ...
    def _build_sound_tab(self): ...
    def _build_space_brain_tab(self): ...
    def _build_monitor_tab(self): ...

    # ── Theme ──
    def _setup_theme(self): ...
    def _setup_preset_themes(self): ...

    # ── Callbacks (grouped by section) ──
    # Transport
    def _on_play(self): ...
    def _on_stop(self): ...
    def _on_reset(self): ...
    # Parameter dispatch
    def _set_param(self, name, value): ...
    # Toggle handlers (mutual exclusion, dependencies)
    def _on_binaural_toggle(self, ...): ...
    def _on_isochronic_toggle(self, ...): ...
    def _on_dialogue_toggle(self, ...): ...
    def _on_bilateral_toggle(self, ...): ...
    # etc.

    # ── Widget helpers ──
    def _slider(self, ...): ...
    def _slider_int(self, ...): ...

    # ── Meters / updates ──
    def _start_meter_thread(self): ...
    def _update_feature_summary(self): ...
    def _update_safety_warnings(self): ...
    def _on_analysis_report(self, report): ...
```

### DearPyGui Collapsing Header API

```python
with dpg.collapsing_header(label="Noise", default_open=True):
    # noise widgets here
with dpg.collapsing_header(label="Grain", default_open=False):
    # grain widgets here
```

### What Does NOT Change

- `audio.py` — untouched, clean architecture
- `presets.py` — untouched
- `analysis.py` — untouched
- `validate.py` — untouched
- `__init__.py` — untouched
- All Rust code — untouched
- Color palette — same phi-based colors
- All parameter names — same tags, same dispatch
- Mutual exclusion logic — same behavior
- Safety warnings — same behavior
- Meter thread — same implementation

### What Changes

1. **14 tabs → 4 tabs** with collapsing headers inside
2. **`run()` 500+ lines → ~30 lines** orchestrator
3. **Layout code extracted** to 4 builder methods
4. **"Feedback" junk drawer eliminated** — contents distributed logically
5. **Tinnitus tab: spectral slope moved to Sound > Spectral**
6. **Preset tab becomes landing tab** (first tab, selected by default)

### Migration Checklist

Every widget tag (`sl_*`, `cb_*`, `txt_*`, `meter_*`, `radio_*`, `combo_*`, `grp_*`) must be preserved exactly. The refactor is layout-only — no parameter renaming, no callback logic change, no engine API change.

## Validation

After refactor:
1. `maturin develop && pytest tests/ -v` — all 194 pass
2. `python -m aureonoise.validate --strict` — 211/212 (1 intentional WARN)
3. Manual: launch app, load each preset, verify all sliders sync, verify meters update
4. Manual: verify mutual exclusion (binaural/isochronic toggle)
5. Manual: verify safety warnings fire on epilepsy-range isochronic
