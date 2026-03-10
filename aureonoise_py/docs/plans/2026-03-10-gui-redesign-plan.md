# GUI Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Refactor `app.py` from 14 tabs to 4 tabs with collapsing headers, extracting layout into builder methods.

**Architecture:** Pure layout refactor — same callbacks, same widget tags, same engine API. The `run()` method becomes a slim orchestrator calling `_build_presets_tab()`, `_build_sound_tab()`, `_build_space_brain_tab()`, `_build_monitor_tab()`. Collapsible sections via `dpg.collapsing_header()`.

**Tech Stack:** Python 3.12, DearPyGui 2.x, existing aureonoise Rust engine (PyO3)

**Design Doc:** `docs/plans/2026-03-10-gui-redesign-design.md`

---

### Task 1: Write app launch smoke test

**Files:**
- Create: `tests/test_app_smoke.py`

**Step 1: Write the test**

```python
"""Smoke test: app instantiation and DPG context creation."""
import pytest


def test_app_instantiation():
    """AureonoiseApp can be created without crashing."""
    from aureonoise.app import AureonoiseApp
    app = AureonoiseApp(sample_rate=44100.0, block_size=512)
    assert app.audio is not None
    assert app.running is False


def test_imports():
    """All app.py imports resolve."""
    from aureonoise.app import COLORS, NOISE_MODES, MODAL_PRESETS, main
    assert len(COLORS) > 0
    assert len(NOISE_MODES) == 6
    assert len(MODAL_PRESETS) == 4
```

**Step 2: Run test**

```bash
source .venv/bin/activate && pytest tests/test_app_smoke.py -v
```

Expected: PASS (2 tests)

**Step 3: Commit**

```bash
git add tests/test_app_smoke.py
git commit -m "test: add app launch smoke test before GUI refactor"
```

---

### Task 2: Rewrite app.py — complete restructure

This is the core task. It rewrites `python/aureonoise/app.py` as a layout-only refactor. **No callback logic changes.** All widget tags preserved exactly.

**Files:**
- Modify: `python/aureonoise/app.py` (full rewrite)

**Critical rules:**
1. Every widget tag (`sl_*`, `cb_*`, `txt_*`, `meter_*`, `radio_*`, `combo_*`, `grp_*`, `btn_*`) MUST be preserved exactly
2. All callback methods (`_on_*`, `_set_param`, `_update_*`, `_apply_preset`, `_sync_gui`, etc.) copy verbatim — no logic changes
3. Module-level helpers (`_safe_set`, `_safe_configure`, `_safe_enable`, `main`) copy verbatim
4. `_slider()`, `_slider_int()` copy verbatim
5. Color palette, imports, constants — unchanged

**Step 1: Read current app.py completely**

Read the entire file to understand all widget tags and their locations.

**Step 2: Write the new app.py**

New file structure (in this exact order):

```
# Docstring + imports (lines 1-25 — UNCHANGED)
# Constants: COLORS, NOISE_MODES, MODAL_PRESETS (lines 27-47 — UNCHANGED)

class AureonoiseApp:
    __init__()                    # UNCHANGED

    def run(self):                # SLIM — ~40 lines
        dpg.create_context()
        self._setup_theme()
        with dpg.window(tag="main", ...):
            # Header (same as before)
            with dpg.tab_bar():
                self._build_presets_tab()
                self._build_sound_tab()
                self._build_space_brain_tab()
                self._build_monitor_tab()
        # Viewport + meter thread + main loop + cleanup (same)

    def _setup_theme(self):       # EXTRACTED from run() lines 67-100
        # Global theme + preset button themes

    # ── Tab builders ────────────────────────────────
    def _build_presets_tab(self):  # NEW — from old "Therapeutic" tab + transport + meters
    def _build_sound_tab(self):   # NEW — merges Timing+Envelope+Timbre+Noise+Spectral+Modal+Stochastic+Polyrhythm
    def _build_space_brain_tab(self):  # NEW — merges Spatial+Bilateral+Binaural+Dialogue+Feedback+Tinnitus
    def _build_monitor_tab(self): # NEW — merges Analysis+System

    # ── Widget helpers ──────────────────────────────
    _slider()                     # UNCHANGED
    _slider_int()                 # UNCHANGED

    # ── Parameter dispatch ──────────────────────────
    _SAFETY_PARAMS                # UNCHANGED
    _set_param()                  # UNCHANGED

    # ── Toggle callbacks ────────────────────────────
    _on_binaural_toggle()         # UNCHANGED
    _on_isochronic_toggle()       # UNCHANGED
    _on_dialogue_toggle()         # UNCHANGED
    _on_bilateral_toggle()        # UNCHANGED
    _on_feedback_toggle()         # UNCHANGED
    _on_coherence_spatial_toggle() # UNCHANGED
    _on_bilateral_nesting_toggle() # UNCHANGED
    _on_polyrhythm_toggle()       # UNCHANGED
    _on_modal_toggle()            # UNCHANGED
    _on_burst_toggle()            # UNCHANGED
    _update_feature_summary()     # UNCHANGED
    _update_safety_warnings()     # UNCHANGED
    _on_noise_mode_change()       # UNCHANGED

    # ── Preset logic ────────────────────────────────
    _make_preset_callback()       # UNCHANGED
    _apply_preset()               # UNCHANGED
    _sync_gui()                   # UNCHANGED

    # ── Audio device ────────────────────────────────
    _get_output_devices()         # UNCHANGED
    _on_device_change()           # UNCHANGED

    # ── Transport ───────────────────────────────────
    _on_play()                    # UNCHANGED
    _on_stop()                    # UNCHANGED
    _on_reset()                   # UNCHANGED

    # ── Session timer callbacks ─────────────────────
    _on_session_warning()         # UNCHANGED
    _on_session_limit()           # UNCHANGED

    # ── Analysis callback ───────────────────────────
    _on_analysis_report()         # UNCHANGED

    # ── Meter thread ────────────────────────────────
    _start_meter_thread()         # UNCHANGED

# Module-level helpers — UNCHANGED
_safe_set()
_safe_configure()
_safe_enable()
main()
```

**Step 2a: `_setup_theme()` method**

Extract from current `run()` lines 66-100:

```python
def _setup_theme(self):
    """Set up global theme and per-preset button themes."""
    with dpg.theme() as global_theme:
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg, COLORS["bg"])
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg, COLORS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_Text, COLORS["text"])
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrab, COLORS["accent"])
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrabActive, COLORS["accent"])
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (40, 40, 48))
            dpg.add_theme_color(dpg.mvThemeCol_Button, COLORS["accent_dim"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, COLORS["accent"])
            dpg.add_theme_color(dpg.mvThemeCol_Tab, (40, 40, 48))
            dpg.add_theme_color(dpg.mvThemeCol_TabActive, COLORS["accent_dim"])
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered, COLORS["accent"])
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, 4)
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 12, 12)
    dpg.bind_theme(global_theme)

    # Per-button themes for therapeutic presets
    self._btn_themes = {}
    for name, (_key, rgba) in FACTORY.items():
        with dpg.theme() as t:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, rgba)
                r, g, b, a = rgba
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered,
                                    (min(r + 30, 255), min(g + 30, 255),
                                     min(b + 30, 255), a))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,
                                    (max(r - 20, 0), max(g - 20, 0),
                                     max(b - 20, 0), a))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (255, 255, 255, 255))
        self._btn_themes[name] = t
```

**Step 2b: `_build_presets_tab()` method**

This is the LANDING TAB. Contains transport + safety + meters + preset grid.

Source widgets:
- Transport: from old `run()` lines 110-118
- Safety panel: from old `run()` lines 122-127
- Meters: from old `run()` lines 131-137
- Preset grid: from old "Therapeutic" tab lines 430-458

```python
def _build_presets_tab(self):
    """Build the Presets tab (landing page)."""
    with dpg.tab(label="Presets"):
        dpg.add_spacer(height=5)

        # Transport
        with dpg.group(horizontal=True):
            dpg.add_button(label="  Play", tag="btn_play", callback=self._on_play)
            dpg.add_button(label="  Stop", tag="btn_stop", callback=self._on_stop)
            dpg.add_button(label="  Reset", callback=self._on_reset)
            dpg.add_spacer(width=20)
            dpg.add_text("", tag="status", color=COLORS["text_dim"])
            dpg.add_spacer(width=20)
            dpg.add_text("00:00", tag="txt_session_timer", color=COLORS["text_dim"])

        dpg.add_spacer(height=4)

        # Safety warning panel
        with dpg.child_window(height=50, border=False, tag="safety_panel"):
            dpg.add_text("", tag="txt_warn_session", color=(220, 200, 60))
            dpg.add_text("", tag="txt_warn_epilepsy", color=(220, 60, 60))
            dpg.add_text("", tag="txt_warn_resonance", color=(220, 200, 60))

        dpg.add_spacer(height=4)

        # Meters
        with dpg.child_window(height=60, border=False):
            with dpg.group(horizontal=True):
                dpg.add_text("L", color=COLORS["meter_l"])
                dpg.add_progress_bar(tag="meter_l", default_value=0, width=200)
                dpg.add_spacer(width=10)
                dpg.add_text("R", color=COLORS["meter_r"])
                dpg.add_progress_bar(tag="meter_r", default_value=0, width=200)

        dpg.add_separator()
        dpg.add_spacer(height=10)

        # Preset grid
        dpg.add_text("Evidence-Based Presets", color=COLORS["section"])
        dpg.add_text(
            "Each preset configures noise, spatial, dialogue and bilateral "
            "parameters for a specific therapeutic context.",
            color=COLORS["text_dim"], wrap=520)
        dpg.add_spacer(height=10)

        row_items = list(FACTORY.items())
        for row_start in range(0, len(row_items), 3):
            row_slice = row_items[row_start:row_start + 3]
            with dpg.group(horizontal=True):
                for name, (_key, _rgba) in row_slice:
                    btn = dpg.add_button(
                        label=f"  {name}  ",
                        callback=self._make_preset_callback(name),
                        height=40,
                        width=170,
                    )
                    dpg.bind_item_theme(btn, self._btn_themes[name])
                    dpg.add_spacer(width=5)
            dpg.add_spacer(height=5)

        dpg.add_separator()
        dpg.add_spacer(height=5)
        dpg.add_text("", tag="txt_preset_active", color=COLORS["accent"])
        dpg.add_text("", tag="txt_feature_summary", color=COLORS["text_dim"])
```

**Step 2c: `_build_sound_tab()` method**

Merges: Noise (205-238) + Timing (146-150) + Envelope (185-191) + Timbre (193-202) + Spectral slope from Tinnitus (350-354) + Modal (276-303) + Stochastic (404-427) + Polyrhythm from Feedback (389-396)

```python
def _build_sound_tab(self):
    """Build the Sound tab — noise, grain, envelope, timbre, texture."""
    with dpg.tab(label="Sound"):
        dpg.add_spacer(height=5)

        # ── Noise ──
        with dpg.collapsing_header(label="Noise", default_open=True):
            dpg.add_text("Noise Mode", color=COLORS["section"])
            dpg.add_radio_button(
                NOISE_MODES,
                default_value="Pink",
                tag="radio_noise_mode",
                callback=self._on_noise_mode_change,
                horizontal=True,
            )
            dpg.add_spacer(height=5)

            # Aureo params (visible when mode=3)
            with dpg.group(tag="grp_aureo", show=False):
                dpg.add_text("Aureo Parameters", color=COLORS["section"])
                self._slider("aureo_decay", "Aureo Decay", 0.0, 1.0, 0.3)
                self._slider("aureo_stride", "Aureo Stride", 0.1, 4.0, 1.0)
                self._slider_int("aureo_harmonics", "Aureo Harmonics", 1, 32, 12)
                dpg.add_spacer(height=5)

            # Quantum params (visible when mode=4)
            with dpg.group(tag="grp_quantum", show=False):
                dpg.add_text("Quantum Parameters", color=COLORS["section"])
                self._slider("quantum_detail", "Quantum Detail", 0.0, 1.0, 0.7)
                self._slider("quantum_base", "Quantum Base (Hz)", 20.0, 2000.0, 220.0, log=True)
                dpg.add_spacer(height=5)

            # Velvet params (visible when mode=5)
            with dpg.group(tag="grp_velvet", show=False):
                dpg.add_text("Velvet Parameters", color=COLORS["section"])
                self._slider("velvet_density", "Velvet Density", 1.0, 96000.0, 2000.0, log=True)
                dpg.add_spacer(height=5)

        # ── Grain ──
        with dpg.collapsing_header(label="Grain", default_open=False):
            self._slider("rate", "Rate (Hz)", 0.1, 60.0, 8.0, log=True)
            self._slider("baselen_ms", "Base Length (ms)", 10.0, 2000.0, 120.0, log=True)
            self._slider("len_phi", "Length phi Spread", 0.0, 1.0, 0.8)

        # ── Envelope ──
        with dpg.collapsing_header(label="Envelope", default_open=False):
            self._slider("env_attack", "Attack", 0.01, 1.0, 0.18)
            self._slider("env_decay", "Decay", 0.01, 1.0, 0.28)
            self._slider("env_sustain", "Sustain", 0.0, 1.0, 0.55)
            self._slider("env_release", "Release", 0.01, 1.0, 0.30)

        # ── Timbre ──
        with dpg.collapsing_header(label="Timbre", default_open=False):
            self._slider("vhs_wow", "VHS Wow", 0.0, 1.0, 0.35)
            self._slider("vhs_flutter", "VHS Flutter", 0.0, 1.0, 0.25)
            self._slider("glitch_mix", "Glitch Mix", 0.0, 1.0, 0.5)
            dpg.add_separator()
            self._slider("srcrush_amt", "Sample Rate Crush", 0.0, 1.0, 0.2)
            self._slider("bitcrush_amt", "Bit Crush", 0.0, 1.0, 0.15)

        # ── Spectral ──
        with dpg.collapsing_header(label="Spectral", default_open=False):
            dpg.add_text(
                "Continuous tilt: 0=white, -1=pink, -2=brown.",
                color=COLORS["text_dim"], wrap=520)
            self._slider("noise_slope", "Slope (dB/oct)", -2.5, 0.5, -1.0)

        # ── Texture ──
        with dpg.collapsing_header(label="Texture", default_open=False):
            # Modal
            dpg.add_text("Modal Resonator", color=COLORS["section"])
            dpg.add_checkbox(
                label="Modal On", default_value=False, tag="cb_modal_on",
                callback=self._on_modal_toggle)
            dpg.add_spacer(height=5)
            dpg.add_combo(
                MODAL_PRESETS,
                default_value="Wood",
                label="Preset",
                tag="combo_modal_preset",
                callback=lambda s, a, u: self._set_param("modal_preset",
                    MODAL_PRESETS.index(a)),
                width=150,
            )
            dpg.add_spacer(height=5)
            self._slider("modal_mix", "Mix", 0.0, 1.0, 0.3)
            self._slider("modal_decay", "Decay", 0.0, 1.0, 0.5)
            self._slider("modal_mirror", "Mirror", 0.0, 1.0, 0.3)
            self._slider("modal_feedback", "Feedback", 0.0, 1.0, 0.1)
            dpg.add_spacer(height=8)
            dpg.add_text("Contralateral Mirror", color=COLORS["section"])
            dpg.add_text(
                "Routes modal response to opposite hemisphere of burst events.",
                color=COLORS["text_dim"],
            )
            self._slider("modal_contralateral", "Contralateral", 0.0, 1.0, 0.0)

            dpg.add_separator()
            dpg.add_spacer(height=5)

            # Stochastic
            dpg.add_text("Stochastic Processes", color=COLORS["section"])
            with dpg.group(horizontal=True):
                dpg.add_checkbox(label="Thermo", default_value=True,
                    tag="cb_thermo",
                    callback=lambda s, a, u: self._set_param(u, a),
                    user_data="thermo")
                dpg.add_checkbox(label="Lattice", default_value=True,
                    tag="cb_lattice",
                    callback=lambda s, a, u: self._set_param(u, a),
                    user_data="lattice")
                dpg.add_checkbox(label="Burst", default_value=True,
                    tag="cb_burst",
                    callback=self._on_burst_toggle,
                    user_data="burst")
            dpg.add_separator()
            self._slider("burst_floor", "Burst Floor", 0.0, 1.0, 0.3)
            self._slider("burst_phi_mix", "Burst Phi Mix", 0.0, 1.0, 0.5)
            dpg.add_separator()
            self._slider("temperature", "Temperature", 0.0, 1.0, 0.45)
            self._slider("lat_rate", "Lattice Rate", 1.0, 2000.0, 250.0, log=True)
            self._slider("lat_eps", "Lattice e", 0.01, 0.5, INV_PHI_CU)
            self._slider("lat_gamma", "Lattice g", 0.5, 3.0, PHI)
            self._slider("lat_sigma", "Lattice s", 0.01, 0.3, 0.06)

            dpg.add_separator()
            dpg.add_spacer(height=5)

            # Polyrhythm
            dpg.add_text("Polyrhythm Clock", color=COLORS["section"])
            dpg.add_checkbox(
                label="Polyrhythm On", default_value=False, tag="cb_polyrhythm_on",
                callback=self._on_polyrhythm_toggle)
            self._slider_int("polyrhythm_p", "P (left)", 2, 8, 3)
            self._slider_int("polyrhythm_q", "Q (right)", 2, 8, 2)
            self._slider("polyrhythm_rate", "Base Rate (Hz)", 0.1, 3.0, 0.5)
            self._slider("polyrhythm_amount", "Amount", 0.0, 1.0, 0.5)
```

**Step 2d: `_build_space_brain_tab()` method**

Merges: Spatial (153-183) + Binaural/Isochronic (306-334) + Dialogue (241-273) + Feedback (357-401 parts) + Tinnitus notch (337-346)

```python
def _build_space_brain_tab(self):
    """Build the Space & Brain tab — spatial, bilateral, entrainment, dialogue."""
    with dpg.tab(label="Space & Brain"):
        dpg.add_spacer(height=5)

        # ── Stereo Field ──
        with dpg.collapsing_header(label="Stereo Field", default_open=True):
            self._slider("width", "Stereo Width", 0.0, 2.0, 1.0)
            self._slider("itd_us", "ITD (us)", 0.0, 800.0, 600.0)
            self._slider("ild_db", "ILD (dB)", 0.0, 12.0, 6.0)
            dpg.add_separator()
            self._slider("hemis_coupling", "Hemisphere Coupling", 0.0, 1.0, 0.6)
            self._slider("spat_ipd", "IPD Amount", 0.0, 1.0, 0.6)
            self._slider("spat_shadow", "Head Shadow", 0.0, 1.0, 0.7)

        # ── Bilateral ──
        with dpg.collapsing_header(label="Bilateral", default_open=False):
            dpg.add_checkbox(
                label="Phi-Pan", default_value=False,
                callback=lambda s, a, u: self._set_param(u, a),
                user_data="phi_pan")
            dpg.add_checkbox(
                label="Bilateral On", default_value=False, tag="cb_bilateral_on",
                callback=self._on_bilateral_toggle)
            self._slider("bilateral_rate", "Bilateral Rate (Hz)", 0.3, 6.0, 1.0)
            self._slider("bilateral_amount", "Bilateral Amount", 0.0, 1.0, 0.8)

        # ── Externalization ──
        with dpg.collapsing_header(label="Externalization", default_open=False):
            self._slider("externalization", "Externalization", 0.0, 1.0, 0.0)
            self._slider("phi_distance", "Phi Distance (m)", 0.0, 10.0, 1.5)
            self._slider("phi_elev", "Phi Elevation (deg)", -90.0, 90.0, 0.0)
            dpg.add_separator()
            dpg.add_spacer(height=5)
            dpg.add_text("Room Reverb", color=COLORS["section"])
            self._slider("room_mix", "Room Mix", 0.0, 0.5, 0.0)

        # ── Entrainment ──
        with dpg.collapsing_header(label="Entrainment", default_open=False):
            dpg.add_text("Binaural Beat Generator", color=COLORS["section"])
            dpg.add_text(
                "Separate sine tones per ear. Beat frequency = "
                "difference between L and R carrier.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_spacer(height=5)
            dpg.add_checkbox(
                label="Binaural On", default_value=False, tag="cb_binaural_on",
                callback=self._on_binaural_toggle)
            self._slider("binaural_carrier_hz", "Carrier (Hz)", 100.0, 500.0, 250.0)
            self._slider("binaural_beat_hz", "Beat (Hz)", 0.5, 40.0, 6.0)
            self._slider("binaural_level", "Level", 0.0, 0.3, 0.08)

            dpg.add_separator()
            dpg.add_spacer(height=5)
            dpg.add_text("Isochronic Tone", color=COLORS["section"])
            dpg.add_text(
                "Pulsed carrier (Tukey-windowed AM). Mono, both ears.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_spacer(height=5)
            dpg.add_checkbox(
                label="Isochronic On", default_value=False, tag="cb_isochronic_on",
                callback=self._on_isochronic_toggle)
            self._slider("isochronic_carrier_hz", "Carrier (Hz)", 100.0, 500.0, 165.0)
            self._slider("isochronic_rate_hz", "Rate (Hz)", 1.0, 40.0, 10.0)
            self._slider("isochronic_duty", "Duty Cycle", 0.2, 0.8, 0.5)
            self._slider("isochronic_level", "Level", 0.0, 0.3, 0.10)

        # ── Dialogue ──
        with dpg.collapsing_header(label="Dialogue", default_open=False):
            dpg.add_text("Interhemispheric Coherence", color=COLORS["section"])
            dpg.add_checkbox(
                label="Dialogue On", default_value=True, tag="cb_dialogue_on",
                callback=self._on_dialogue_toggle)
            dpg.add_spacer(height=5)
            self._slider("dialogue_strength", "Strength", 0.0, 1.0, 0.6)
            self._slider("dialogue_memory", "Memory", 0.0, 1.0, 0.5)
            self._slider("dialogue_phi_mix", "Phi Mix", 0.0, 1.0, 0.75)

            dpg.add_separator()
            dpg.add_spacer(height=5)
            dpg.add_text("Engine Metrics", color=COLORS["section"])

            with dpg.group(horizontal=True):
                dpg.add_text("Coherence:", color=COLORS["text_dim"])
                dpg.add_progress_bar(
                    tag="meter_coherence", default_value=0, width=180,
                    overlay="0.000")

            with dpg.group(horizontal=True):
                dpg.add_text("Handshakes:", color=COLORS["text_dim"])
                dpg.add_text("0", tag="txt_handshake_count")

            with dpg.group(horizontal=True):
                dpg.add_text("Handshake Ratio:", color=COLORS["text_dim"])
                dpg.add_text("0.000", tag="txt_handshake_ratio")

            with dpg.group(horizontal=True):
                dpg.add_text("Mean Coherence:", color=COLORS["text_dim"])
                dpg.add_text("0.000", tag="txt_coherence_mean")

        # ── Feedback ──
        with dpg.collapsing_header(label="Feedback", default_open=False):
            dpg.add_text("Coherence Feedback Loop", color=COLORS["section"])
            dpg.add_text(
                "BAC-inspired closed-loop: high coherence calms, "
                "low coherence explores.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_spacer(height=5)
            dpg.add_checkbox(
                label="Feedback On", default_value=False, tag="cb_feedback_on",
                callback=self._on_feedback_toggle)
            self._slider("temp_ramp_sec", "Temp Ramp (s)", 0.0, 60.0, 0.0)

            dpg.add_separator()
            dpg.add_spacer(height=5)
            dpg.add_text("Theta-Gamma Nesting", color=COLORS["section"])
            dpg.add_checkbox(
                label="Bilateral Nesting", default_value=False, tag="cb_bilateral_nesting",
                callback=self._on_bilateral_nesting_toggle)

            dpg.add_separator()
            dpg.add_spacer(height=5)
            dpg.add_text("Coherence Spatial Morphing", color=COLORS["section"])
            dpg.add_text(
                "Width/ITD modulated by dialogue coherence.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_checkbox(
                label="Coherence Spatial", default_value=False, tag="cb_coherence_spatial",
                callback=self._on_coherence_spatial_toggle)

        # ── Tinnitus ──
        with dpg.collapsing_header(label="Tinnitus", default_open=False):
            dpg.add_text("Tinnitus Notch Filter", color=COLORS["section"])
            dpg.add_text(
                "4th-order Butterworth notch at your tinnitus frequency. "
                "Set to 0 to disable.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_spacer(height=5)
            self._slider("tinnitus_notch_hz", "Center Freq (Hz)", 0.0, 12000.0, 0.0)
            self._slider("tinnitus_notch_q", "Q Factor", 1.0, 20.0, 6.0)
```

**Step 2e: `_build_monitor_tab()` method**

Merges: Analysis (461-512) + System (515-560)

```python
def _build_monitor_tab(self):
    """Build the Monitor tab — analysis and system info."""
    with dpg.tab(label="Monitor"):
        dpg.add_spacer(height=5)

        # ── Analysis ──
        with dpg.collapsing_header(label="Analysis", default_open=True):
            dpg.add_text("Real-Time Signal Analysis", color=COLORS["section"])
            dpg.add_text(
                "Spectral, stereo, and therapeutic quality metrics. "
                "Updates every ~1 second while playing.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_spacer(height=8)

            with dpg.group(horizontal=True):
                dpg.add_text("Verdict:", color=COLORS["text_dim"])
                dpg.add_text("--", tag="txt_analysis_verdict",
                             color=COLORS["accent"])
            dpg.add_spacer(height=5)

            dpg.add_text("Spectral", color=COLORS["section"])
            dpg.add_text("Slope: -- dB/oct", tag="txt_an_slope",
                         color=COLORS["text"])
            dpg.add_text("Centroid: -- Hz", tag="txt_an_centroid",
                         color=COLORS["text"])

            dpg.add_separator()
            dpg.add_spacer(height=3)

            dpg.add_text("Stereo", color=COLORS["section"])
            dpg.add_text("Correlation: --", tag="txt_an_stereo",
                         color=COLORS["text"])
            dpg.add_text("ITD: -- us", tag="txt_an_itd",
                         color=COLORS["text"])
            dpg.add_text("ILD: 1k=-- dB  4k=-- dB", tag="txt_an_ild",
                         color=COLORS["text"])
            dpg.add_text("Bilateral sym: --", tag="txt_an_sym",
                         color=COLORS["text"])

            dpg.add_separator()
            dpg.add_spacer(height=3)

            dpg.add_text("Level", color=COLORS["section"])
            dpg.add_text("RMS: -- dBFS  Peak: -- dBFS",
                         tag="txt_an_level", color=COLORS["text"])

            dpg.add_separator()
            dpg.add_spacer(height=3)

            dpg.add_text("Warnings", color=COLORS["section"])
            dpg.add_text("--", tag="txt_an_warnings",
                         color=COLORS["meter_peak"], wrap=520)

        # ── System ──
        with dpg.collapsing_header(label="System", default_open=False):
            dpg.add_text("Audio Output", color=COLORS["section"])
            self._output_devices = self._get_output_devices()
            device_labels = [f"{idx}: {name}" for idx, name in self._output_devices]
            default_dev = sd.default.device[1] if sd else 0
            default_label = ""
            for idx, name in self._output_devices:
                if idx == default_dev:
                    default_label = f"{idx}: {name}"
                    break
            if not default_label and device_labels:
                default_label = device_labels[0]
            dpg.add_combo(
                device_labels,
                default_value=default_label,
                label="Output Device",
                tag="combo_audio_device",
                callback=self._on_device_change,
                width=350,
            )
            dpg.add_text(
                "Change device while stopped for best results.",
                color=COLORS["text_dim"], wrap=520)
            dpg.add_spacer(height=10)
            dpg.add_separator()

            dpg.add_spacer(height=5)
            dpg.add_input_int(
                label="Seed",
                default_value=20251010,
                callback=lambda s, a, u: self._set_param(u, a),
                user_data="seed",
                width=150,
            )
            dpg.add_spacer(height=10)
            dpg.add_text(f"Sample Rate: {self.audio.sample_rate} Hz", color=COLORS["text_dim"])
            dpg.add_text(f"Block Size: {self.audio.block_size}", color=COLORS["text_dim"])
            dpg.add_spacer(height=10)
            dpg.add_text("phi Constants:", color=COLORS["accent"])
            dpg.add_text(f"  phi = {PHI:.6f}", color=COLORS["text_dim"])
            dpg.add_text(f"  1/phi = {INV_PHI:.6f}", color=COLORS["text_dim"])
            dpg.add_text(f"  1/phi^2 = {INV_PHI_SQ:.6f}", color=COLORS["text_dim"])
            dpg.add_text(f"  1/phi^3 = {INV_PHI_CU:.6f}", color=COLORS["text_dim"])
```

**Step 2f: New slim `run()` method**

```python
def run(self):
    """Run the application."""
    dpg.create_context()
    self._setup_theme()

    with dpg.window(label="aureonoise", tag="main", no_title_bar=True):
        # Header
        dpg.add_text("aureonoise", color=COLORS["accent"])
        dpg.add_text("phi-based granular texture generator", color=COLORS["text_dim"])
        dpg.add_spacer(height=10)

        # Tab bar with 4 tabs
        with dpg.tab_bar():
            self._build_presets_tab()
            self._build_sound_tab()
            self._build_space_brain_tab()
            self._build_monitor_tab()

    # Viewport setup
    dpg.create_viewport(
        title="aureonoise",
        width=680,
        height=820,
        resizable=True,
    )
    dpg.setup_dearpygui()
    dpg.set_primary_window("main", True)
    dpg.show_viewport()

    # Start meter update thread
    self._start_meter_thread()

    # Main loop
    while dpg.is_dearpygui_running():
        dpg.render_dearpygui_frame()

    # Cleanup
    self.running = False
    self.audio.stop()
    dpg.destroy_context()
```

**Step 2g: Assemble the complete file**

Write the full `app.py` with:
1. Imports + constants (unchanged)
2. `AureonoiseApp.__init__()` (unchanged)
3. New slim `run()` (Step 2f)
4. `_setup_theme()` (Step 2a)
5. `_build_presets_tab()` (Step 2b)
6. `_build_sound_tab()` (Step 2c)
7. `_build_space_brain_tab()` (Step 2d)
8. `_build_monitor_tab()` (Step 2e)
9. All callbacks and helpers (lines 587-1125 — COPY VERBATIM, zero changes)
10. Module-level helpers + `main()` (lines 1127-1158 — COPY VERBATIM)

**Step 3: Run smoke test**

```bash
source .venv/bin/activate && pytest tests/test_app_smoke.py -v
```

Expected: PASS

---

### Task 3: Verify existing tests pass

**Step 1: Run full pytest suite**

```bash
source .venv/bin/activate && pytest tests/ -v
```

Expected: 196 passed (194 + 2 new smoke), 2 xfailed

**Step 2: Run validation suite**

```bash
source .venv/bin/activate && python -m aureonoise.validate --verbose
```

Expected: 211/212 passed, 0 failures, 1 warning (Delta Reset 3Hz absence seizure — intentional)

---

### Task 4: Commit the refactor

**Step 1: Stage and commit**

```bash
git add python/aureonoise/app.py tests/test_app_smoke.py
git commit -m "refactor: GUI 14 tabs → 4 tabs with collapsing headers

Restructure app.py layout:
- Presets (landing), Sound, Space & Brain, Monitor
- Collapsing headers for expandable parameter sections
- Extract _setup_theme(), 4 builder methods from monolithic run()
- No callback/logic changes, all widget tags preserved

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Widget Tag Checklist

Every tag from the original must exist in the new file:

**Sliders (sl_*):** rate, baselen_ms, len_phi, width, itd_us, ild_db, hemis_coupling, spat_ipd, spat_shadow, externalization, phi_distance, phi_elev, bilateral_rate, bilateral_amount, env_attack, env_decay, env_sustain, env_release, vhs_wow, vhs_flutter, glitch_mix, srcrush_amt, bitcrush_amt, aureo_decay, aureo_stride, aureo_harmonics, quantum_detail, quantum_base, velvet_density, dialogue_strength, dialogue_memory, dialogue_phi_mix, modal_mix, modal_decay, modal_mirror, modal_feedback, modal_contralateral, binaural_carrier_hz, binaural_beat_hz, binaural_level, isochronic_carrier_hz, isochronic_rate_hz, isochronic_duty, isochronic_level, tinnitus_notch_hz, tinnitus_notch_q, noise_slope, burst_floor, burst_phi_mix, temperature, lat_rate, lat_eps, lat_gamma, lat_sigma, polyrhythm_p, polyrhythm_q, polyrhythm_rate, polyrhythm_amount, temp_ramp_sec, room_mix

**Checkboxes (cb_*):** dialogue_on, bilateral_on, modal_on, binaural_on, isochronic_on, feedback_on, bilateral_nesting, coherence_spatial, polyrhythm_on, thermo, lattice, burst

**Text (txt_*):** session_timer, warn_session, warn_epilepsy, warn_resonance, preset_active, feature_summary, analysis_verdict, an_slope, an_centroid, an_stereo, an_itd, an_ild, an_sym, an_level, an_warnings, handshake_count, handshake_ratio, coherence_mean

**Meters:** meter_l, meter_r, meter_coherence

**Groups (grp_*):** aureo, quantum, velvet

**Other:** status, btn_play, btn_stop, radio_noise_mode, combo_modal_preset, combo_audio_device, safety_panel, main
