"""
aureonoise - GUI Application

DearPyGui-based interface for real-time control.
"""

from typing import Optional
import math
import threading
import time

try:
    import dearpygui.dearpygui as dpg
except ImportError:
    dpg = None

try:
    import sounddevice as sd
except ImportError:
    sd = None

from aureonoise import Params, PHI, INV_PHI, INV_PHI_SQ, INV_PHI_CU
from aureonoise.audio import AudioEngine
from aureonoise.presets import PresetBank, FACTORY


# phi-based color palette
COLORS = {
    "bg": (18, 18, 22),
    "panel": (28, 28, 34),
    "accent": (180, 140, 90),      # Golden
    "accent_dim": (120, 95, 60),
    "text": (220, 215, 205),
    "text_dim": (140, 135, 125),
    "meter_l": (90, 180, 140),
    "meter_r": (180, 140, 90),
    "meter_peak": (220, 80, 80),
    "coherence": (130, 180, 220),
    "section": (160, 140, 100),
}

# Noise mode labels (index -> display name)
NOISE_MODES = ["White", "Pink", "Brown", "Aureo", "Quantum", "Velvet"]

# Modal preset labels (index -> display name)
MODAL_PRESETS = ["Off", "Wood", "Metal", "Glass"]


class AureonoiseApp:
    """Main GUI application."""

    def __init__(self, sample_rate: float = 44100.0, block_size: int = 512):
        if dpg is None:
            raise ImportError("DearPyGui required: pip install dearpygui")

        self.audio = AudioEngine(sample_rate, block_size)
        self.running = False
        self._meter_thread: Optional[threading.Thread] = None
        self._preset_bank = PresetBank()
        self._log_sliders: dict[str, tuple[float, float]] = {}

    def run(self):
        """Run the application."""
        dpg.create_context()

        # Theme setup
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

        # Per-button themes for therapeutic presets (created once, applied later)
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

        # Main window
        with dpg.window(label="aureonoise", tag="main", no_title_bar=True):

            # Header
            dpg.add_text("aureonoise", color=COLORS["accent"])
            dpg.add_text("phi-based granular texture generator", color=COLORS["text_dim"])
            dpg.add_spacer(height=10)

            # Transport
            with dpg.group(horizontal=True):
                dpg.add_button(label="  Play", tag="btn_play", callback=self._on_play)
                dpg.add_button(label="  Stop", tag="btn_stop", callback=self._on_stop)
                dpg.add_button(label="  Reset", callback=self._on_reset)
                dpg.add_spacer(width=20)
                dpg.add_text("", tag="status", color=COLORS["text_dim"])

            dpg.add_spacer(height=10)

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

            # Parameters in tabs
            with dpg.tab_bar():

                # ── Timing tab ──────────────────────────────────────
                with dpg.tab(label="Timing"):
                    dpg.add_spacer(height=5)
                    self._slider("rate", "Rate (Hz)", 0.1, 60.0, 8.0, log=True)
                    self._slider("baselen_ms", "Base Length (ms)", 10.0, 2000.0, 120.0, log=True)
                    self._slider("len_phi", "Length phi Spread", 0.0, 1.0, 0.8)

                # ── Spatial tab (expanded) ──────────────────────────
                with dpg.tab(label="Spatial"):
                    dpg.add_spacer(height=5)
                    dpg.add_text("Stereo Field", color=COLORS["section"])
                    self._slider("width", "Stereo Width", 0.0, 2.0, 1.0)
                    self._slider("itd_us", "ITD (us)", 0.0, 800.0, 600.0)
                    self._slider("ild_db", "ILD (dB)", 0.0, 12.0, 6.0)
                    dpg.add_separator()
                    self._slider("hemis_coupling", "Hemisphere Coupling", 0.0, 1.0, 0.6)
                    self._slider("spat_ipd", "IPD Amount", 0.0, 1.0, 0.6)
                    self._slider("spat_shadow", "Head Shadow", 0.0, 1.0, 0.7)

                    dpg.add_separator()
                    dpg.add_spacer(height=5)
                    dpg.add_text("Externalization", color=COLORS["section"])
                    self._slider("externalization", "Externalization", 0.0, 1.0, 0.0)
                    self._slider("phi_distance", "Phi Distance (m)", 0.0, 10.0, 1.5)
                    self._slider("phi_elev", "Phi Elevation (deg)", -90.0, 90.0, 0.0)

                    dpg.add_separator()
                    dpg.add_spacer(height=5)
                    dpg.add_text("Phi-Pan + Bilateral", color=COLORS["section"])
                    dpg.add_checkbox(
                        label="Phi-Pan", default_value=False,
                        callback=lambda s, a, u: self._set_param(u, a),
                        user_data="phi_pan")
                    dpg.add_checkbox(
                        label="Bilateral On", default_value=False, tag="cb_bilateral_on",
                        callback=self._on_bilateral_toggle)
                    self._slider("bilateral_rate", "Bilateral Rate (Hz)", 0.5, 2.0, 1.0)
                    self._slider("bilateral_amount", "Bilateral Amount", 0.0, 1.0, 0.8)

                # ── Envelope tab ────────────────────────────────────
                with dpg.tab(label="Envelope"):
                    dpg.add_spacer(height=5)
                    self._slider("env_attack", "Attack", 0.01, 1.0, 0.18)
                    self._slider("env_decay", "Decay", 0.01, 1.0, 0.28)
                    self._slider("env_sustain", "Sustain", 0.0, 1.0, 0.55)
                    self._slider("env_release", "Release", 0.01, 1.0, 0.30)

                # ── Timbre tab ──────────────────────────────────────
                with dpg.tab(label="Timbre"):
                    dpg.add_spacer(height=5)
                    dpg.add_text("Noise color controlled in Noise tab", color=COLORS["text_dim"])
                    dpg.add_separator()
                    self._slider("vhs_wow", "VHS Wow", 0.0, 1.0, 0.35)
                    self._slider("vhs_flutter", "VHS Flutter", 0.0, 1.0, 0.25)
                    self._slider("glitch_mix", "Glitch Mix", 0.0, 1.0, 0.5)
                    dpg.add_separator()
                    self._slider("srcrush_amt", "Sample Rate Crush", 0.0, 1.0, 0.2)
                    self._slider("bitcrush_amt", "Bit Crush", 0.0, 1.0, 0.15)

                # ── Noise tab (extended 6 modes) ────────────────────
                with dpg.tab(label="Noise"):
                    dpg.add_spacer(height=5)
                    dpg.add_text("Noise Mode", color=COLORS["section"])
                    dpg.add_radio_button(
                        NOISE_MODES,
                        default_value="Pink",
                        tag="radio_noise_mode",
                        callback=self._on_noise_mode_change,
                        horizontal=True,
                    )

                    dpg.add_separator()
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

                # ── Dialogue tab ────────────────────────────────────
                with dpg.tab(label="Dialogue"):
                    dpg.add_spacer(height=5)
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

                    # Coherence meter
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

                # ── Modal tab ───────────────────────────────────────
                with dpg.tab(label="Modal"):
                    dpg.add_spacer(height=5)
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

                # ── Binaural tab ───────────────────────────────────
                with dpg.tab(label="Binaural"):
                    dpg.add_spacer(height=5)
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

                # ── Tinnitus / Notch tab ──────────────────────────
                with dpg.tab(label="Tinnitus"):
                    dpg.add_spacer(height=5)
                    dpg.add_text("Tinnitus Notch Filter", color=COLORS["section"])
                    dpg.add_text(
                        "4th-order Butterworth notch at your tinnitus frequency. "
                        "Set to 0 to disable.",
                        color=COLORS["text_dim"], wrap=520)
                    dpg.add_spacer(height=5)
                    self._slider("tinnitus_notch_hz", "Center Freq (Hz)", 0.0, 12000.0, 0.0)
                    self._slider("tinnitus_notch_q", "Q Factor", 1.0, 20.0, 6.0)

                    dpg.add_separator()
                    dpg.add_spacer(height=5)
                    dpg.add_text("Spectral Slope", color=COLORS["section"])
                    dpg.add_text(
                        "Continuous tilt: 0=white, -1=pink, -2=brown.",
                        color=COLORS["text_dim"], wrap=520)
                    self._slider("noise_slope", "Slope (dB/oct)", -2.5, 0.5, -1.0)

                # ── Feedback / Spatial tab ─────────────────────────
                with dpg.tab(label="Feedback"):
                    dpg.add_spacer(height=5)
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

                    dpg.add_separator()
                    dpg.add_spacer(height=5)
                    dpg.add_text("Polyrhythm Clock", color=COLORS["section"])
                    dpg.add_checkbox(
                        label="Polyrhythm On", default_value=False, tag="cb_polyrhythm_on",
                        callback=self._on_polyrhythm_toggle)
                    self._slider_int("polyrhythm_p", "P (left)", 2, 8, 3)
                    self._slider_int("polyrhythm_q", "Q (right)", 2, 8, 2)
                    self._slider("polyrhythm_rate", "Base Rate (Hz)", 0.1, 3.0, 0.5)
                    self._slider("polyrhythm_amount", "Amount", 0.0, 1.0, 0.5)

                    dpg.add_separator()
                    dpg.add_spacer(height=5)
                    dpg.add_text("Room Reverb", color=COLORS["section"])
                    self._slider("room_mix", "Room Mix", 0.0, 0.5, 0.0)

                # ── Stochastic tab ──────────────────────────────────
                with dpg.tab(label="Stochastic"):
                    dpg.add_spacer(height=5)
                    with dpg.group(horizontal=True):
                        dpg.add_checkbox(label="Thermo", default_value=True,
                            callback=lambda s, a, u: self._set_param(u, a),
                            user_data="thermo")
                        dpg.add_checkbox(label="Lattice", default_value=True,
                            callback=lambda s, a, u: self._set_param(u, a),
                            user_data="lattice")
                        dpg.add_checkbox(label="Burst", default_value=True,
                            callback=lambda s, a, u: self._set_param(u, a),
                            user_data="burst")
                    dpg.add_separator()
                    self._slider("temperature", "Temperature", 0.0, 1.0, 0.45)
                    self._slider("lat_rate", "Lattice Rate", 1.0, 2000.0, 250.0, log=True)
                    self._slider("lat_eps", "Lattice e", 0.01, 0.5, INV_PHI_CU)
                    self._slider("lat_gamma", "Lattice g", 0.5, 3.0, PHI)
                    self._slider("lat_sigma", "Lattice s", 0.01, 0.3, 0.06)

                # ── Therapeutic tab ──────────────────────────────────
                with dpg.tab(label="Therapeutic"):
                    dpg.add_spacer(height=5)
                    dpg.add_text("Evidence-Based Presets", color=COLORS["section"])
                    dpg.add_text(
                        "Each preset configures noise, spatial, dialogue and bilateral "
                        "parameters for a specific therapeutic context.",
                        color=COLORS["text_dim"], wrap=520)
                    dpg.add_spacer(height=10)

                    # 3x2 button grid
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

                # ── System tab ──────────────────────────────────────
                with dpg.tab(label="System"):
                    dpg.add_spacer(height=5)

                    # Audio device selector
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

    # ── Widget helpers ──────────────────────────────────────────────

    def _slider(self, param: str, label: str, vmin: float, vmax: float,
                default: float, log: bool = False):
        """Create a parameter slider."""
        fmt = "%.3f" if vmax <= 2.0 else "%.1f"

        if log:
            self._log_sliders[param] = (vmin, vmax)

            def log_cb(sender, app_data, user_data):
                p, lo, hi = user_data
                v = lo * ((hi / lo) ** app_data)
                self._set_param(p, v)
                dpg.set_value(f"val_{p}", f"{v:.2f}")

            with dpg.group(horizontal=True):
                dpg.add_slider_float(
                    tag=f"sl_{param}",
                    label="",
                    default_value=0.5,
                    min_value=0.0,
                    max_value=1.0,
                    callback=log_cb,
                    user_data=(param, vmin, vmax),
                    width=200,
                )
                dpg.add_text(f"{default:.2f}", tag=f"val_{param}", color=COLORS["text_dim"])
                dpg.add_text(label, color=COLORS["text"])
        else:
            def lin_cb(sender, app_data, user_data):
                self._set_param(user_data, app_data)

            dpg.add_slider_float(
                tag=f"sl_{param}",
                label=label,
                default_value=default,
                min_value=vmin,
                max_value=vmax,
                callback=lin_cb,
                user_data=param,
                format=fmt,
                width=200,
            )

    def _slider_int(self, param: str, label: str, vmin: int, vmax: int,
                    default: int):
        """Create an integer parameter slider."""
        def int_cb(sender, app_data, user_data):
            self._set_param(user_data, int(app_data))

        dpg.add_slider_int(
            tag=f"sl_{param}",
            label=label,
            default_value=default,
            min_value=vmin,
            max_value=vmax,
            callback=int_cb,
            user_data=param,
            width=200,
        )

    # ── Parameter dispatch ──────────────────────────────────────────

    def _set_param(self, name: str, value):
        """Set parameter on audio engine."""
        self.audio.set_param(name, value)

    # ── Mutual exclusion / dependency toggle handlers ────────────────

    def _on_binaural_toggle(self, sender, value, user_data):
        self._set_param("binaural_on", value)
        if value:
            self._set_param("isochronic_on", False)
            _safe_set("cb_isochronic_on", False)
        self._update_feature_summary()

    def _on_isochronic_toggle(self, sender, value, user_data):
        self._set_param("isochronic_on", value)
        if value:
            self._set_param("binaural_on", False)
            _safe_set("cb_binaural_on", False)
        self._update_feature_summary()

    def _on_dialogue_toggle(self, sender, value, user_data):
        self._set_param("dialogue_on", value)
        if not value:
            self._set_param("feedback_on", False)
            _safe_set("cb_feedback_on", False)
            self._set_param("coherence_spatial", False)
            _safe_set("cb_coherence_spatial", False)
        _safe_enable("cb_feedback_on", value)
        _safe_enable("cb_coherence_spatial", value)
        self._update_feature_summary()

    def _on_bilateral_toggle(self, sender, value, user_data):
        self._set_param("bilateral_on", value)
        if not value:
            self._set_param("bilateral_nesting", False)
            _safe_set("cb_bilateral_nesting", False)
        _safe_enable("cb_bilateral_nesting", value)
        self._update_feature_summary()

    def _on_feedback_toggle(self, sender, value, user_data):
        self._set_param("feedback_on", value)
        self._update_feature_summary()

    def _on_coherence_spatial_toggle(self, sender, value, user_data):
        self._set_param("coherence_spatial", value)
        self._update_feature_summary()

    def _on_bilateral_nesting_toggle(self, sender, value, user_data):
        self._set_param("bilateral_nesting", value)
        self._update_feature_summary()

    def _on_polyrhythm_toggle(self, sender, value, user_data):
        self._set_param("polyrhythm_on", value)
        self._update_feature_summary()

    def _on_modal_toggle(self, sender, value, user_data):
        self._set_param("modal_on", value)
        self._update_feature_summary()

    def _update_feature_summary(self):
        """Update the active features summary text."""
        parts = []
        try:
            if dpg.get_value("cb_bilateral_on"):
                try:
                    rate = self.audio.engine.params.bilateral_rate
                    parts.append(f"Bilateral {rate}Hz")
                except Exception:
                    parts.append("Bilateral")
            if dpg.get_value("cb_binaural_on"):
                try:
                    beat = dpg.get_value("sl_binaural_beat_hz")
                    parts.append(f"Binaural {beat:.0f}Hz")
                except Exception:
                    parts.append("Binaural")
            if dpg.get_value("cb_isochronic_on"):
                parts.append("Isochronic")
            if dpg.get_value("cb_dialogue_on"):
                parts.append("Dialogue")
            if dpg.get_value("cb_feedback_on"):
                parts.append("Feedback")
            if dpg.get_value("cb_coherence_spatial"):
                parts.append("CohSpatial")
            if dpg.get_value("cb_bilateral_nesting"):
                parts.append("Nesting")
            if dpg.get_value("cb_polyrhythm_on"):
                parts.append("Polyrhythm")
            if dpg.get_value("cb_modal_on"):
                parts.append("Modal")
        except Exception:
            pass
        summary = " + ".join(parts) if parts else "No active features"
        _safe_set("txt_feature_summary", summary)

    # ── Noise mode visibility toggle ────────────────────────────────

    def _on_noise_mode_change(self, sender, app_data):
        """Handle noise mode radio button change."""
        idx = NOISE_MODES.index(app_data) if app_data in NOISE_MODES else 1
        self._set_param("noise_mode", idx)

        # Show/hide conditional parameter groups
        dpg.configure_item("grp_aureo", show=(idx == 3))
        dpg.configure_item("grp_quantum", show=(idx == 4))
        dpg.configure_item("grp_velvet", show=(idx == 5))

    # ── Preset application ──────────────────────────────────────────

    def _make_preset_callback(self, name: str):
        """Return a callback that loads a therapeutic preset."""
        def cb():
            self._apply_preset(name)
        return cb

    def _apply_preset(self, name: str):
        """Apply a named preset from the bank, sync all GUI controls."""
        factory_key, _ = FACTORY.get(name, (name, None))
        preset = self._preset_bank.get(factory_key)
        if preset is None:
            return

        # Apply every param from preset dict
        for k, v in preset.params.items():
            self._set_param(k, v)

        # Sync GUI sliders and controls to reflect new values
        self._sync_gui(preset.params)

        dpg.set_value("txt_preset_active", f"Active: {name}")

    def _sync_gui(self, params: dict):
        """Synchronize all GUI widgets to match a parameter dict.

        Best-effort: skips any widget that does not exist yet (tags created
        lazily on first tab visit in DPG).
        """
        for k, v in params.items():
            sl_tag = f"sl_{k}"
            val_tag = f"val_{k}"

            # Log sliders: reverse-compute 0-1 normalized value
            if k in self._log_sliders:
                lo, hi = self._log_sliders[k]
                try:
                    clamped = max(lo, min(float(v), hi))
                    norm = math.log(clamped / lo) / math.log(hi / lo)
                    dpg.set_value(sl_tag, max(0.0, min(norm, 1.0)))
                except Exception:
                    pass
            else:
                # Regular sliders (float and int)
                try:
                    dpg.set_value(sl_tag, v)
                except Exception:
                    pass

            # Log-slider companion text
            try:
                dpg.set_value(val_tag, f"{v:.2f}")
            except Exception:
                pass

        # Checkboxes
        _safe_set("cb_dialogue_on", params.get("dialogue_on", True))
        _safe_set("cb_bilateral_on", params.get("bilateral_on", False))
        _safe_set("cb_modal_on", params.get("modal_on", False))
        _safe_set("cb_binaural_on", params.get("binaural_on", False))
        _safe_set("cb_isochronic_on", params.get("isochronic_on", False))
        _safe_set("cb_feedback_on", params.get("feedback_on", False))
        _safe_set("cb_bilateral_nesting", params.get("bilateral_nesting", False))
        _safe_set("cb_coherence_spatial", params.get("coherence_spatial", False))
        _safe_set("cb_polyrhythm_on", params.get("polyrhythm_on", False))

        # Noise mode radio
        nm = params.get("noise_mode", 1)
        if 0 <= nm < len(NOISE_MODES):
            _safe_set("radio_noise_mode", NOISE_MODES[nm])
            # Trigger visibility
            _safe_configure("grp_aureo", show=(nm == 3))
            _safe_configure("grp_quantum", show=(nm == 4))
            _safe_configure("grp_velvet", show=(nm == 5))

        # Modal preset combo
        mp = params.get("modal_preset", 1)
        if 0 <= mp < len(MODAL_PRESETS):
            _safe_set("combo_modal_preset", MODAL_PRESETS[mp])

        # Enforce mutual exclusion and dependencies
        dialogue_on = params.get("dialogue_on", True)
        bilateral_on = params.get("bilateral_on", False)
        _safe_enable("cb_feedback_on", dialogue_on)
        _safe_enable("cb_coherence_spatial", dialogue_on)
        _safe_enable("cb_bilateral_nesting", bilateral_on)
        self._update_feature_summary()

    # ── Audio device ─────────────────────────────────────────────────

    @staticmethod
    def _get_output_devices():
        """Return list of (index, name) for output-capable devices."""
        if sd is None:
            return []
        devices = sd.query_devices()
        return [
            (i, d["name"])
            for i, d in enumerate(devices)
            if d["max_output_channels"] > 0
        ]

    def _on_device_change(self, sender, app_data):
        """Handle audio device combo change."""
        try:
            dev_idx = int(app_data.split(":")[0])
        except (ValueError, IndexError):
            return
        was_running = self.audio.is_running()
        if was_running:
            self.audio.stop()
        self.audio.output_device = dev_idx
        if was_running:
            self.audio.start()
            dpg.set_value("status", "Playing")

    # ── Transport ───────────────────────────────────────────────────

    def _on_play(self):
        """Play button handler."""
        if not self.audio.is_running():
            self.audio.start()
            dpg.set_value("status", "Playing")

    def _on_stop(self):
        """Stop button handler."""
        self.audio.stop()
        dpg.set_value("status", "Stopped")
        dpg.set_value("meter_l", 0)
        dpg.set_value("meter_r", 0)

    def _on_reset(self):
        """Reset button handler."""
        self.audio.reset()
        dpg.set_value("status", "Reset — DSP state cleared")
        # Clear dialogue metrics display
        _safe_set("meter_coherence", 0)
        _safe_configure("meter_coherence", overlay="0.000")
        _safe_set("txt_handshake_count", "0")
        _safe_set("txt_handshake_ratio", "0.000")
        _safe_set("txt_coherence_mean", "0.000")

    # ── Meter / coherence update thread ─────────────────────────────

    def _start_meter_thread(self):
        """Start meter update thread."""
        self.running = True

        def update_meters():
            while self.running:
                if self.audio.is_running():
                    # Audio level meters
                    pl, pr, rl, rr = self.audio.get_meters()
                    try:
                        dpg.set_value("meter_l", min(pl, 1.0))
                        dpg.set_value("meter_r", min(pr, 1.0))
                    except Exception:
                        pass

                    # Dialogue engine metrics (read-only from Rust engine)
                    try:
                        eng = self.audio.engine
                        coh = eng.coherence()
                        dpg.set_value("meter_coherence", min(coh, 1.0))
                        dpg.configure_item("meter_coherence",
                                           overlay=f"{coh:.3f}")
                        dpg.set_value("txt_handshake_count",
                                      str(eng.handshake_count()))
                        dpg.set_value("txt_handshake_ratio",
                                      f"{eng.handshake_ratio():.3f}")
                        dpg.set_value("txt_coherence_mean",
                                      f"{eng.coherence_mean():.3f}")
                    except Exception:
                        pass

                time.sleep(0.033)  # ~30 fps

        self._meter_thread = threading.Thread(target=update_meters, daemon=True)
        self._meter_thread.start()


# ── Module-level helpers (avoid closure capture issues) ─────────────

def _safe_set(tag: str, value):
    """Set DPG item value, silently ignoring missing tags."""
    try:
        dpg.set_value(tag, value)
    except Exception:
        pass


def _safe_configure(tag: str, **kwargs):
    """Configure DPG item, silently ignoring missing tags."""
    try:
        dpg.configure_item(tag, **kwargs)
    except Exception:
        pass


def _safe_enable(tag: str, enabled: bool):
    """Enable or disable a DPG item."""
    try:
        dpg.configure_item(tag, enabled=enabled)
    except Exception:
        pass


def main():
    """Main entry point."""
    app = AureonoiseApp()
    app.run()


if __name__ == "__main__":
    main()
