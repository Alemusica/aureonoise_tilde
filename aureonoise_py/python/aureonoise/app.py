"""
aureonoise - GUI Application

DearPyGui-based interface for real-time control.
"""

import sys
from typing import Optional
import threading
import time

try:
    import dearpygui.dearpygui as dpg
except ImportError:
    dpg = None

from aureonoise import Params, PHI, INV_PHI, INV_PHI_SQ, INV_PHI_CU
from aureonoise.audio import AudioEngine


# φ-based color palette
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
}


class AureonoiseApp:
    """Main GUI application."""
    
    def __init__(self, sample_rate: float = 44100.0, block_size: int = 512):
        if dpg is None:
            raise ImportError("DearPyGui required: pip install dearpygui")
        
        self.audio = AudioEngine(sample_rate, block_size)
        self.running = False
        self._meter_thread: Optional[threading.Thread] = None
        
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
                dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
                dpg.add_theme_style(dpg.mvStyleVar_GrabRounding, 4)
                dpg.add_theme_style(dpg.mvStyleVar_WindowPadding, 12, 12)
        
        dpg.bind_theme(global_theme)
        
        # Main window
        with dpg.window(label="aureonoise", tag="main", no_title_bar=True):
            
            # Header
            dpg.add_text("aureonoise", color=COLORS["accent"])
            dpg.add_text("φ-based granular texture generator", color=COLORS["text_dim"])
            dpg.add_spacer(height=10)
            
            # Transport
            with dpg.group(horizontal=True):
                dpg.add_button(label="▶ Play", tag="btn_play", callback=self._on_play)
                dpg.add_button(label="■ Stop", tag="btn_stop", callback=self._on_stop)
                dpg.add_button(label="↺ Reset", callback=self._on_reset)
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
                
                # Timing tab
                with dpg.tab(label="Timing"):
                    dpg.add_spacer(height=5)
                    self._slider("rate", "Rate (Hz)", 0.1, 60.0, 8.0, log=True)
                    self._slider("baselen_ms", "Base Length (ms)", 10.0, 2000.0, 120.0, log=True)
                    self._slider("len_phi", "Length φ Spread", 0.0, 1.0, 0.8)
                
                # Spatial tab
                with dpg.tab(label="Spatial"):
                    dpg.add_spacer(height=5)
                    self._slider("width", "Stereo Width", 0.0, 2.0, 1.0)
                    self._slider("itd_us", "ITD (μs)", 0.0, 800.0, 600.0)
                    self._slider("ild_db", "ILD (dB)", 0.0, 12.0, 6.0)
                    dpg.add_separator()
                    self._slider("hemis_coupling", "Hemisphere Coupling", 0.0, 1.0, 0.6)
                    self._slider("spat_ipd", "IPD Amount", 0.0, 1.0, 0.6)
                    self._slider("spat_shadow", "Head Shadow", 0.0, 1.0, 0.7)
                
                # Envelope tab
                with dpg.tab(label="Envelope"):
                    dpg.add_spacer(height=5)
                    self._slider("env_attack", "Attack", 0.01, 1.0, 0.18)
                    self._slider("env_decay", "Decay", 0.01, 1.0, 0.28)
                    self._slider("env_sustain", "Sustain", 0.0, 1.0, 0.55)
                    self._slider("env_release", "Release", 0.01, 1.0, 0.30)
                
                # Timbre tab  
                with dpg.tab(label="Timbre"):
                    dpg.add_spacer(height=5)
                    dpg.add_combo(
                        ["White", "Pink", "Brown"],
                        default_value="Pink",
                        label="Noise Color",
                        callback=lambda s, a: self._set_param("noise_color", 
                            {"White": 0, "Pink": 1, "Brown": 2}.get(a, 1)),
                        width=150,
                    )
                    self._slider("color_amt", "Color Amount", 0.0, 1.0, 0.65)
                    dpg.add_separator()
                    self._slider("vhs_wow", "VHS Wow", 0.0, 1.0, 0.35)
                    self._slider("vhs_flutter", "VHS Flutter", 0.0, 1.0, 0.25)
                    self._slider("glitch_mix", "Glitch Mix", 0.0, 1.0, 0.5)
                    dpg.add_separator()
                    self._slider("srcrush_amt", "Sample Rate Crush", 0.0, 1.0, 0.2)
                    self._slider("bitcrush_amt", "Bit Crush", 0.0, 1.0, 0.15)
                
                # Stochastic tab
                with dpg.tab(label="Stochastic"):
                    dpg.add_spacer(height=5)
                    with dpg.group(horizontal=True):
                        dpg.add_checkbox(label="Thermo", default_value=True, 
                            callback=lambda s, a: self._set_param("thermo", a))
                        dpg.add_checkbox(label="Lattice", default_value=True,
                            callback=lambda s, a: self._set_param("lattice", a))
                        dpg.add_checkbox(label="Burst", default_value=True,
                            callback=lambda s, a: self._set_param("burst", a))
                    dpg.add_separator()
                    self._slider("temperature", "Temperature", 0.0, 1.0, 0.45)
                    self._slider("lat_rate", "Lattice Rate", 1.0, 2000.0, 250.0, log=True)
                    self._slider("lat_eps", "Lattice ε", 0.01, 0.5, INV_PHI_CU)
                    self._slider("lat_gamma", "Lattice γ", 0.5, 3.0, PHI)
                    self._slider("lat_sigma", "Lattice σ", 0.01, 0.3, 0.06)
                
                # System tab
                with dpg.tab(label="System"):
                    dpg.add_spacer(height=5)
                    dpg.add_input_int(
                        label="Seed",
                        default_value=20251010,
                        callback=lambda s, a: self._set_param("seed", a),
                        width=150,
                    )
                    dpg.add_spacer(height=10)
                    dpg.add_text(f"Sample Rate: {self.audio.sample_rate} Hz", color=COLORS["text_dim"])
                    dpg.add_text(f"Block Size: {self.audio.block_size}", color=COLORS["text_dim"])
                    dpg.add_spacer(height=10)
                    dpg.add_text("φ Constants:", color=COLORS["accent"])
                    dpg.add_text(f"  φ = {PHI:.6f}", color=COLORS["text_dim"])
                    dpg.add_text(f"  1/φ = {INV_PHI:.6f}", color=COLORS["text_dim"])
                    dpg.add_text(f"  1/φ² = {INV_PHI_SQ:.6f}", color=COLORS["text_dim"])
                    dpg.add_text(f"  1/φ³ = {INV_PHI_CU:.6f}", color=COLORS["text_dim"])
        
        # Viewport setup
        dpg.create_viewport(
            title="aureonoise",
            width=600,
            height=700,
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
    
    def _slider(self, param: str, label: str, vmin: float, vmax: float, default: float, log: bool = False):
        """Create a parameter slider."""
        fmt = "%.3f" if vmax <= 2.0 else "%.1f"
        
        if log:
            # Logarithmic slider
            def callback(sender, app_data):
                # Map 0-1 to log range
                t = app_data
                v = vmin * ((vmax/vmin) ** t)
                self._set_param(param, v)
                dpg.set_value(f"val_{param}", f"{v:.2f}")
            
            t0 = (default - vmin) / (vmax - vmin) if not log else 0.5
            with dpg.group(horizontal=True):
                dpg.add_slider_float(
                    tag=f"sl_{param}",
                    label="",
                    default_value=0.5,
                    min_value=0.0,
                    max_value=1.0,
                    callback=callback,
                    width=200,
                )
                dpg.add_text(f"{default:.2f}", tag=f"val_{param}", color=COLORS["text_dim"])
                dpg.add_text(label, color=COLORS["text"])
        else:
            def callback(sender, app_data):
                self._set_param(param, app_data)
            
            dpg.add_slider_float(
                tag=f"sl_{param}",
                label=label,
                default_value=default,
                min_value=vmin,
                max_value=vmax,
                callback=callback,
                format=fmt,
                width=200,
            )
    
    def _set_param(self, name: str, value):
        """Set parameter on audio engine."""
        self.audio.set_param(name, value)
    
    def _on_play(self):
        """Play button handler."""
        if not self.audio.is_running():
            self.audio.start()
            dpg.set_value("status", "● Playing")
    
    def _on_stop(self):
        """Stop button handler."""
        self.audio.stop()
        dpg.set_value("status", "○ Stopped")
        dpg.set_value("meter_l", 0)
        dpg.set_value("meter_r", 0)
    
    def _on_reset(self):
        """Reset button handler."""
        self.audio.reset()
        dpg.set_value("status", "↺ Reset")
    
    def _start_meter_thread(self):
        """Start meter update thread."""
        self.running = True
        
        def update_meters():
            while self.running:
                if self.audio.is_running():
                    pl, pr, rl, rr = self.audio.get_meters()
                    try:
                        dpg.set_value("meter_l", min(pl, 1.0))
                        dpg.set_value("meter_r", min(pr, 1.0))
                    except:
                        pass
                time.sleep(0.033)  # ~30 fps
        
        self._meter_thread = threading.Thread(target=update_meters, daemon=True)
        self._meter_thread.start()


def main():
    """Main entry point."""
    app = AureonoiseApp()
    app.run()


if __name__ == "__main__":
    main()
