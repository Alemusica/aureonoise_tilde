"""
aureonoise - Preset management

Save and load parameter presets.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict

from aureonoise import Params, PHI, INV_PHI, INV_PHI_SQ, INV_PHI_CU


@dataclass
class Preset:
    """A named preset with parameters."""
    name: str
    description: str
    author: str
    params: Dict
    
    def to_params(self) -> Params:
        """Convert to Params object."""
        p = Params()
        for k, v in self.params.items():
            if hasattr(p, k):
                setattr(p, k, v)
        return p
    
    @classmethod
    def from_params(cls, name: str, params: Params, description: str = "", author: str = ""):
        """Create from Params object."""
        d = {
            # Timing
            "rate": params.rate,
            "baselen_ms": params.baselen_ms,
            "len_phi": params.len_phi,
            # Spatial
            "width": params.width,
            "itd_us": params.itd_us,
            "ild_db": params.ild_db,
            "hemis_coupling": params.hemis_coupling,
            "spat_min_deg": params.spat_min_deg,
            "spat_min_ms": params.spat_min_ms,
            "spat_ipd": params.spat_ipd,
            "spat_shadow": params.spat_shadow,
            # Envelope
            "env_attack": params.env_attack,
            "env_decay": params.env_decay,
            "env_sustain": params.env_sustain,
            "env_release": params.env_release,
            # Timbre
            "noise_color": params.noise_color,
            "color_amt": params.color_amt,
            "vhs_wow": params.vhs_wow,
            "vhs_flutter": params.vhs_flutter,
            "glitch_mix": params.glitch_mix,
            "srcrush_amt": params.srcrush_amt,
            "bitcrush_amt": params.bitcrush_amt,
            # Stochastic
            "thermo": params.thermo,
            "lattice": params.lattice,
            "burst": params.burst,
            "temperature": params.temperature,
            "lat_rate": params.lat_rate,
            "lat_eps": params.lat_eps,
            "lat_gamma": params.lat_gamma,
            "lat_sigma": params.lat_sigma,
            "burst_floor": params.burst_floor,
            "burst_phi_mix": params.burst_phi_mix,
            # Externalization
            "externalization": params.externalization,
            # Dialogue
            "dialogue_on": params.dialogue_on,
            "dialogue_strength": params.dialogue_strength,
            "dialogue_memory": params.dialogue_memory,
            "dialogue_phi_mix": params.dialogue_phi_mix,
            # Phi-Pan + Bilateral
            "phi_pan": params.phi_pan,
            "bilateral_on": params.bilateral_on,
            "bilateral_rate": params.bilateral_rate,
            "bilateral_amount": params.bilateral_amount,
            # Noise (extended modes)
            "noise_mode": params.noise_mode,
            "aureo_decay": params.aureo_decay,
            "aureo_stride": params.aureo_stride,
            "aureo_harmonics": params.aureo_harmonics,
            "quantum_detail": params.quantum_detail,
            "quantum_base": params.quantum_base,
            "velvet_density": params.velvet_density,
            # Modal
            "modal_on": params.modal_on,
            "modal_mix": params.modal_mix,
            "modal_decay": params.modal_decay,
            "modal_preset": params.modal_preset,
            "modal_mirror": params.modal_mirror,
            "modal_feedback": params.modal_feedback,
            # Binaural beat
            "binaural_on": params.binaural_on,
            "binaural_carrier_hz": params.binaural_carrier_hz,
            "binaural_beat_hz": params.binaural_beat_hz,
            "binaural_level": params.binaural_level,
            # Isochronic
            "isochronic_on": params.isochronic_on,
            "isochronic_carrier_hz": params.isochronic_carrier_hz,
            "isochronic_rate_hz": params.isochronic_rate_hz,
            "isochronic_duty": params.isochronic_duty,
            "isochronic_level": params.isochronic_level,
            # Spectral slope
            "noise_slope": params.noise_slope,
            # Tinnitus notch
            "tinnitus_notch_hz": params.tinnitus_notch_hz,
            "tinnitus_notch_q": params.tinnitus_notch_q,
            # Feedback / nesting
            "bilateral_nesting": params.bilateral_nesting,
            "feedback_on": params.feedback_on,
            "temp_ramp_sec": params.temp_ramp_sec,
            # Phi model
            "phi_distance": params.phi_distance,
            "phi_elev": params.phi_elev,
            # System
            "seed": params.seed,
        }
        return cls(name=name, description=description, author=author, params=d)
    
    def to_json(self) -> str:
        """Serialize to JSON."""
        return json.dumps(asdict(self), indent=2)
    
    @classmethod
    def from_json(cls, data: str) -> "Preset":
        """Deserialize from JSON."""
        d = json.loads(data)
        return cls(**d)
    
    def save(self, path: Path):
        """Save to file."""
        path.write_text(self.to_json())
    
    @classmethod
    def load(cls, path: Path) -> "Preset":
        """Load from file."""
        return cls.from_json(path.read_text())


class PresetBank:
    """Collection of presets."""
    
    def __init__(self, directory: Optional[Path] = None):
        self.directory = directory or (Path.home() / ".aureonoise" / "presets")
        self.presets: Dict[str, Preset] = {}
        self._init_defaults()
    
    def _init_defaults(self):
        """Initialize factory presets."""
        self.presets["Default"] = Preset(
            name="Default",
            description="Default φ-balanced texture",
            author="aureonoise",
            params={
                "rate": 8.0,
                "baselen_ms": 120.0,
                "len_phi": 0.8,
                "width": 1.0,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.6,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.6,
                "spat_shadow": 0.7,
                "env_attack": 0.18,
                "env_decay": 0.28,
                "env_sustain": 0.55,
                "env_release": 0.30,
                "noise_color": 1,
                "color_amt": 0.65,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.5,
                "srcrush_amt": 0.2,
                "bitcrush_amt": 0.15,
                "thermo": True,
                "lattice": True,
                "burst": True,
                "temperature": 0.45,
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "seed": 20251010,
            }
        )
        
        self.presets["Calm Rain"] = Preset(
            name="Calm Rain",
            description="Gentle rain-like texture",
            author="aureonoise",
            params={
                "rate": 12.0,
                "baselen_ms": 80.0,
                "len_phi": 0.5,
                "width": 1.2,
                "itd_us": 500.0,
                "ild_db": 4.0,
                "hemis_coupling": 0.4,
                "spat_min_deg": 15.0,
                "spat_min_ms": 40.0,
                "spat_ipd": 0.5,
                "spat_shadow": 0.6,
                "env_attack": 0.25,
                "env_decay": 0.35,
                "env_sustain": 0.4,
                "env_release": 0.4,
                "noise_color": 2,  # Brown
                "color_amt": 0.7,
                "vhs_wow": 0.1,
                "vhs_flutter": 0.1,
                "glitch_mix": 0.1,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                "thermo": True,
                "lattice": True,
                "burst": False,
                "temperature": 0.3,
                "lat_rate": 150.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.04,
                "seed": 20251010,
            }
        )
        
        self.presets["Digital Storm"] = Preset(
            name="Digital Storm",
            description="Aggressive glitch texture",
            author="aureonoise",
            params={
                "rate": 20.0,
                "baselen_ms": 60.0,
                "len_phi": 1.0,
                "width": 1.5,
                "itd_us": 700.0,
                "ild_db": 8.0,
                "hemis_coupling": 0.8,
                "spat_min_deg": 10.0,
                "spat_min_ms": 25.0,
                "spat_ipd": 0.8,
                "spat_shadow": 0.8,
                "env_attack": 0.05,
                "env_decay": 0.15,
                "env_sustain": 0.7,
                "env_release": 0.15,
                "noise_color": 0,  # White
                "color_amt": 0.5,
                "vhs_wow": 0.6,
                "vhs_flutter": 0.5,
                "glitch_mix": 0.9,
                "srcrush_amt": 0.5,
                "bitcrush_amt": 0.4,
                "thermo": True,
                "lattice": True,
                "burst": True,
                "temperature": 0.8,
                "lat_rate": 400.0,
                "lat_eps": INV_PHI_SQ,
                "lat_gamma": PHI * 1.2,
                "lat_sigma": 0.1,
                "seed": 20251010,
            }
        )
        
        self.presets["Vinyl Crackle"] = Preset(
            name="Vinyl Crackle",
            description="Lo-fi vinyl simulation",
            author="aureonoise",
            params={
                "rate": 15.0,
                "baselen_ms": 40.0,
                "len_phi": 0.6,
                "width": 0.8,
                "itd_us": 400.0,
                "ild_db": 3.0,
                "hemis_coupling": 0.5,
                "spat_min_deg": 8.0,
                "spat_min_ms": 30.0,
                "spat_ipd": 0.4,
                "spat_shadow": 0.5,
                "env_attack": 0.02,
                "env_decay": 0.1,
                "env_sustain": 0.3,
                "env_release": 0.2,
                "noise_color": 1,  # Pink
                "color_amt": 0.8,
                "vhs_wow": 0.4,
                "vhs_flutter": 0.3,
                "glitch_mix": 0.3,
                "srcrush_amt": 0.1,
                "bitcrush_amt": 0.2,
                "thermo": True,
                "lattice": False,
                "burst": True,
                "temperature": 0.5,
                "lat_rate": 200.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.05,
                "seed": 19780815,
            }
        )
        
        self.presets["Deep Space"] = Preset(
            name="Deep Space",
            description="Sparse cosmic drones",
            author="aureonoise",
            params={
                "rate": 2.0,
                "baselen_ms": 800.0,
                "len_phi": 0.9,
                "width": 2.0,
                "itd_us": 800.0,
                "ild_db": 10.0,
                "hemis_coupling": 0.9,
                "spat_min_deg": 20.0,
                "spat_min_ms": 50.0,
                "spat_ipd": 0.7,
                "spat_shadow": 0.9,
                "env_attack": 0.4,
                "env_decay": 0.5,
                "env_sustain": 0.6,
                "env_release": 0.5,
                "noise_color": 2,  # Brown
                "color_amt": 0.9,
                "vhs_wow": 0.2,
                "vhs_flutter": 0.1,
                "glitch_mix": 0.2,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                "thermo": True,
                "lattice": True,
                "burst": False,
                "temperature": 0.2,
                "lat_rate": 50.0,
                "lat_eps": INV_PHI_CU * 0.5,
                "lat_gamma": PHI * 0.8,
                "lat_sigma": 0.02,
                "seed": 20010101,
            }
        )

        # ── Therapeutic presets ──────────────────────────────────────────

        self.presets["EMDR Bilateral"] = Preset(
            name="EMDR Bilateral",
            description="L-R alternation for EMDR reprocessing therapy",
            author="aureonoise",
            params={
                # Timing
                "rate": 8.0,
                "baselen_ms": 120.0,
                "len_phi": 0.7,
                # Spatial
                "width": 1.0,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.90,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.6,
                "spat_shadow": 0.7,
                # Envelope
                "env_attack": 0.18,
                "env_decay": 0.28,
                "env_sustain": 0.55,
                "env_release": 0.30,
                # Timbre
                "noise_color": 1,  # Pink
                "color_amt": 0.7,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.0,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                # Stochastic
                "thermo": True,
                "lattice": True,
                "burst": False,
                "temperature": 0.35,
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "burst_floor": 0.35,
                "burst_phi_mix": 0.6,
                # Externalization
                "externalization": 0.5,
                # Dialogue
                "dialogue_on": True,
                "dialogue_strength": 0.8,
                "dialogue_memory": 0.6,
                "dialogue_phi_mix": 0.8,
                # Phi-Pan + Bilateral
                "phi_pan": False,
                "bilateral_on": True,
                "bilateral_rate": 1.5,
                "bilateral_amount": 0.85,
                # Noise modes
                "noise_mode": 1,  # Pink
                "aureo_decay": 0.3,
                "aureo_stride": 1.0,
                "aureo_harmonics": 12,
                "quantum_detail": 0.7,
                "quantum_base": 220.0,
                "velvet_density": 2000.0,
                # Modal
                "modal_on": False,
                "modal_mix": 0.3,
                "modal_decay": 0.5,
                "modal_preset": 0,
                "modal_mirror": 0.3,
                "modal_feedback": 0.1,
                # Phi model
                "phi_distance": 1.5,
                "phi_elev": 0.0,
                # System
                "seed": 20251010,
            }
        )

        self.presets["ASMR Intimate"] = Preset(
            name="ASMR Intimate",
            description="Micro-transients with proximity sensation for ASMR",
            author="aureonoise",
            params={
                # Timing
                "rate": 6.0,
                "baselen_ms": 200.0,
                "len_phi": 0.5,
                # Spatial
                "width": 0.6,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.6,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.70,
                "spat_shadow": 0.7,
                # Envelope
                "env_attack": 0.25,
                "env_decay": 0.35,
                "env_sustain": 0.4,
                "env_release": 0.4,
                # Timbre
                "noise_color": 2,  # Brown
                "color_amt": 0.8,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.0,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                # Stochastic
                "thermo": True,
                "lattice": True,
                "burst": True,
                "temperature": 0.25,
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "burst_floor": 0.2,
                "burst_phi_mix": 0.7,
                # Externalization
                "externalization": 0.1,
                # Dialogue
                "dialogue_on": False,
                "dialogue_strength": 0.6,
                "dialogue_memory": 0.5,
                "dialogue_phi_mix": 0.75,
                # Phi-Pan + Bilateral
                "phi_pan": False,
                "bilateral_on": False,
                "bilateral_rate": 1.0,
                "bilateral_amount": 0.8,
                # Noise modes
                "noise_mode": 2,  # Brown
                "aureo_decay": 0.3,
                "aureo_stride": 1.0,
                "aureo_harmonics": 12,
                "quantum_detail": 0.7,
                "quantum_base": 220.0,
                "velvet_density": 2000.0,
                # Modal
                "modal_on": False,
                "modal_mix": 0.3,
                "modal_decay": 0.5,
                "modal_preset": 0,
                "modal_mirror": 0.3,
                "modal_feedback": 0.1,
                # Phi model
                "phi_distance": 1.5,
                "phi_elev": 0.0,
                # System
                "seed": 20251010,
            }
        )

        self.presets["Sleep Pink"] = Preset(
            name="Sleep Pink",
            description="Slow-wave sleep promotion with pink noise",
            author="aureonoise",
            params={
                # Timing
                "rate": 3.0,
                "baselen_ms": 500.0,
                "len_phi": 0.4,
                # Spatial
                "width": 0.8,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.6,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.6,
                "spat_shadow": 0.7,
                # Envelope
                "env_attack": 0.3,
                "env_decay": 0.3,
                "env_sustain": 0.7,
                "env_release": 0.5,
                # Timbre
                "noise_color": 1,  # Pink
                "color_amt": 0.8,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.0,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                # Stochastic
                "thermo": True,
                "lattice": False,
                "burst": False,
                "temperature": 0.15,
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "burst_floor": 0.35,
                "burst_phi_mix": 0.6,
                # Externalization
                "externalization": 0.0,
                # Dialogue
                "dialogue_on": False,
                "dialogue_strength": 0.6,
                "dialogue_memory": 0.5,
                "dialogue_phi_mix": 0.75,
                # Phi-Pan + Bilateral
                "phi_pan": False,
                "bilateral_on": False,
                "bilateral_rate": 1.0,
                "bilateral_amount": 0.8,
                # Noise modes
                "noise_mode": 1,  # Pink
                "aureo_decay": 0.3,
                "aureo_stride": 1.0,
                "aureo_harmonics": 12,
                "quantum_detail": 0.7,
                "quantum_base": 220.0,
                "velvet_density": 2000.0,
                # Modal
                "modal_on": False,
                "modal_mix": 0.3,
                "modal_decay": 0.5,
                "modal_preset": 0,
                "modal_mirror": 0.3,
                "modal_feedback": 0.1,
                # Phi model
                "phi_distance": 1.5,
                "phi_elev": 0.0,
                # System
                "seed": 20251010,
            }
        )

        self.presets["Focus Brown"] = Preset(
            name="Focus Brown",
            description="Attention and working memory enhancement with brown noise",
            author="aureonoise",
            params={
                # Timing
                "rate": 10.0,
                "baselen_ms": 100.0,
                "len_phi": 0.6,
                # Spatial
                "width": 0.9,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.6,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.6,
                "spat_shadow": 0.7,
                # Envelope
                "env_attack": 0.18,
                "env_decay": 0.28,
                "env_sustain": 0.55,
                "env_release": 0.30,
                # Timbre
                "noise_color": 2,  # Brown
                "color_amt": 0.75,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.0,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                # Stochastic
                "thermo": True,
                "lattice": True,
                "burst": False,
                "temperature": 0.30,
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "burst_floor": 0.35,
                "burst_phi_mix": 0.6,
                # Externalization
                "externalization": 0.3,
                # Dialogue
                "dialogue_on": False,
                "dialogue_strength": 0.6,
                "dialogue_memory": 0.5,
                "dialogue_phi_mix": 0.75,
                # Phi-Pan + Bilateral
                "phi_pan": False,
                "bilateral_on": False,
                "bilateral_rate": 1.0,
                "bilateral_amount": 0.8,
                # Noise modes
                "noise_mode": 2,  # Brown
                "aureo_decay": 0.3,
                "aureo_stride": 1.0,
                "aureo_harmonics": 12,
                "quantum_detail": 0.7,
                "quantum_base": 220.0,
                "velvet_density": 2000.0,
                # Modal
                "modal_on": False,
                "modal_mix": 0.3,
                "modal_decay": 0.5,
                "modal_preset": 0,
                "modal_mirror": 0.3,
                "modal_feedback": 0.1,
                # Phi model
                "phi_distance": 1.5,
                "phi_elev": 0.0,
                # System
                "seed": 20251010,
            }
        )

        self.presets["Theta Drift"] = Preset(
            name="Theta Drift",
            description="Bridge EMDR+ASMR with theta-band oscillation",
            author="aureonoise",
            params={
                # Timing
                "rate": 4.0,
                "baselen_ms": 250.0,
                "len_phi": 0.8,
                # Spatial
                "width": 1.0,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.6,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.6,
                "spat_shadow": 0.7,
                # Envelope
                "env_attack": 0.18,
                "env_decay": 0.28,
                "env_sustain": 0.55,
                "env_release": 0.30,
                # Timbre
                "noise_color": 1,  # Pink
                "color_amt": 0.7,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.0,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                # Stochastic
                "thermo": True,
                "lattice": True,
                "burst": True,
                "temperature": 0.35,
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "burst_floor": 0.3,
                "burst_phi_mix": 0.65,
                # Externalization
                "externalization": 0.2,
                # Dialogue
                "dialogue_on": True,
                "dialogue_strength": 0.5,
                "dialogue_memory": 0.4,
                "dialogue_phi_mix": 0.75,
                # Phi-Pan + Bilateral
                "phi_pan": True,
                "bilateral_on": False,
                "bilateral_rate": 1.0,
                "bilateral_amount": 0.8,
                # Noise modes
                "noise_mode": 1,  # Pink
                "aureo_decay": 0.3,
                "aureo_stride": 1.0,
                "aureo_harmonics": 12,
                "quantum_detail": 0.7,
                "quantum_base": 220.0,
                "velvet_density": 2000.0,
                # Modal
                "modal_on": False,
                "modal_mix": 0.3,
                "modal_decay": 0.5,
                "modal_preset": 0,
                "modal_mirror": 0.3,
                "modal_feedback": 0.1,
                # Phi model
                "phi_distance": 1.5,
                "phi_elev": 0.0,
                # System
                "seed": 20251010,
            }
        )

        self.presets["Hemispheric Bridge"] = Preset(
            name="Hemispheric Bridge",
            description="Alpha-band bilateral stimulation for corpus callosum synchronization",
            author="aureonoise",
            params={
                # Timing
                "rate": 8.0,
                "baselen_ms": 100.0,
                "len_phi": 0.7,
                # Spatial
                "width": 1.0,
                "itd_us": 600.0,
                "ild_db": 6.0,
                "hemis_coupling": 0.85,
                "spat_min_deg": 12.0,
                "spat_min_ms": 35.0,
                "spat_ipd": 0.6,
                "spat_shadow": 0.7,
                # Envelope
                "env_attack": 0.15,
                "env_decay": 0.25,
                "env_sustain": 0.6,
                "env_release": 0.25,
                # Timbre — Pink (1/f = brain's spectral structure)
                "noise_color": 1,
                "color_amt": 0.7,
                "vhs_wow": 0.35,
                "vhs_flutter": 0.25,
                "glitch_mix": 0.0,
                "srcrush_amt": 0.0,
                "bitcrush_amt": 0.0,
                # Stochastic
                "thermo": True,
                "lattice": True,
                "burst": False,
                "temperature": 0.35,  # peak stochastic resonance
                "lat_rate": 250.0,
                "lat_eps": INV_PHI_CU,
                "lat_gamma": PHI,
                "lat_sigma": 0.06,
                "burst_floor": 0.35,
                "burst_phi_mix": 0.6,
                # Externalization
                "externalization": 0.35,
                # Dialogue (interhemispheric coherence)
                "dialogue_on": True,
                "dialogue_strength": 0.7,
                "dialogue_memory": 0.6,
                "dialogue_phi_mix": 0.8,
                # Phi-Pan + Bilateral
                "phi_pan": False,
                "bilateral_on": True,
                "bilateral_rate": 1.0,
                "bilateral_amount": 0.75,
                # Noise modes — Pink (1/f)
                "noise_mode": 1,
                "aureo_decay": 0.3,
                "aureo_stride": 1.0,
                "aureo_harmonics": 12,
                "quantum_detail": 0.7,
                "quantum_base": 220.0,
                "velvet_density": 2000.0,
                # Modal
                "modal_on": False,
                "modal_mix": 0.3,
                "modal_decay": 0.5,
                "modal_preset": 0,
                "modal_mirror": 0.3,
                "modal_feedback": 0.1,
                # Phi model
                "phi_distance": 1.5,
                "phi_elev": 0.0,
                # System
                "seed": 20251010,
            }
        )
    
        # ── New therapeutic presets (Sprint 4/5) ──────────────────────

        _THERAPEUTIC_BASE = {
            "vhs_wow": 0.0, "vhs_flutter": 0.0, "glitch_mix": 0.0,
            "srcrush_amt": 0.0, "bitcrush_amt": 0.0,
            "spat_min_deg": 12.0, "spat_min_ms": 35.0,
            "spat_ipd": 0.6, "spat_shadow": 0.7,
            "aureo_decay": 0.3, "aureo_stride": 1.0, "aureo_harmonics": 12,
            "quantum_detail": 0.7, "quantum_base": 220.0, "velvet_density": 2000.0,
            "modal_on": False, "modal_mix": 0.3, "modal_decay": 0.5,
            "modal_preset": 0, "modal_mirror": 0.3, "modal_feedback": 0.1,
            "phi_distance": 1.5, "phi_elev": 0.0, "seed": 20251010,
        }

        def _therapeutic(name, desc, overrides):
            p = dict(_THERAPEUTIC_BASE)
            p.update(overrides)
            self.presets[name] = Preset(
                name=name, description=desc, author="aureonoise", params=p,
            )

        _therapeutic("Sleep Delta Binaural",
            "Delta binaural beat (2.5 Hz) in pink noise for deep sleep (Jirakittayakorn 2017)",
            {
                "rate": 3.0, "baselen_ms": 500.0, "len_phi": 0.4,
                "width": 0.8, "itd_us": 600.0, "ild_db": 4.0,
                "hemis_coupling": 0.5, "noise_color": 1, "color_amt": 0.8,
                "noise_mode": 1, "noise_slope": -1.0,
                "env_attack": 0.3, "env_decay": 0.3, "env_sustain": 0.7, "env_release": 0.5,
                "thermo": True, "lattice": False, "burst": False,
                "temperature": 0.10,
                "externalization": 0.0,
                "dialogue_on": False, "dialogue_strength": 0.5,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": False,
                "bilateral_rate": 1.0, "bilateral_amount": 0.8,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 30.0,
                "binaural_on": True, "binaural_carrier_hz": 200.0,
                "binaural_beat_hz": 2.5, "binaural_level": 0.08,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 10.0, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("Theta Meditation",
            "Theta binaural beat (6 Hz) in brown noise for meditation (Lavallee 2011)",
            {
                "rate": 4.0, "baselen_ms": 300.0, "len_phi": 0.6,
                "width": 0.9, "itd_us": 600.0, "ild_db": 5.0,
                "hemis_coupling": 0.6, "noise_color": 2, "color_amt": 0.75,
                "noise_mode": 2, "noise_slope": -2.0,
                "env_attack": 0.25, "env_decay": 0.3, "env_sustain": 0.6, "env_release": 0.4,
                "thermo": True, "lattice": True, "burst": False,
                "temperature": 0.18,
                "externalization": 0.15,
                "dialogue_on": True, "dialogue_strength": 0.4,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": False,
                "bilateral_rate": 1.0, "bilateral_amount": 0.8,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 20.0,
                "binaural_on": True, "binaural_carrier_hz": 250.0,
                "binaural_beat_hz": 6.0, "binaural_level": 0.06,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 10.0, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("Alpha Relax",
            "Alpha binaural beat (10 Hz) in pink noise for relaxation (Solca 2016)",
            {
                "rate": 6.0, "baselen_ms": 200.0, "len_phi": 0.6,
                "width": 1.0, "itd_us": 600.0, "ild_db": 5.0,
                "hemis_coupling": 0.6, "noise_color": 1, "color_amt": 0.7,
                "noise_mode": 1, "noise_slope": -1.0,
                "env_attack": 0.2, "env_decay": 0.28, "env_sustain": 0.55, "env_release": 0.35,
                "thermo": True, "lattice": True, "burst": False,
                "temperature": 0.20,
                "externalization": 0.2,
                "dialogue_on": True, "dialogue_strength": 0.5,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": False,
                "bilateral_rate": 1.0, "bilateral_amount": 0.8,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 15.0,
                "binaural_on": True, "binaural_carrier_hz": 300.0,
                "binaural_beat_hz": 10.0, "binaural_level": 0.07,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 10.0, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("Gamma Focus",
            "40 Hz isochronic entrainment in brown noise (MIT GENUS: 69% reduced atrophy)",
            {
                "rate": 10.0, "baselen_ms": 100.0, "len_phi": 0.6,
                "width": 0.9, "itd_us": 600.0, "ild_db": 6.0,
                "hemis_coupling": 0.6, "noise_color": 2, "color_amt": 0.75,
                "noise_mode": 2, "noise_slope": -2.0,
                "env_attack": 0.15, "env_decay": 0.25, "env_sustain": 0.6, "env_release": 0.25,
                "thermo": True, "lattice": True, "burst": False,
                "temperature": 0.22,
                "externalization": 0.3,
                "dialogue_on": False, "dialogue_strength": 0.6,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": False,
                "bilateral_rate": 1.0, "bilateral_amount": 0.8,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 0.0,
                "binaural_on": False, "binaural_carrier_hz": 250.0,
                "binaural_beat_hz": 6.0, "binaural_level": 0.0,
                "isochronic_on": True, "isochronic_carrier_hz": 400.0,
                "isochronic_rate_hz": 40.0, "isochronic_duty": 0.5, "isochronic_level": 0.12,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("Tinnitus Relief",
            "Notch-filtered pink noise at customizable tinnitus frequency",
            {
                "rate": 5.0, "baselen_ms": 300.0, "len_phi": 0.5,
                "width": 0.8, "itd_us": 500.0, "ild_db": 4.0,
                "hemis_coupling": 0.5, "noise_color": 1, "color_amt": 0.75,
                "noise_mode": 1, "noise_slope": -1.0,
                "env_attack": 0.25, "env_decay": 0.3, "env_sustain": 0.6, "env_release": 0.4,
                "thermo": True, "lattice": False, "burst": False,
                "temperature": 0.15,
                "externalization": 0.0,
                "dialogue_on": False, "dialogue_strength": 0.5,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": False,
                "bilateral_rate": 1.0, "bilateral_amount": 0.8,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 0.0,
                "binaural_on": False, "binaural_carrier_hz": 250.0,
                "binaural_beat_hz": 6.0, "binaural_level": 0.0,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 10.0, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 4000.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("CC Gentle",
            "Gentle corpus callosum stimulation (bilateral 0.7 Hz, dialogue 0.5)",
            {
                "rate": 8.0, "baselen_ms": 120.0, "len_phi": 0.7,
                "width": 1.0, "itd_us": 600.0, "ild_db": 5.0,
                "hemis_coupling": 0.7, "noise_color": 1, "color_amt": 0.7,
                "noise_mode": 1, "noise_slope": -1.0,
                "env_attack": 0.18, "env_decay": 0.28, "env_sustain": 0.55, "env_release": 0.30,
                "thermo": True, "lattice": True, "burst": False,
                "temperature": 0.20,
                "externalization": 0.3,
                "dialogue_on": True, "dialogue_strength": 0.5,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": True,
                "bilateral_rate": 0.7, "bilateral_amount": 0.6,
                "bilateral_nesting": True, "feedback_on": True, "temp_ramp_sec": 30.0,
                "binaural_on": False, "binaural_carrier_hz": 250.0,
                "binaural_beat_hz": 6.0, "binaural_level": 0.0,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 10.0, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("CC Maximum",
            "Maximum corpus callosum drive (bilateral 1.0 Hz, dialogue 0.85, feedback loop)",
            {
                "rate": 8.0, "baselen_ms": 100.0, "len_phi": 0.7,
                "width": 1.0, "itd_us": 600.0, "ild_db": 6.0,
                "hemis_coupling": 0.85, "noise_color": 1, "color_amt": 0.7,
                "noise_mode": 1, "noise_slope": -1.0,
                "env_attack": 0.15, "env_decay": 0.25, "env_sustain": 0.6, "env_release": 0.25,
                "thermo": True, "lattice": True, "burst": False,
                "temperature": 0.22,
                "externalization": 0.35,
                "dialogue_on": True, "dialogue_strength": 0.85,
                "dialogue_memory": 0.65, "dialogue_phi_mix": 0.85,
                "phi_pan": False, "bilateral_on": True,
                "bilateral_rate": 1.0, "bilateral_amount": 0.9,
                "bilateral_nesting": True, "feedback_on": True, "temp_ramp_sec": 20.0,
                "binaural_on": False, "binaural_carrier_hz": 250.0,
                "binaural_beat_hz": 6.0, "binaural_level": 0.0,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 10.0, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("Delta Reset 3Hz",
            "3 Hz delta isochronic for prayer state / deep reset (Slezin 2003)",
            {
                "rate": 3.0, "baselen_ms": 400.0, "len_phi": 0.5,
                "width": 0.8, "itd_us": 500.0, "ild_db": 4.0,
                "hemis_coupling": 0.5, "noise_color": 2, "color_amt": 0.8,
                "noise_mode": 2, "noise_slope": -1.5,
                "env_attack": 0.3, "env_decay": 0.35, "env_sustain": 0.6, "env_release": 0.45,
                "thermo": True, "lattice": False, "burst": False,
                "temperature": 0.12,
                "externalization": 0.0,
                "dialogue_on": False, "dialogue_strength": 0.5,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": False,
                "bilateral_rate": 1.0, "bilateral_amount": 0.8,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 60.0,
                "binaural_on": False, "binaural_carrier_hz": 200.0,
                "binaural_beat_hz": 3.0, "binaural_level": 0.0,
                "isochronic_on": True, "isochronic_carrier_hz": 150.0,
                "isochronic_rate_hz": 3.0, "isochronic_duty": 0.4, "isochronic_level": 0.10,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

        _therapeutic("Schumann 7.83Hz",
            "7.83 Hz Schumann resonance at theta/alpha border (Earth's electromagnetic pulse)",
            {
                "rate": 6.0, "baselen_ms": 200.0, "len_phi": 0.6,
                "width": 1.0, "itd_us": 600.0, "ild_db": 5.0,
                "hemis_coupling": 0.6, "noise_color": 1, "color_amt": 0.7,
                "noise_mode": 1, "noise_slope": -1.0,
                "env_attack": 0.2, "env_decay": 0.28, "env_sustain": 0.55, "env_release": 0.35,
                "thermo": True, "lattice": True, "burst": False,
                "temperature": 0.20,
                "externalization": 0.15,
                "dialogue_on": True, "dialogue_strength": 0.5,
                "dialogue_memory": 0.5, "dialogue_phi_mix": 0.75,
                "phi_pan": False, "bilateral_on": True,
                "bilateral_rate": 0.5, "bilateral_amount": 0.5,
                "bilateral_nesting": False, "feedback_on": False, "temp_ramp_sec": 15.0,
                "binaural_on": True, "binaural_carrier_hz": 220.0,
                "binaural_beat_hz": 7.83, "binaural_level": 0.06,
                "isochronic_on": False, "isochronic_carrier_hz": 165.0,
                "isochronic_rate_hz": 7.83, "isochronic_duty": 0.5, "isochronic_level": 0.0,
                "tinnitus_notch_hz": 0.0, "tinnitus_notch_q": 6.0,
                "burst_floor": 0.35, "burst_phi_mix": 0.6,
                "lat_rate": 250.0, "lat_eps": INV_PHI_CU, "lat_gamma": PHI, "lat_sigma": 0.06,
            })

    def list(self) -> List[str]:
        """List all preset names."""
        return list(self.presets.keys())
    
    def get(self, name: str) -> Optional[Preset]:
        """Get preset by name."""
        return self.presets.get(name)
    
    def add(self, preset: Preset):
        """Add a preset."""
        self.presets[preset.name] = preset
    
    def remove(self, name: str):
        """Remove a preset."""
        if name in self.presets:
            del self.presets[name]
    
    def save_all(self):
        """Save all presets to directory."""
        self.directory.mkdir(parents=True, exist_ok=True)
        for name, preset in self.presets.items():
            safe_name = "".join(c if c.isalnum() or c in "._- " else "_" for c in name)
            preset.save(self.directory / f"{safe_name}.json")
    
    def load_all(self):
        """Load all presets from directory."""
        if not self.directory.exists():
            return
        for path in self.directory.glob("*.json"):
            try:
                preset = Preset.load(path)
                self.presets[preset.name] = preset
            except Exception as e:
                print(f"Failed to load preset {path}: {e}")


# ── FACTORY: therapeutic preset names for quick GUI access ─────────
# Maps display name -> (preset key in PresetBank, button color RGBA)
FACTORY = {
    # Original therapeutic
    "EMDR Bilateral":    ("EMDR Bilateral",    (70, 130, 210, 255)),   # blue
    "ASMR Intimate":     ("ASMR Intimate",     (210, 140, 70, 255)),   # warm/orange
    "Sleep Pink":        ("Sleep Pink",         (80, 70, 160, 255)),    # dark blue/purple
    "Focus Brown":       ("Focus Brown",        (70, 170, 100, 255)),   # green
    "Theta Drift":       ("Theta Drift",        (70, 180, 170, 255)),   # teal
    "Hemispheric Bridge":("Hemispheric Bridge", (200, 170, 70, 255)),   # gold
    # New binaural/isochronic/notch
    "Sleep Delta Binaural": ("Sleep Delta Binaural", (60, 50, 140, 255)),  # deep purple
    "Theta Meditation":  ("Theta Meditation",   (100, 80, 180, 255)),   # purple
    "Alpha Relax":       ("Alpha Relax",        (90, 170, 220, 255)),   # sky blue
    "Gamma Focus":       ("Gamma Focus",        (220, 180, 50, 255)),   # amber
    "Tinnitus Relief":   ("Tinnitus Relief",    (180, 100, 100, 255)),  # muted red
    # Corpus callosum
    "CC Gentle":         ("CC Gentle",          (140, 190, 80, 255)),   # lime
    "CC Maximum":        ("CC Maximum",         (220, 120, 50, 255)),   # orange-red
    # Entrainment
    "Delta Reset 3Hz":   ("Delta Reset 3Hz",    (70, 60, 120, 255)),    # dark indigo
    "Schumann 7.83Hz":   ("Schumann 7.83Hz",    (100, 160, 130, 255)), # sage
}
