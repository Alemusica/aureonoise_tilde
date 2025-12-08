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
            "rate": params.rate,
            "baselen_ms": params.baselen_ms,
            "len_phi": params.len_phi,
            "width": params.width,
            "itd_us": params.itd_us,
            "ild_db": params.ild_db,
            "hemis_coupling": params.hemis_coupling,
            "spat_min_deg": params.spat_min_deg,
            "spat_min_ms": params.spat_min_ms,
            "spat_ipd": params.spat_ipd,
            "spat_shadow": params.spat_shadow,
            "env_attack": params.env_attack,
            "env_decay": params.env_decay,
            "env_sustain": params.env_sustain,
            "env_release": params.env_release,
            "noise_color": params.noise_color,
            "color_amt": params.color_amt,
            "vhs_wow": params.vhs_wow,
            "vhs_flutter": params.vhs_flutter,
            "glitch_mix": params.glitch_mix,
            "srcrush_amt": params.srcrush_amt,
            "bitcrush_amt": params.bitcrush_amt,
            "thermo": params.thermo,
            "lattice": params.lattice,
            "burst": params.burst,
            "temperature": params.temperature,
            "lat_rate": params.lat_rate,
            "lat_eps": params.lat_eps,
            "lat_gamma": params.lat_gamma,
            "lat_sigma": params.lat_sigma,
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
