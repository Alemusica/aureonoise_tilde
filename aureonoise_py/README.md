# aureonoise

**φ-based granular noise/glitch texture generator**

A Python application with a high-performance Rust DSP core for generating procedural noise and glitch textures using golden ratio (φ) mathematics.

## Features

- **φ-Based Timing**: Grain scheduling and durations follow golden ratio proportions
- **Binaural Spatialization**: ITD, ILD, IPD, and head shadow modeling
- **Stochastic Processes**: Ornstein-Uhlenbeck, 3D Lattice dynamics, Hawkes burst clustering
- **Colored Noise**: White, pink, brown noise with φ-filtered states
- **Glitch Effects**: VHS drop, stutter, sample rate crush, bit crush
- **Real-time GUI**: DearPyGui interface with live parameter control
- **Preset System**: Save/load parameter configurations

## Installation

### From Source (requires Rust toolchain)

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Clone and build
git clone https://github.com/yourusername/aureonoise.git
cd aureonoise/aureonoise_py

# Install with maturin
pip install maturin
maturin develop --release

# Install Python dependencies
pip install sounddevice dearpygui numpy
```

### From Wheel (when available)

```bash
pip install aureonoise[gui]
```

## Usage

### GUI Application

```bash
aureonoise
```

Or from Python:

```python
from aureonoise.app import main
main()
```

### Python API

```python
from aureonoise import Engine, Params, PHI, INV_PHI

# Create engine
engine = Engine(sample_rate=44100.0)

# Set parameters
params = Params()
params.rate = 10.0
params.baselen_ms = 100.0
params.noise_color = 1  # Pink
engine.set_params(params)

# Render audio
import numpy as np
left, right = engine.process(44100)  # 1 second

# Save to file
from aureonoise.audio import save_wav
save_wav("output.wav", np.array(left), np.array(right))
```

### Real-time Audio

```python
from aureonoise.audio import AudioEngine

audio = AudioEngine(sample_rate=44100.0, block_size=512)
audio.set_param("rate", 12.0)
audio.set_param("glitch_mix", 0.7)
audio.start()

# ...

audio.stop()
```

### Presets

```python
from aureonoise.presets import PresetBank

bank = PresetBank()
print(bank.list())  # ['Default', 'Calm Rain', 'Digital Storm', ...]

preset = bank.get("Digital Storm")
params = preset.to_params()
engine.set_params(params)
```

## Parameters

### Timing
| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| `rate` | 0.1-60 Hz | 8.0 | Grain event rate |
| `baselen_ms` | 10-2000 ms | 120.0 | Base grain length |
| `len_phi` | 0-1 | 0.8 | φ-spread for duration |

### Spatial
| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| `width` | 0-2 | 1.0 | Stereo width |
| `itd_us` | 0-800 μs | 600.0 | Interaural time difference |
| `ild_db` | 0-12 dB | 6.0 | Interaural level difference |
| `hemis_coupling` | 0-1 | 0.6 | Hemisphere anti-correlation |
| `spat_ipd` | 0-1 | 0.6 | Interaural phase difference |
| `spat_shadow` | 0-1 | 0.7 | Head shadow filtering |

### Envelope
| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| `env_attack` | 0.01-1 | 0.18 | Attack time (normalized) |
| `env_decay` | 0.01-1 | 0.28 | Decay time |
| `env_sustain` | 0-1 | 0.55 | Sustain level |
| `env_release` | 0.01-1 | 0.30 | Release time |

### Timbre
| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| `noise_color` | 0,1,2 | 1 | White, Pink, Brown |
| `color_amt` | 0-1 | 0.65 | Color filter amount |
| `vhs_wow` | 0-1 | 0.35 | VHS wow modulation |
| `vhs_flutter` | 0-1 | 0.25 | VHS flutter modulation |
| `glitch_mix` | 0-1 | 0.5 | Glitch effect probability |
| `srcrush_amt` | 0-1 | 0.2 | Sample rate reduction |
| `bitcrush_amt` | 0-1 | 0.15 | Bit depth reduction |

### Stochastic
| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| `thermo` | bool | true | Enable Ornstein-Uhlenbeck |
| `lattice` | bool | true | Enable 3D lattice dynamics |
| `burst` | bool | true | Enable Hawkes clustering |
| `temperature` | 0-1 | 0.45 | OU process intensity |
| `lat_rate` | 1-2000 Hz | 250.0 | Lattice update rate |
| `lat_eps` | 0.01-0.5 | 1/φ³ | Lattice step size |
| `lat_gamma` | 0.5-3 | φ | Lattice coupling strength |
| `lat_sigma` | 0.01-0.3 | 0.06 | Lattice noise level |

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        Python Layer                              │
│  ┌─────────┐  ┌─────────────┐  ┌─────────────┐  ┌────────────┐ │
│  │   GUI   │  │ AudioEngine │  │   Presets   │  │    I/O     │ │
│  │ DearPy  │  │ sounddevice │  │   JSON      │  │   WAV      │ │
│  └────┬────┘  └──────┬──────┘  └──────┬──────┘  └─────┬──────┘ │
│       │              │                │               │         │
│       └──────────────┴────────────────┴───────────────┘         │
│                              │                                   │
│                        ┌─────▼─────┐                            │
│                        │  Engine   │   ← PyO3 binding           │
└────────────────────────┴─────┬─────┴────────────────────────────┘
                               │
┌──────────────────────────────▼──────────────────────────────────┐
│                        Rust DSP Core                             │
│  ┌───────┐  ┌───────┐  ┌───────┐  ┌───────┐  ┌───────────────┐ │
│  │  Rng  │  │ Weyl  │  │ Noise │  │ Ring  │  │  Stochastic   │ │
│  │xorsh64│  │ φ,φ²  │  │ W/P/B │  │Buffer │  │ OU/Lat/Hawkes │ │
│  └───────┘  └───────┘  └───────┘  └───────┘  └───────────────┘ │
│                              │                                   │
│  ┌───────────────────────────▼───────────────────────────────┐  │
│  │                      GrainPool                             │  │
│  │   Grain × 32: env, pan, itd, ild, ipd, crush, glitch     │  │
│  └───────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

## φ Constants

The engine uses golden ratio constants throughout:

| Constant | Value | Usage |
|----------|-------|-------|
| φ | 1.6180339887... | Scaling, timing |
| 1/φ | 0.6180339887... | Probability, decay |
| 1/φ² | 0.3819660113... | Amplitude, spacing |
| 1/φ³ | 0.2360679775... | Fine modulation |
| √2 | 1.4142135624... | Weyl sequence alpha |
| ρ | 1.3247179572... | Plastic ratio |

## License

MIT License - see LICENSE file.

## Credits

Original Max/MSP external: [aureonoise~](../README.md)

Based on concepts from:
- Granular synthesis (Roads, Gabor)  
- Binaural audio (Blauert)
- Stochastic processes (Uhlenbeck, Ornstein, Hawkes)
- φ mathematics (Fibonacci, Penrose)
