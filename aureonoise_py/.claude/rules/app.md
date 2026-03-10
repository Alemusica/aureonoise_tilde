---
globs:
  - "python/aureonoise/app.py"
  - "python/aureonoise/audio.py"
  - "python/aureonoise/__init__.py"
---

# Application & GUI — Regole di contesto

## DearPyGui 2.x — API

DearPyGui 2.x usa `dpg.` namespace (non `dearpygui.dearpygui as dpg`):
```python
import dearpygui.dearpygui as dpg

dpg.create_context()
with dpg.window(label="Aureonoise"):
    dpg.add_button(label="Play", callback=on_play)
    dpg.add_slider_float(label="Volume", default_value=0.7, min_value=0, max_value=1)
dpg.create_viewport(title="Aureonoise", width=800, height=600)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()
```

## Thread safety

Il callback sounddevice gira su un thread audio separato. Regola assoluta:

```python
# CORRETTO — protetto da lock
def set_param(self, name, value):
    with self._lock:
        setattr(self._params, name, value)
        self.engine.set_params(self._params)

# SBAGLIATO — race condition
def set_param(self, name, value):
    setattr(self._params, name, value)
    self.engine.set_params(self._params)
```

## SR auto-detection

```python
dev_info = sd.query_devices(self.output_device, 'output')
actual_sr = dev_info['default_samplerate']

if actual_sr != self.sample_rate:
    self.sample_rate = actual_sr
    self.engine = Engine(actual_sr)  # ricreare con SR nativa
    self.engine.set_params(self._params)
```

Non rimuovere. PortAudio su macOS fallisce con SR mismatch.

## Metering — exponential decay

```python
self._peak_l = max(self._peak_l * 0.95, np.abs(left_np).max())
```

Il fattore 0.95 produce decay visuale naturale. Non usare peak istantaneo senza decay.

## Analysis callback — interval

```python
def on_analysis(self, callback, interval_blocks=86):
    # 86 blocks × 512 samples = 44032 samples ≈ 1 secondo a 44100 Hz
```

Non ridurre sotto 40 blocks (~0.5s) — l'analisi Welch PSD ha bisogno di abbastanza campioni per essere significativa.

## Mutual exclusion — preset conflicts

Binaural beats e isochronic tones sono meccanismi di entrainment diversi. Attivarli simultaneamente crea interferenza.

La GUI deve:
1. Disabilitare slider isochronic quando binaural e attivo (e viceversa)
2. I preset gestiscono questo automaticamente (ogni preset ha _FULL_DEFAULTS)
3. Se l'utente modifica manualmente, la GUI deve impedire combinazioni invalide

## Color palette — phi-based

I colori della GUI usano hue spacing basato su phi:
```python
hue_step = 360.0 / PHI  # ≈ 222.5°
```

Non usare colori arbitrari. Mantenere la palette phi-ratio per coerenza visiva.

## Volume limiter (DA IMPLEMENTARE)

Serve un ceiling indipendente dal gain utente:
```python
# Proposta
MAX_SAFE_AMPLITUDE = 0.89  # -1 dBFS
output = np.clip(output, -MAX_SAFE_AMPLITUDE, MAX_SAFE_AMPLITUDE)
```

Il soft_tanh Rust previene clipping ma non limita il volume percepito.
