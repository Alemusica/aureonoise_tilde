---
globs:
  - "src/dvf.rs"
  - "src/isochronic.rs"
  - "src/modal.rs"
  - "src/polyrhythm.rs"
  - "src/room.rs"
  - "src/tinnitus.rs"
  - "src/external.rs"
---

# Effects & Modulation — Regole di contesto

## Signal chain order

```
grain mix → modal → binaural/isochronic add → dvf → external → room → soft_tanh
```

Rispettare l'ordine. Il modal processa il mix dei grani, poi binaural/isochronic vengono aggiunti (non processati dal modal).

## Body resonance — frequenze pericolose

| Frequenza | Organo | Rischio |
|-----------|--------|---------|
| 5-8 Hz | Cavita toracica | Risonanza meccanica |
| 19 Hz | Bulbo oculare | Vibrazione retina |
| 0.5-1.5 Hz | Stomaco | Nausea |

I moduli isochronic, polyrhythm e modal possono tutti produrre energia concentrata a queste frequenze. Implementare guardie.

## Isochronic — duty cycle semantics

```
duty = 0.1 → impulsi molto corti (10% on, 90% off) — piu percussivo
duty = 0.5 → square wave — massimo contenuto armonico
duty = 0.9 → quasi continuo (90% on, 10% off) — piu tonale
```

Per entrainment gamma (40 Hz MIT GENUS): duty 0.3-0.5 e ottimale.

## Isochronic — rischio epilessia

Frequenze 8-25 Hz con isochronic possono causare auditory driving analogo al photic driving. Questo e un rischio reale per soggetti fotosensibili/acusticamente sensibili.

Mitigation: warning in GUI (responsabilita app-specialist), ma il modulo deve esporre un flag `is_seizure_risk_range()`.

## Tinnitus notch — Butterworth 4th order

```
2 × biquad cascaded
  biquad1: 2nd-order notch
  biquad2: 2nd-order notch (same coefficients)
Result: 4th-order rolloff, >20 dB notch depth
```

I coefficienti biquad si ricalcolano SOLO quando `center_freq` o `q` cambiano. Non ricalcolare ad ogni sample.

## Modal — contralateral mirror

Il modal produce un segnale "wet" contralaterale (mono). Questo segnale viene processato in lib.rs attraverso il phi head pipeline:
1. ITD delay ring
2. ILD
3. Head shadow
4. Crossfeed

NON applicare spazializzazione dentro modal.rs — la spazializzazione avviene in lib.rs.

## Room — delay times phi-ratio

Tutti i delay times del reverb sono scalati per phi:
- Early reflections: base_delay, base_delay * PHI, base_delay * PHI^2, base_delay * PHI^3
- Comb filters: idem

Non usare rapporti interi (1:2:3:4) — producono colorazione metallica.

## Biquad — regola dirty flag

```rust
if self.dirty {
    self.compute_coefficients();
    self.dirty = false;
}
// process sample with current coefficients
```

Ricalcolare coefficienti ad ogni sample e uno spreco e puo causare artefatti.

## Test

```bash
maturin develop && pytest tests/test_therapeutic.py -v -k "modal or isochronic or tinnitus or room"
```
