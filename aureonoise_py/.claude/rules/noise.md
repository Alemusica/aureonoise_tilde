---
globs:
  - "src/noise.rs"
  - "src/stoch.rs"
  - "src/burst.rs"
  - "tests/test_noise_slope.py"
  - "tests/test_burst.py"
---

# Noise & Stochastic — Regole di contesto

## Slope spettrale — definizioni fisiche

| Noise | Slope target | PSD ∝ | Tolleranza |
|-------|-------------|-------|------------|
| White | 0 dB/oct | f^0 | ±0.5 |
| Pink | -1 dB/oct | 1/f | ±1.0 |
| Brown | -2 dB/oct | 1/f^2 | ±1.5 |

Misurata con Welch PSD, linear regression su log2(freq) vs dB(PSD), range 100-8000 Hz.

La tolleranza larga per brown (±1.5) tiene conto del grain envelope windowing e output saturation a valle.

## PinkFilter — non toccare

Il Paul Kellet 6-stage IIR produce -1 dB/oct ±0.5 dB. E il gold standard per pink noise generation. Se la slope pink misurata e sbagliata, il problema e a valle (envelope, soft_tanh), non nel filtro.

```rust
// Kellet coefficients — DO NOT MODIFY
b0 += (white - b0) * 0.99886;
b1 += (white - b1) * 0.99332;
b2 += (white - b2) * 0.96900;
b3 += (white - b3) * 0.86650;
b4 += (white - b4) * 0.55000;
b5 += (white - b5) * 0.16073;
```

## SpectralTilt brown — BUG NOTO

Il path brown del SpectralTilt usa un doppio integratore leaky:
```rust
z1 = z1 * 0.995 + x * 0.005;    // primo polo
z2 = z2 * 0.985 + z1 * 0.015;   // secondo polo
out = soft_tanh(z2 * 2.4) * 0.5; // saturazione
```

Questo produce -6/-10 dB/oct invece del target -2 dB/oct. Root cause: doppia integrazione + saturazione tanh.

Fix proposto: singolo polo a ~20 Hz per -6 dB/oct, poi compensare con crossfade white per raggiungere -2 dB/oct. Oppure: IIR dedicato.

## Hawkes — rate in Hz

```rust
pub fn tick(&mut self, dt: f64) -> bool {
    self.intensity = self.base_rate + ...;
    // base_rate e in Hz (eventi/secondo)
    // dt e in secondi (1.0 / sample_rate)
}
```

Se cambi base_rate, verifica con `tests/test_burst.py` che la frequenza misurata corrisponde.

## OU — calibrazione theta

Per modulazione a frequenza F Hz: theta ≈ 2*pi*F.
- 0.1 Hz macro-modulation (heart-brain coherence): theta ≈ 0.628
- 1.0 Hz bilateral: theta ≈ 6.28

## Stochastic resonance

Non ancora implementata. Quando si implementa:
- Livello: -15 a -20 dB sotto soglia di percezione del segnale coerente
- Tipo: noise bianco (flat spectrum)
- Scopo: amplificare pattern coerenti, esporre pattern isolati (Collins 1995)

## Test

```bash
maturin develop && pytest tests/test_noise_slope.py tests/test_burst.py -v
```
