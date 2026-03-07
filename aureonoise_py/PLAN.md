# AUREONOISE — Piano di Implementazione Completo
## Sessione 2026-03-07 — Sintesi di 12+ agenti di ricerca e analisi

---

## FASE 0: BUG CRITICI (blockers — senza questi il dialogue system opera su dati corrotti)

### T0.1 — Fix BilateralOscillator phase advance [CRITICO]
- **File:** `src/dialogue.rs:483`
- **Bug:** `tick()` avanza fase di `rate/sr` (1 sample) ma chiamato 1 volta per GRANO, non per sample
- **Effetto:** A 8 grains/sec, 1 Hz bilateral impiega 92 MINUTI per ciclo
- **Fix:** `self.phase += self.rate * elapsed_samples / sr` oppure spostare tick nel loop per-sample
- **Complessità:** Bassa (5 righe)

### T0.2 — Riordinare pan chain: dialogue PRIMA di burst [CRITICO]
- **File:** `src/lib.rs:787-832` (spawn_grain)
- **Bug:** BurstEngine modifica pan PRIMA che DialogueSystem lo valuti
- **Effetto:** center_pull del burst distrugge i rapporti phi. Pan 0.8 → 0.53, ratio score crolla da 0.9 a 0.15. Handshake sistematicamente perso durante burst Hawkes.
- **Fix:** Spostare dialogue.evaluate() prima di burst_engine.apply_position(). 15 righe da riordinare.
- **Ordine nuovo:** raw → lattice → OU → hemisphere → DIALOGUE → burst → phi_pan → bilateral → ITD/ILD
- **Complessità:** Bassa (riordinare codice esistente)

### T0.3 — Fix temperature default [TRIVIALE]
- **File:** `src/lib.rs:226` (Params default)
- **Bug:** Default temperature=0.45, ottimale per stochastic resonance è 0.18-0.25
- **Effetto:** A 0.45 il noise OU sul pan è +/-0.126, troppo alto — phi-ratio drowning nel rumore
- **Fix:** Cambiare default a 0.22. Aggiornare preset terapeutici a range 0.18-0.25
- **Complessità:** Triviale (1 costante + preset values)

---

## FASE 1: CONNESSIONI MANCANTI (codice esistente mai collegato)

### T1.1 — Connettere phi_model.rs al motore [ALTO]
- **File:** `src/phi_model.rs` (742 righe), `src/lib.rs` (spawn_grain + process_block)
- **Stato:** Modello geometrico testa ellissoidale con pinna 5-tap FIR, torso 2-tap, distanza, assorbimento aria — tutto codificato, mai chiamato
- **Task:** Wire `PhiModel::compute()` nel percorso ITD/ILD di spawn_grain. Aggiungere `torso_mix: f64` a Params (default 0.0, terapeutico 0.15-0.25)
- **Complessità:** Media

### T1.2 — Consumare fib_gap_queue [ALTO]
- **File:** `src/dialogue.rs:336-358` (writer), `src/lib.rs:908-945` (schedule_gap_samples)
- **Bug:** fib_gap_queue scritta ma fib_gap_queue_pos mai incrementato, mai letto
- **Fix:** Aggiungere `pub fn pop_queued_gap(&mut self) -> Option<f64>` a DialogueSystem. In schedule_gap_samples: `if let Some(g) = self.dialogue.pop_queued_gap() { return g.round().max(1.0) as i32; }`
- **Effetto:** Dopo handshake, prossimi 2-3 grani a intervalli Fibonacci = "eco temporale"
- **Complessità:** Bassa-Media

### T1.3 — Aggiungere per-grain ring_start [MEDIO]
- **File:** `src/grain.rs`, `src/lib.rs` (spawn_grain + process_block)
- **Bug:** Tutti i grani leggono dalla stessa posizione del ring buffer → zero decorrelazione inter-grano
- **Fix:** Aggiungere `ring_start: u32` a Grain. Al spawn, `grain.ring_start = self.ring.write_pos() - offset` dove offset è quasi-random (Weyl)
- **Complessità:** Bassa-Media

---

## FASE 2: BILATERAL ENGINE v2 (sostituzione completa della sinusoide)

### T2.1 — Raised-cosine trajectory con dwell [ALTO]
- **File:** `src/dialogue.rs:442-509` (BilateralOscillator)
- **Stato attuale:** `(TWO_PI * self.phase).sin()` — pan sinusoidale continuo
- **Ricerca:** Raised-cosine con 15% dwell > sinusoide pura (EMDR+, van den Hout 2011)
- **Fix:** Struttura ciclo: dwell_L(15%) → transit_L→R(35%) → dwell_R(15%) → transit_R→L(35%). Onset < 5ms (Tukey taper alpha=0.1-0.2). Amplitude coupling: `amp_mod = 0.85 + 0.15 * abs(bilateral_pan)` (più forte alle estremità)
- **Nuovi params:** `bilateral_dwell: f64` (default 0.15), `bilateral_crossfade: f64` (default 0.7)
- **Complessità:** Media

### T2.2 — Estendere range bilateral rate [BASSO]
- **File:** `src/dialogue.rs:469` (clamp)
- **Stato:** [0.5, 2.0] Hz
- **Ricerca:** Sleep/delta necessita 0.3 Hz, theta moving sounds funzionano fino a 6 Hz
- **Fix:** Estendere a [0.3, 6.0] Hz
- **Complessità:** Triviale

### T2.3 — Theta-gamma nesting [MEDIO]
- **Ricerca:** Lisman-Jensen 2013. Ciclo theta contiene 4-8 gamma. Grain rate deve essere 4:1 a 8:1 del bilateral rate.
- **Fix:** Quando bilateral_on, grain rate opzionalmente agganciata come multiplo intero del bilateral rate. Param: `bilateral_nesting: bool` (default false)
- **Complessità:** Bassa

---

## FASE 3: HANDSHAKE ENGINE v2 (multi-target + Kaplan)

### T3.1 — Multi-target ratio scoring [MEDIO]
- **File:** `src/dialogue.rs:375` (ratio_score)
- **Bug:** Solo 1 target Fibonacci (walk +/-1). PHI^2 (2.618) non riconosciuto
- **Fix:** Array SACRED_RATIOS [0.382, 0.500, 0.618, 0.667, 1.0, 1.5, 1.618, 2.0, 2.618, 4.236]. Valutare contro tutti, prendere max score.
- **Complessità:** Bassa

### T3.2 — Ribilanciare pesi handshake [BASSO]
- **File:** `src/dialogue.rs:247`
- **Stato:** 0.55*pan + 0.15*amp + 0.15*dur + 0.15*gap (pan over-pesato)
- **Fix:** 0.40*pan + 0.12*amp + 0.12*dur + 0.12*gap + 0.24*temporal_score (dimensione Kaplan: ritmo degli handshake)
- **Ricerca:** Kaplan ISS (1997) — la consistenza temporale degli handshake è più importante della singola qualità
- **Complessità:** Bassa

### T3.3 — PLV-like metric [MEDIO]
- **File:** `src/dialogue.rs` (nuovo)
- **Ricerca:** Distinguere "high coherence because frequent" da "high coherence because RHYTHMIC"
- **Fix:** `handshake_plv = |mean(exp(j*2pi*t_k/T_bilateral))| over last 8-16 handshakes`. Esporre via Engine.
- **Complessità:** Media

### T3.4 — Hawkes-aware scoring [BASSO]
- **File:** `src/dialogue.rs:375`
- **Fix:** Durante burst (hawkes.lambda >> base), allargare tolerance: `tolerance = 0.4 + 0.3 * hawkes_weight`
- **Complessità:** Bassa

### T3.5 — Fibonacci walk con salti [BASSO]
- **File:** `src/dialogue.rs:190-210`
- **Stato:** Walk solo +/-1 nel FIB array
- **Fix:** Permettere +/-2, +/-3 con probabilità decrescente. Allarga spazio dei ratio raggiungibili.
- **Complessità:** Bassa

---

## FASE 4: CLOSED-LOOP FEEDBACK (BAC-inspired)

### T4.1 — Coherence feedback loop [ALTO]
- **File:** `src/lib.rs` (process_block, nuovi campi Engine)
- **Stato:** Coherence calcolata ma MAI letta dal codice Rust. Sistema completamente open-loop.
- **Design:**
  - `coh_norm = clamp01((coherence - 0.6) / 1.2)`
  - Slow EMA (tau ~5s): `coh_slow = 0.999 * coh_slow + 0.001 * coh_norm`
  - High (>0.7): temp verso 0.22, bilateral rate -30%, OU sigma -25%
  - Low (<0.3): temp +10%, rate +30%, amount +40%
- **Nuovi campi Engine:** `effective_temp: f64`, `coh_slow: f64`, `feedback_enabled: bool`
- **Param:** `feedback_on: bool` (default false — opt-in per preset terapeutici)
- **Attenzione:** Feedback deve essere PIÙ LENTO del smoothing della coerenza per evitare oscillazioni
- **Complessità:** Media-Alta

### T4.2 — Temperature ramp nei preset [BASSO]
- **Ricerca:** BAC approach — esplorazione poi convergenza
- **Fix:** Ogni preset terapeutico specifica `temp_start` e `temp_target` con ramp time 30-60s
- **Complessità:** Bassa

---

## FASE 5: NUOVI MODULI DSP

### T5.1 — Binaural Beat Generator [ALTO — colma lacuna maggiore]
- **File:** `src/noise.rs` o nuovo `src/binaural.rs`
- **Stato:** aureonoise NON genera binaural beats. Gap significativo per applicazione terapeutica.
- **Design:** Due oscillatori sinusoidali, L=carrier, R=carrier+beat. Mix a -20/-10 dB sotto il noise floor.
- **Nuovi params:** `binaural_on: bool`, `binaural_beat_hz: f64` (0.5-40), `binaural_carrier_hz: f64` (200-500), `binaural_level: f64` (0-1)
- **CRITICO:** Beat SEPARATO dal noise path. Beat embeddato in pink noise = ZERO effetti (confermato dalla ricerca)
- **Preset:** Sleep Delta (2.5/200), Theta Meditation (6/250), Alpha Relax (10/300), Gamma Focus (40/400)
- **Complessità:** Media

### T5.2 — Continuous Spectral Slope [MEDIO]
- **File:** `src/noise.rs`
- **Stato:** Discreto noise_color (0=White, 1=Pink, 2=Brown)
- **Fix:** `noise_slope: f64` (-2.0 to +0.5). Julius O. Smith spectral tilt filter (cascade N shelving). Backward compat: White=0.0, Pink=-1.0, Brown=-2.0.
- **Complessità:** Media

### T5.3 — Tinnitus Notch Filter [BASSO-MEDIO]
- **File:** `src/noise.rs` o `src/lib.rs`
- **Ricerca:** 10 Hz AM a freq tinnitus: 19/28 pazienti showed suppression (p<0.0001). Butterworth order 4, 1 octave bandwidth.
- **Nuovi params:** `tinnitus_notch_hz: f64` (0=off, >0=center freq), `tinnitus_notch_q: f64` (default 6.0)
- **Complessità:** Bassa

### T5.4 — Isochronic Tone Generator [BASSO]
- **File:** nuovo `src/isochronic.rs`
- **Design:** Sine a carrier (150-180 Hz), AM da Tukey-windowed pulse a target rate (1-40 Hz), duty cycle 0.3-0.7
- **Nuovi params:** `isochronic_on: bool`, `isochronic_rate_hz: f64`, `isochronic_carrier_hz: f64`, `isochronic_duty: f64`
- **Complessità:** Bassa

### T5.5 — PolyrhythmClock [ENHANCEMENT]
- **File:** nuovo `src/polyrhythm.rs`
- **Design:** Due stream impulsi a ratio p:q per emisfero. Coincidenza = handshake strutturale. Layer sopra scheduling, non sostituzione.
- **Rate:** Calm 3:2 @0.5Hz, Focus 5:3 @1.0Hz, EMDR 8:5 @1.5Hz
- **Complessità:** Media

---

## FASE 6: SPATIAL AUDIO AVANZATO

### T6.1 — Near-field DVF per ASMR [MEDIO]
- **File:** `src/phi_model.rs`
- **Ricerca:** DVF a 30cm aumenta ILD di 8-15 dB alle alte frequenze
- **Fix:** 2nd-order parametric shelving filter per orecchio, cutoff ~625 Hz, gain proporzionale a 1/r^2. Estendere phi_distance range a 0.15-4.0m.
- **Complessità:** Media

### T6.2 — Room divergence minimale [BASSO]
- **File:** `src/external.rs`
- **Fix:** Schroeder allpass reverb (2 allpass + 1 comb) con delay times phi-ratio. RT60 ~300ms, 4 early reflections. Per esternalizzazione in cuffia.
- **Complessità:** Media

### T6.3 — Coherence-driven spatial morphing [BASSO]
- **Ricerca:** Width, ITD, ILD seguono curva di coerenza del dialogue
- **Fix:** Quando coherence alta → aumenta separazione spaziale (enfasi bilaterale)
- **Dipende da:** T4.1 (feedback loop)
- **Complessità:** Bassa

---

## FASE 7: AUDIO ANALYSIS & QUALITY FEEDBACK

### T7.1 — Modulo analisi offline (Python) [MEDIO]
- **File:** nuovo `python/aureonoise/analysis.py`
- **Reuse da:** `~/AbletonScripts/analyzer/master_analyzer.py` (LUFS, stereo correlation, spectral)
- **Metriche:**
  - Spectral slope (Welch PSD, R²) — già in test_noise_slope.py
  - Stereo correlation (Pearson L/R)
  - ITD measurement (cross-correlation L/R, ±1.5ms window) — da analyze_sweep_csv.py
  - ILD per 4 bande (1-2k, 2-4k, 4-8k, 8-12k Hz)
  - LUFS momentary (ITU-R BS.1770-4)
  - Handshake rate, coherence mean, PLV (da Engine API)
  - Spectral centroid/spread
  - Onset detection (rise time measurement)
  - Bilateral symmetry score
- **Output:** JSON report + opzionalmente WAV annotato
- **Complessità:** Media

### T7.2 — Test suite per qualità bilaterale [MEDIO]
- **File:** `tests/test_bilateral_quality.py`
- **Tests:**
  - Bilateral cycle time matches params (regression per T0.1)
  - Rise time < 5ms su onset bilaterale (dopo T2.1)
  - Dwell time 15% ± 2% del ciclo
  - Handshake rate > 0 dopo 5s con hemispheric bridge preset
  - Stereo correlation inversamente proporzionale a bilateral_amount
  - Spectral slope matches noise_color per tutti i modi
  - No NaN/Inf su 60s di processing per ogni preset terapeutico
- **Complessità:** Media

### T7.3 — Real-time analysis hook in AudioEngine [BASSO]
- **File:** `python/aureonoise/audio.py`
- **Fix:** Nel callback audio, opzionalmente calcola spectral features ogni N blocchi. Esporre via meter callback esteso.
- **Complessità:** Bassa

---

## FASE 8: PRESET & GUI

### T8.1 — Nuovi preset terapeutici [BASSO]
- **File:** `python/aureonoise/presets.py`
- **Preset da aggiungere:**
  - 3 Hz Delta Reset (Slezin prayer state)
  - 7.83 Hz Schumann (theta/alpha border)
  - CC Gentle (bilateral 0.7Hz, amount 0.6, dialogue 0.5)
  - CC Maximum (bilateral 1.0Hz, amount 0.9, dialogue 0.85, raised-cosine+dwell)
  - Sleep Delta (binaural 2.5/200 + pink noise, temp 0.10)
  - Theta Meditation (binaural 6/250 + brown noise)
  - Alpha Relax (binaural 10/300 + pink noise)
  - Gamma Focus (isochronic 40/400 + brown noise)
  - Tinnitus Notch (personalizzabile)
- **Temperature corretta:** tutti i preset terapeutici a 0.18-0.25

### T8.2 — GUI tabs per nuovi moduli [MEDIO]
- **File:** `python/aureonoise/app.py`
- **Fix:** Tab per binaural beat, isochronic, tinnitus notch, analysis display
- **Complessità:** Media

### T8.3 — Dead params cleanup [TRIVIALE]
- **File:** `src/lib.rs` (Params)
- **Fix:** Rimuovere o collegare: spat_min_deg, spat_min_ms, phi_distance, phi_elev
- **Complessità:** Triviale

---

## DIPENDENZE

```
T0.1 ──────────────────────────┐
T0.2 ──────────────────────────┤
T0.3 ──────────────────────────┤
                               ▼
T1.1 (phi_model) ─────────── T6.1 (DVF)
T1.2 (fib_gap_queue) ──┐
T1.3 (ring_start)      │
                        ▼
T2.1 (raised-cosine) ─── T7.2 (test bilateral quality)
T2.2 (rate range)
T2.3 (theta-gamma)

T3.1 (multi-ratio) ───┐
T3.2 (pesi)            │
T3.3 (PLV) ────────────┤
T3.4 (hawkes-aware)    ▼
T3.5 (fib jumps)     T4.1 (feedback loop) ─── T6.3 (spatial morphing)
                                               T4.2 (temp ramp)

T5.1 (binaural beat) ──┐
T5.2 (spectral slope)  │
T5.3 (tinnitus notch)  ├── T8.1 (nuovi preset) ── T8.2 (GUI)
T5.4 (isochronic)      │
T5.5 (polyrhythm) ─────┘

T7.1 (analysis module) ── T7.3 (real-time hook)
```

---

## ORDINE DI ESECUZIONE CONSIGLIATO

### Sprint 1: Bug critici + connessioni (1-2 giorni)
1. T0.1 — Fix bilateral phase advance
2. T0.2 — Riordinare pan chain
3. T0.3 — Temperature default
4. T1.2 — Consumare fib_gap_queue
5. T1.3 — Per-grain ring_start
6. T8.3 — Dead params cleanup

### Sprint 2: Bilateral Engine v2 (1-2 giorni)
7. T2.1 — Raised-cosine con dwell
8. T2.2 — Estendere rate range
9. T3.1 — Multi-target ratio scoring
10. T3.2 — Ribilanciare pesi handshake
11. T3.4 — Hawkes-aware scoring
12. T3.5 — Fibonacci walk con salti

### Sprint 3: Feedback + Spatial (2-3 giorni)
13. T1.1 — Connettere phi_model.rs
14. T4.1 — Coherence feedback loop
15. T4.2 — Temperature ramp preset
16. T3.3 — PLV metric
17. T2.3 — Theta-gamma nesting

### Sprint 4: Nuovi moduli DSP (2-3 giorni)
18. T5.1 — Binaural beat generator
19. T5.2 — Continuous spectral slope
20. T5.3 — Tinnitus notch
21. T5.4 — Isochronic tones

### Sprint 5: Quality + Polish (1-2 giorni)
22. T7.1 — Analysis module
23. T7.2 — Test suite bilateral quality
24. T7.3 — Real-time analysis hook
25. T8.1 — Nuovi preset terapeutici
26. T8.2 — GUI tabs

### Sprint 6: Enhancement (opzionale)
27. T5.5 — PolyrhythmClock
28. T6.1 — DVF near-field
29. T6.2 — Room divergence
30. T6.3 — Coherence-driven spatial

---

## CONTESTO MONADE (per recovery post-compaction)

Pacchetti liofilizzati:
- `aureonoise-bilateral-8agent-analysis` — bug + architettura
- `aureonoise-audio-analysis-inventory` — strumenti analisi
- `aureonoise-session-state-2026-03-07` — stato sessione
- `aureonoise-corpus-callosum-research` — ricerca CC (sommario)
- `aureonoise-stochastic-handshake-coupling` — analisi accoppiamento (sommario)
- `aureonoise-corpus-callosum-detailed-v2` — ricerca CC (dettaglio completo)
- `aureonoise-stochastic-handshake-detailed-v2` — analisi accoppiamento (dettaglio)
- `aureonoise-agent1-therapeutic-algorithms` — algoritmi terapeutici SOTA
- `aureonoise-agent2-codebase-spatial` — analisi codebase spaziale
- `aureonoise-agent3-practitioner-forum` — conoscenza practitioner
- `aureonoise-agent4-neurofeedback-entrainment` — neurofeedback
- `aureonoise-agent5-russian-soviet` — ricerca sovietica
- `aureonoise-agent6-fringe-solfeggio` — solfeggio/alternative
- `aureonoise-agent7-military-intelligence` — programmi militari/intelligence
- `aureonoise-agent8-esoteric-traditional` — tradizioni non-occidentali

Ricerca primaria: `research/SYNTHESIS.md` (251 righe, 4 tier di evidenza)

---

## FONTI CHIAVE

| Dominio | Fonte | Impatto |
|---------|-------|---------|
| Kaplan ISS | Fingelkurts 1997, 2001, 2003 | Framework del dialogue system |
| EMDR speed | van den Hout 2011 (PMC4387929) | 1.2 Hz > 0.8 Hz |
| EMDR+ | PMC10377614 | 432 Hz, 0.4 Hz, d=-6.1 |
| Theta-gamma | Lisman-Jensen 2013 (Cell/Neuron) | 4-8 gamma per theta |
| Alpha coherence | Solca 2016 | Binaural beats → alpha coherence |
| Gamma bilateral | Preisig 2021 (PNAS) | Anti-phase gamma |
| SR optimale | PMC3954722 | Pink noise, -15/-20 dB |
| CC fibre | Caminiti 2013 | 4-12ms, event-triggered |
| 40 Hz gamma | MIT GENUS | 69% reduced atrophy |
| BAC | Konstantinov/Pavlov Institute | Closed-loop neurofeedback |
| Rise time | PMC3756151 | N1 driven by onset dP/dt |
