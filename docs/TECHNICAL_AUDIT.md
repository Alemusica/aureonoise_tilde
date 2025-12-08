# aureonoise~ Technical Audit
## Signal Flow Analysis & Parameter Interaction Matrix

**Data**: 2025-12-08  
**Tipo di Audit**: System Architecture Review (SAR) + Data Flow Analysis (DFA)

---

## 📊 1. Signal Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           aureonoise~ DSP SIGNAL FLOW                           │
└─────────────────────────────────────────────────────────────────────────────────┘

                          ┌─────────────────────┐
                          │    NOISE SOURCE     │
                          │  ┌───────────────┐  │
                          │  │     RNG       │  │
                          │  │  (xorshift)   │  │
                          │  └───────┬───────┘  │
                          │          │          │
                          │  ┌───────▼───────┐  │
                          │  │ NoiseColor    │◄─┼── [color] [color_amt]
                          │  │ White/Pink/   │  │
                          │  │ Brown         │  │
                          │  └───────┬───────┘  │
                          │          │          │
                          │  ┌───────▼───────┐  │
                          │  │  soft_tanh    │  │
                          │  │   (×1.2)      │  │
                          │  └───────┬───────┘  │
                          └──────────┼──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │    RING BUFFER      │
                          │   (131072 samples)  │
                          └──────────┬──────────┘
                                     │
           ┌─────────────────────────┼─────────────────────────┐
           │                         │                         │
           ▼                         ▼                         ▼
    ┌──────────────┐          ┌──────────────┐          ┌──────────────┐
    │   GRAIN 0    │          │   GRAIN 1    │   ...    │   GRAIN 31   │
    └──────┬───────┘          └──────┬───────┘          └──────┬───────┘
           │                         │                         │
           └─────────────────────────┼─────────────────────────┘
                                     │
┌────────────────────────────────────┼────────────────────────────────────────────┐
│                           PER-GRAIN PROCESSING                                  │
│  ┌─────────────────────────────────┼─────────────────────────────────────────┐  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 1. ITD READ (Lagrange interpolation)                                │  │  │
│  │  │    ring_read_stereo_itd_frac(ring, wi, itd) → sL, sR                │  │  │
│  │  │    ◄── [itd_us] [vhs_mod] [lattice]                                 │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 2. PAN + ILD GAIN                                                   │  │  │
│  │  │    sL *= panL * gL    sR *= panR * gR                               │  │  │
│  │  │    ◄── [width] [ild_db] [hemis_coupling]                            │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 3. CROSSFEED (binaural)                                             │  │  │
│  │  │    sL = sL + crossfeed * sR                                         │  │  │
│  │  │    sR = sR + crossfeed * sL                                         │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 4. SAMPLE-RATE CRUSH (lo-fi)                                        │  │  │
│  │  │    Hold samples for sr_holdN samples                                │  │  │
│  │  │    ◄── [srcrush_amt]                                                │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 5. BIT CRUSH (quantization)                                         │  │  │
│  │  │    sL = round(sL * q_levels) / q_levels                             │  │  │
│  │  │    ◄── [bitcrush_amt]                                               │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 6. GLITCH KIND PROCESSING                                           │  │  │
│  │  │    VhsDrop:  *= 0.5 + 0.5*(1-|vhs_mod|)                             │  │  │
│  │  │    Stutter:  *= 0.2 every 8 samples                                 │  │  │
│  │  │    Aliaser:  (no extra processing)                                  │  │  │
│  │  │    ◄── [glitch_mix]                                                 │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 7. IPD (Interaural Phase Decorrelation)                             │  │  │
│  │  │    All-pass filter: sL = allpass(sL, +coeff)                        │  │  │
│  │  │                     sR = allpass(sR, -coeff)                        │  │  │
│  │  │    ◄── [spat_ipd]                                                   │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 8. HEAD SHADOW (contralateral LP filter)                            │  │  │
│  │  │    if pan > 0: sL = LP(sL)                                          │  │  │
│  │  │    if pan < 0: sR = LP(sR)                                          │  │  │
│  │  │    ◄── [spat_shadow]                                                │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  │                                 │                                         │  │
│  │                                 ▼                                         │  │
│  │  ┌─────────────────────────────────────────────────────────────────────┐  │  │
│  │  │ 9. ENVELOPE                                                         │  │  │
│  │  │    env = ADSR(phase)                                                │  │  │
│  │  │    sL *= amp * env    sR *= amp * env                               │  │  │
│  │  │    ◄── [env_attack] [env_decay] [env_sustain] [env_release]         │  │  │
│  │  └─────────────────────────────────────────────────────────────────────┘  │  │
│  └─────────────────────────────────┬─────────────────────────────────────────┘  │
└────────────────────────────────────┼────────────────────────────────────────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │   GRAIN SUMMATION   │
                          │   yL = Σ grain_L    │
                          │   yR = Σ grain_R    │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │   PINNA NOTCH       │
                          │   (optional)        │
                          │   Biquad notch L/R  │◄── [pinna_on] [pinna_depth]
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │   OUTPUT SOFT CLIP  │
                          │   tanh(y * 1.2)/1.2 │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │     OUT L / OUT R   │
                          └─────────────────────┘
```

---

## 🔄 2. Control Flow: Event Scheduling

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         GRAIN SCHEDULING FLOW                                   │
└─────────────────────────────────────────────────────────────────────────────────┘

                    ┌─────────────────┐
                    │  samples_to_next│
                    │     counter     │
                    └────────┬────────┘
                             │
                             ▼
              ┌──────────────────────────────┐
              │   counter <= 0 ?             │
              └──────────────┬───────────────┘
                             │ YES
                             ▼
              ┌──────────────────────────────┐
              │   find_free_grain()          │
              │   (scan grains[0..31])       │
              └──────────────┬───────────────┘
                             │
                             ▼
         ┌───────────────────────────────────────┐
         │         WEYL SEQUENCES                │
         │   u1 = w_phi.next()    (durata)       │◄── step = 1/φ
         │   u2 = w_s2.next()     (ampiezza)     │◄── step = 1/√2 or 1/φ²
         │   u3 = w_pl.next()     (pan)          │◄── step = 1/ρ or 1/φ³
         │   u4,u5,u6 = rng.uni01() (extra)      │
         └───────────────────┬───────────────────┘
                             │
                             ▼
         ┌───────────────────────────────────────┐
         │    LATTICE/OU MODULATION              │
         │    (if AUREO_THERMO_LATTICE)          │
         │                                       │
         │    lat_u = probe(w_phi.next())        │◄── [lattice]
         │    oup = ou_pan.y                     │◄── [thermo] [T]
         │    oua = map_phi_range(1/φ, φ, ou_amp)│
         │    oui = ou_itd.y                     │
         └───────────────────┬───────────────────┘
                             │
                             ▼
         ┌───────────────────────────────────────┐
         │    HEMISPHERE COUPLING                │
         │                                       │
         │    coupling = p_hemis_coupling        │◄── [hemis_coupling]
         │    time_weight = f(ratio_gap)         │
         │    hemi = coupling * time_weight      │
         │                                       │
         │    pan = (1-hemi)*pan - hemi*prev_pan │
         │    itd = (1-hemi)*itd - hemi*prev_itd │
         │    ild = (1-hemi)*ild - hemi*prev_ild │
         └───────────────────┬───────────────────┘
                             │
                             ▼
         ┌───────────────────────────────────────┐
         │    POISSON SPATIAL ENFORCEMENT        │
         │    (minimum pan distance check)       │
         │                                       │◄── [spat_min_deg] [spat_min_ms]
         │    Avoid spatial clustering           │
         └───────────────────┬───────────────────┘
                             │
                             ▼
         ┌───────────────────────────────────────┐
         │    SCHEDULE NEXT EVENT                │
         │                                       │
         │    schedule_gap_samples():            │◄── [rate]
         │    λ = rate * (1 + 0.2*sin(2π*t/φ))   │
         │    if burst: λ += 0.3 * hawkes.λ      │◄── [burst]
         │    gap = -ln(U) / λ                   │    (Poisson)
         └───────────────────────────────────────┘
```

---

## 🔗 3. Parameter Interaction Matrix

### Legenda Dipendenze
- **→** Influenza diretta
- **⟷** Influenza bidirezionale
- **⇢** Influenza indiretta/modulata

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│                      PARAMETER DEPENDENCY MATRIX                                 │
├──────────────────┬───────────────────────────────────────────────────────────────┤
│   PARAMETRO      │  INFLUENZA / DIPENDE DA                                       │
├──────────────────┼───────────────────────────────────────────────────────────────┤
│                  │                                                               │
│  rate            │  → samples_to_next, λ (Poisson)                               │
│                  │  ⇢ thermo → ou_rate modula rate ±φ                            │
│                  │  ⇢ burst → hawkes.λ aumenta densità                           │
│                  │                                                               │
│  baselen_ms      │  → grain.dur                                                  │
│                  │  ⟷ len_phi (esponente φ su durata)                            │
│                  │                                                               │
│  len_phi         │  → grain.dur tramite: dur = base * φ^((2u-1)*len_phi)         │
│                  │  ⟷ u1 da w_phi                                                │
│                  │                                                               │
│  width           │  → panL, panR (equal-power)                                   │
│                  │  → pinna_freq_left/right (spread)                             │
│                  │  → binaural.focus, crossfeed                                  │
│                  │                                                               │
│  itd_us          │  → ITD in campioni = itd_us * 1e-6 * sr                       │
│                  │  ⇢ vhs_mod modula ±25%                                        │
│                  │  ⇢ lattice modula ±18%                                        │
│                  │  ⇢ hemis_coupling inverte segno                               │
│                  │                                                               │
│  ild_db          │  → gL, gR gain (db_to_lin)                                    │
│                  │  ⇢ hemis_coupling inverte segno                               │
│                  │                                                               │
│  hemis_coupling  │  → pan, itd, ild (inversione rispetto a precedente)           │
│                  │  ⟷ ratio_gap (tempo tra eventi)                               │
│                  │                                                               │
│  spat_min_deg    │  → Poisson enforcement (min distanza angolare)                │
│  spat_min_ms     │  → Poisson enforcement (min distanza temporale)               │
│                  │                                                               │
│  spat_ipd        │  → ipd_coeff per all-pass decorrelation                       │
│                  │  ⟷ grain.focus                                                │
│                  │                                                               │
│  spat_shadow     │  → shadow_a (LP cutoff contralaterale)                        │
│                  │  ⟷ |pan| (modulazione basata su lateralità)                   │
│                  │                                                               │
│  env_attack      │  → env.attackEnd                                              │
│  env_decay       │  → env.decayEnd                                               │
│  env_sustain     │  → env.sustainLevel                                           │
│  env_release     │  → env.releaseStart                                           │
│                  │  ⟷ gap_samples, dur_samples (dinamico)                        │
│                  │  ⟷ center_distance (|pan|)                                    │
│                  │                                                               │
│  color           │  → noise.color (White/Pink/Brown)                             │
│  color_amt       │  → noise.amount → colorazione filtro                          │
│                  │                                                               │
│  vhs_wow         │  → wowHz (0.1..1.5 Hz) → vhs_mod                              │
│  vhs_flutter     │  → fltHz (7..12 Hz) → vhs_mod                                 │
│                  │  vhs_mod → itd jitter, VhsDrop attenuation                    │
│                  │                                                               │
│  glitch_mix      │  → probabilità GrainKind (VhsDrop/Stutter/Aliaser)            │
│                  │                                                               │
│  srcrush_amt     │  → sr_holdN (1..64 samples hold)                              │
│                  │  ⇢ AUREO_SR_PRIME_SNAP → prime snapping                       │
│                  │                                                               │
│  bitcrush_amt    │  → q_levels (4..16 bit depth)                                 │
│                  │                                                               │
│  seed            │  → RNG state, Weyl sequences x0                               │
│                  │                                                               │
│  pinna_on        │  → abilita pinna notch processing                             │
│  pinna_depth     │  → depth del notch (0..24 dB)                                 │
│                  │  ⟷ width (frequenze notch L/R)                                │
│                  │                                                               │
│  thermo          │  → abilita OU processes                                       │
│                  │  → ou_pan.sigma, ou_itd.sigma, ou_amp.sigma, ou_rate.sigma    │
│                  │  ⟷ T (temperatura stocastica)                                 │
│                  │                                                               │
│  lattice         │  → abilita lattice probe                                      │
│                  │  → lat_u → amp, pan, itd modulation                           │
│                  │                                                               │
│  burst           │  → abilita Hawkes process                                     │
│                  │  → λ boost su eventi cluster                                  │
│                  │                                                               │
│  T               │  → sigma di tutti gli OU processes                            │
│                  │  ⟷ thermo (deve essere attivo)                                │
│                  │                                                               │
│  lat_rate        │  → frequenza step lattice (1..2000 Hz)                        │
│  lat_eps         │  → coupling tra celle lattice                                 │
│  lat_gamma       │  → nonlinearità tanh(γ*x)                                     │
│  lat_sigma       │  → rumore aggiunto al lattice                                 │
│  lat_x/y/z       │  → dimensioni lattice 3D                                      │
│                  │                                                               │
└──────────────────┴───────────────────────────────────────────────────────────────┘
```

---

## 🔀 4. Parameter Groups & Clusters

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         PARAMETER CLUSTERS                                      │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────┐     ┌─────────────────────────┐
│   TIMING & DENSITY      │     │   SPATIAL / BINAURAL    │
│   ─────────────────     │     │   ──────────────────    │
│   • rate                │     │   • width               │
│   • baselen_ms          │────▶│   • itd_us              │
│   • len_phi             │     │   • ild_db              │
│   • thermo → ou_rate    │     │   • hemis_coupling      │
│   • burst               │     │   • spat_min_deg/ms     │
└─────────────────────────┘     │   • spat_ipd            │
           │                    │   • spat_shadow         │
           │                    │   • pinna_on/depth      │
           ▼                    └─────────────────────────┘
┌─────────────────────────┐                │
│   ENVELOPE / SHAPE      │                │
│   ─────────────────     │                ▼
│   • env_attack          │     ┌─────────────────────────┐
│   • env_decay           │     │   TIMBRE / LO-FI        │
│   • env_sustain         │     │   ──────────────        │
│   • env_release         │     │   • color               │
└─────────────────────────┘     │   • color_amt           │
                                │   • vhs_wow             │
┌─────────────────────────┐     │   • vhs_flutter         │
│   STOCHASTIC / THERMO   │     │   • glitch_mix          │
│   ────────────────────  │     │   • srcrush_amt         │
│   • thermo (master)     │     │   • bitcrush_amt        │
│   • lattice             │     └─────────────────────────┘
│   • burst               │
│   • T                   │     ┌─────────────────────────┐
│   • lat_rate            │     │   SYSTEM / INIT         │
│   • lat_eps/gamma/sigma │     │   ─────────────         │
│   • lat_x/y/z           │     │   • seed                │
└─────────────────────────┘     └─────────────────────────┘
```

---

## ⚡ 5. Critical Paths & Performance Analysis

### Complessità per Sample (nel perform loop)

| Operazione | Complessità | Note |
|------------|-------------|------|
| Noise generation | O(1) | Costante, filtro IIR semplice |
| Ring write | O(1) | Singola scrittura |
| LFO update | O(1) | 2 seni |
| Lattice step | O(X×Y×Z) | **Critico** se grande (default 8×8×4=256) |
| OU step | O(4) | 4 processi OU |
| Grain scheduling | O(G) | G = grains attivi (max 32) |
| Per-grain DSP | O(G) | ~15 operazioni per grano |
| Pinna notch | O(1) | 2 biquad |
| Output clip | O(1) | 2 tanh |

### Bottleneck Identificati

1. **Lattice step** (se lat_rate alto e dimensioni grandi)
   - Mitigato da: trylock, step batch
   - Rischio: glitch audio se mutex bloccato

2. **Grain loop** con 32 grani tutti attivi
   - Ogni grano: ITD read (4 accessi ring), pan, ILD, crossfeed, SR/bit crush, IPD, shadow, envelope
   - ~60-80 operazioni floating point per grano

3. **`map_phi_range`** chiamato frequentemente
   - 2× std::log + 1× std::pow per chiamata
   - Potrebbe essere LUT-izzato

---

## 🎯 6. State Machine: Grain Lifecycle

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         GRAIN STATE MACHINE                                     │
└─────────────────────────────────────────────────────────────────────────────────┘

        ┌─────────┐
        │  FREE   │◄──────────────────────────────────────┐
        │ (on=0)  │                                       │
        └────┬────┘                                       │
             │                                            │
             │ find_free_grain() &&                       │
             │ samples_to_next <= 0 &&                    │
             │ poisson_enforce() OK                       │
             │                                            │
             ▼                                            │
        ┌─────────┐                                       │
        │  INIT   │                                       │
        │         │                                       │
        │ • dur = map_len_samples(u1)                     │
        │ • amp = f(u2, lat_u, oua)                       │
        │ • pan = f(u3, lat_u, hemi, prev_pan)            │
        │ • setup ITD, ILD, IPD, shadow                   │
        │ • env = make_envelope(...)                      │
        │ • kind = choose_kind(glitch_mix, u4)            │
        │ • sr_holdN = prime_snap(srcrush)                │
        │ • q_levels = map_bits(bitcrush)                 │
        │ • on = true, age = 0                            │
        └────┬────┘                                       │
             │                                            │
             ▼                                            │
        ┌─────────┐                                       │
        │ ACTIVE  │                                       │
        │ (on=1)  │                                       │
        │         │                                       │
        │ Per sample:                                     │
        │ • phase = age / dur                             │
        │ • env = envelope(phase)                         │
        │ • read ring with ITD                            │
        │ • apply pan, ILD, crossfeed                     │
        │ • apply SR crush, bit crush                     │
        │ • apply glitch kind effect                      │
        │ • apply IPD, shadow                             │
        │ • accumulate to output                          │
        │ • age++                                         │
        └────┬────┘                                       │
             │                                            │
             │ age >= dur                                 │
             │                                            │
             ▼                                            │
        ┌─────────┐                                       │
        │  END    │───────────────────────────────────────┘
        │ on = 0  │
        └─────────┘
```

---

## 📈 7. φ (Golden Ratio) Usage Audit

### Utilizzo Corrente

| Dove | Formula | φ-Coerente? |
|------|---------|-------------|
| `map_phi_range` | `vmin * pow(φ, K*u)` | ✅ |
| `map_len_samples` | `base * pow(φ, exp)` | ✅ |
| `w_phi` step | `1/φ` | ✅ |
| `w_s2` step (default) | `1/√2` | ❌ (non φ) |
| `w_pl` step (default) | `1/ρ` | ❌ (non φ) |
| Rate modulation | `sin(2π * t/φ)` | ✅ |
| Amp mapping | `[1/φ, φ]` | ✅ |
| OU amp | `map_phi_range(1/φ, φ, ...)` | ✅ |
| Hawkes boost | `0.7` | ⚠️ (~1/φ ma non esatto) |
| Time weight | `0.35, 0.65` | ⚠️ (~1/φ², 1/φ ma non esatti) |
| Lattice eps default | `0.18` | ❌ (dovrebbe essere 1/φ³≈0.236) |

### Raccomandazioni

1. Abilitare `AUREO_WD_PHI_POWERS=1` per coerenza
2. Definire costanti:
   ```cpp
   inline constexpr double kPhiSq    = kPhi * kPhi;        // ≈ 2.618
   inline constexpr double kInvPhiSq = 1.0 / kPhiSq;       // ≈ 0.382
   inline constexpr double kInvPhiCu = kInvPhi * kInvPhiSq; // ≈ 0.236
   ```
3. Sostituire magic numbers con costanti φ-based

---

## 🔒 8. Thread Safety Analysis

| Risorsa | Protezione | Rischio |
|---------|------------|---------|
| `report_log` | `report_mu` (trylock) | Basso - eventi persi ma non crash |
| `lat` (Lattice) | `lat_mu` (trylock) | Medio - step saltati possibili |
| `ring` | Nessuna (single writer) | Nessuno - solo DSP thread scrive |
| `grains` | Nessuna (single thread) | Nessuno - solo DSP thread |
| Parametri `p_*` | Atomici impliciti (double) | Basso - eventual consistency OK |

---

## 📋 9. Audit Summary

### Punti di Forza
- Architettura modulare ben separata
- Signal flow chiaro e documentabile
- Uso efficace di Weyl sequences per quasi-randomness
- Sistema di envelope adattivo intelligente
- Thread safety adeguata per contesto audio

### Aree di Miglioramento
- [ ] Coerenza φ incompleta (~60%)
- [ ] Magic numbers non documentati
- [ ] `map_phi_range` potrebbe usare LUT
- [ ] Mancano unit test per aureo_core
- [ ] Documentazione inline insufficiente

### Rischi
- Lattice con dimensioni grandi + lat_rate alto = CPU spike
- 32 grani simultanei con tutti gli effetti = ~2500 ops/sample

---

*Audit generato il 2025-12-08*
