# Agent Architecture Design — Aureonoise

**Date:** 2026-03-10
**Status:** Implemented
**Pattern:** Sintetizzatore MediaPack three-layer (CLAUDE.md + .claude/agents/ + .claude/rules/)

## Decision

Approach A — 6 specialists with exclusive file ownership, replicating the sintetizzatore mediapack pattern.

### Alternatives Considered

- **Approach B (4 specialists)**: Merged engine+effects and spatial+noise. Rejected — loses critical separation between spatial (therapeutic core) and noise (spectral correctness).
- **Approach C (3 specialists)**: Rust-core / python-therapeutic / python-app. Rejected — too coarse, single agent would own 18 files with fundamentally different concerns.

## Agent Registry

| Agent | Files Owned | Therapeutic Responsibility |
|-------|-------------|---------------------------|
| engine-specialist | lib.rs, grain.rs, ring.rs, envelope.rs, rng.rs, weyl.rs, math.rs, constants.rs | Real-time safety, zero-allocation |
| spatial-specialist | dialogue.rs, phi_model.rs, binaural.rs, phit.rs | EMDR bilateral, corpus callosum onset <5ms, head model |
| noise-specialist | noise.rs, stoch.rs, burst.rs | Spectral correctness (-1/-2 dB/oct), stochastic resonance |
| effects-specialist | dvf.rs, isochronic.rs, modal.rs, polyrhythm.rs, room.rs, tinnitus.rs, external.rs | 40Hz gamma, tinnitus notch, body resonance safety |
| validation-specialist | validate.py, analysis.py, presets.py, tests/ | Scientific guardian, literature-derived thresholds |
| app-specialist | app.py, audio.py, __init__.py | GUI, session management, safety warnings, volume limiter |

## Three-Layer Architecture

1. **CLAUDE.md** (root): Project overview, architecture key facts, agent registry, code conventions, critical issues, learned rules
2. **.claude/agents/*.md**: Agent definition with YAML frontmatter (name, description, tools, model:opus), Scope (Read/Write, Read-only, NON toccare), Architecture diagrams, Regole non negoziabili
3. **.claude/rules/*.md**: Glob-scoped context rules, auto-injected when touching matching files

## Cross-Agent Constraints

1. Every Rust change must pass `python -m aureonoise.validate --strict`
2. No agent modifies files owned by another agent without explicit declaration
3. Therapeutic thresholds in validate.py require paper citation to modify
4. Zero allocation in process_block() (engine-specialist mandate)
5. phi_model.rs has 4 disconnected functions — spatial-specialist has the wiring mandate

## Priority Implementation Work

### P0 — Safety (non-negotiable)
- Body resonance amplitude limiting (effects-specialist + engine-specialist)
- Epilepsy contraindication warnings (app-specialist + effects-specialist)
- Volume limiter (app-specialist)
- Session dose control (app-specialist + validation-specialist)

### P1 — Correctness
- SpectralTilt brown fix (noise-specialist)
- Envelope Hann/raised-cosine option (engine-specialist)
- Coherence threshold fix in validate.py (validation-specialist)
- EMDR rate 1.5→1.2 Hz (validation-specialist)

### P2 — Completeness
- phi_model.rs wiring (spatial-specialist)
- 0.1 Hz macro-modulation (noise-specialist + engine-specialist)
- BAC feedback loop (spatial-specialist + engine-specialist)
