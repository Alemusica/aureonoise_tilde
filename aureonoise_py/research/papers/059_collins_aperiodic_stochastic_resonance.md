# 059 — Collins 1995: Aperiodic Stochastic Resonance in Excitable Systems

## Metadata
- **Title:** Aperiodic stochastic resonance in excitable systems
- **Authors:** JJ Collins, CC Chow, TT Imhoff
- **Year:** 1995
- **Journal:** Physical Review E
- **Volume/Pages:** 52(4):R3321-R3324
- **DOI:** 10.1103/PhysRevE.52.R3321
- **PubMed:** 9963950
- **URL:** https://journals.aps.org/pre/abstract/10.1103/PhysRevE.52.R3321

## Abstract

Stochastic resonance (SR) is a phenomenon wherein the response of a nonlinear system to a weak periodic input signal is optimized by the presence of a particular level of noise. The authors present a method and theory for characterizing SR-type behavior in excitable systems with aperiodic (non-periodic) inputs. They demonstrate aperiodic stochastic resonance (ASR) in the FitzHugh-Nagumo neuronal model, showing that noise can enhance the transmission of irregular, biologically realistic signals — not just periodic ones.

## Key Findings
- **Coined "aperiodic stochastic resonance" (ASR):** Extended SR theory from periodic-only signals to arbitrary aperiodic inputs
- Demonstrated ASR in the FitzHugh-Nagumo neuronal model
- Noise-enhanced signal detection works for biologically realistic (non-periodic) neural spike trains
- Optimal noise level exists: too little noise doesn't help, too much noise drowns the signal (inverted-U function)
- The cross-correlation between input and output is maximized at an intermediate noise level

## Relevance to Aureonoise

This is the foundational paper for Aureonoise's stochastic resonance mechanism. Aureonoise's granular noise texture is aperiodic by design (Hawkes process timing ensures no two moments are identical). Collins' ASR theory predicts that adding an optimal level of this aperiodic noise to the brain's input can enhance detection of weak, coherent neural signals while leaving isolated/rigid patterns exposed. The optimal noise level (-15 to -20 dB below signal, from auditory psychophysics literature) is the target for Aureonoise's temperature parameter. The OU process in stoch.rs naturally produces mean-reverting noise that matches the aperiodic character Collins describes.

## See Also
- 060 — McDonnell 2009, "What is Stochastic Resonance?" review
- 030 — Vazquez-Rodriguez 2017, SR at criticality in human cortex
- 031 — Matthews 2024, SR in sensory systems
- 032 — van der Groen 2016, transcranial random noise and SR
