# 068 — Kanda 2023: Hawkes Process Modeling of Cortical Neuron Firing

## Metadata
- **Title:** Hawkes process modeling quantifies complicated firing behaviors of cortical neurons during sleep and wakefulness
- **Authors:** Takeshi Kanda, Toshimitsu Aritake, Kaoru Ohyama, Kaspar E Vogt, Yuichi Makino, Thomas J McHugh, Hideitsu Hino, Shotaro Akaho, Noboru Murata
- **Year:** 2023
- **Journal:** bioRxiv (preprint)
- **DOI:** 10.1101/2023.07.29.550297
- **URL:** https://www.biorxiv.org/content/10.1101/2023.07.29.550297v1

## Abstract

The study applies Hawkes process modeling to quantify complicated firing behaviors of cortical neurons during sleep and wakefulness. The Hawkes process models sequential random events exhibiting temporal clusters — events that self-excite, where one firing increases the probability of subsequent firings within a time window. The magnitude of firing intensity was inversely proportional to the time constant of firing intensity. Non-REM sleep increased the magnitude of firing intensity and decreased the time constant of firing intensity. Repetitive firing is ordered to become high frequency and short term during non-REM sleep, while unregulated components of firing are independent of the sleep/wake state in the cortex.

## Key Findings
- **Hawkes process accurately models cortical neuron firing patterns** — capturing temporal clustering and self-excitation
- Neural firing exhibits self-exciting dynamics: each spike increases the probability of subsequent spikes (the hallmark of Hawkes processes)
- Sleep and wakefulness produce distinct Hawkes process parameters
- Non-REM sleep: higher magnitude, shorter time constant (sharp, intense bursts)
- Wake state: lower magnitude, longer time constant (more sustained, lower-frequency firing)

## Relevance to Aureonoise

This paper validates Aureonoise's use of the Hawkes process (burst.rs) for grain timing. Neural spike trains ARE Hawkes processes — the brain's own temporal structure follows self-exciting point process dynamics. By using a Hawkes process to schedule grain events, Aureonoise produces temporal patterns that match the brain's native temporal statistics. The sleep/wake parameter differences suggest that Aureonoise could adjust Hawkes parameters for different therapeutic modes: higher base_rate + shorter time constant for sleep/delta presets, lower base_rate + longer time constant for awake/alert presets.

## See Also
- 044 — Zhou 2012, pink noise and brain synchronization
- 045 — Hesse 2014, 1/f neural noise and criticality
