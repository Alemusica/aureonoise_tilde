# aureonoise~ Project Analysis

## Overview

aureonoise~ è un external audio per Cycling '74 Max progettato per generare trame di rumore e glitch stereofoniche basate su sequenze quasi-periodiche e modelli stocastici controllati. Il codice sorgente è organizzato attorno a un modulo core riutilizzabile (`aureo_core`) e al file principale dell'external (`aureonoise~/aureonoise_tilde.cpp`) che collega la logica DSP all'SDK di Max.【F:README.md†L1-L37】【F:source/aureonoise~/aureonoise_tilde.cpp†L1-L81】

## Build system

Il progetto usa CMake (>=3.21) e richiede il Max SDK 8.5+. Il `CMakeLists.txt` configura un target `MODULE` senza prefisso che produce il bundle `aureonoise~.mxo`, include gli script `max-pretarget.cmake`/`max-posttarget.cmake` e definisce opzioni per ottimizzazioni native, LTO e warning estesi. La configurazione gestisce anche l'installazione opzionale nel percorso `MAX_PACKAGE_DIR` ed è pensata per macOS (architettura `arm64`, deployment target 11.0).【F:CMakeLists.txt†L1-L95】【F:CMakeLists.txt†L109-L136】

## Struttura del codice

- `source/aureonoise~/aureonoise_tilde.cpp`: implementa il lifecycle dell'oggetto Max (costruzione, liberazione, assistenza, reset) e il perform DSP a 64 bit. Gestisce sequenziamento dei grani, modulazioni wow/flutter, miscelazione binaurale/pinna, gestione del ring buffer e scheduling pseudo-poissoniano degli eventi con modulazioni OU e lattice stocastico opzionali.【F:source/aureonoise~/aureonoise_tilde.cpp†L20-L192】【F:source/aureonoise~/aureonoise_tilde.cpp†L193-L304】
- `source/aureonoise~/aureonoise_state.hpp`: definisce la struttura `t_aureonoise` con parametri d'istanza, stati dei moduli (RNG, sequenze Weyl, ring buffer, filtri pinna, Ornstein-Uhlenbeck, lattice termico, processo di Hawkes) e helper inline per aggiornare i filtri pinna.【F:source/aureonoise~/aureonoise_state.hpp†L1-L123】【F:source/aureonoise~/aureonoise_state.hpp†L124-L196】
- `source/aureonoise_attributes.cpp`: dichiara gli attributi Max esposti all'utente (rate, durate, parametri binaurali, rumore, effetti lo-fi, seed, opzioni termiche) e fornisce utility come la generazione dei primi da usare nello scheduling dei grani.【F:source/aureonoise_attributes.cpp†L1-L108】【F:source/aureonoise_attributes.cpp†L109-L196】
- `source/aureo_core/`: contiene le primitive DSP e matematiche condivise, come costanti e flag di configurazione (`aureo_config.hpp`), RNG xoroshiro-like (`aureo_rng.hpp`), sequenze di Weyl, generatori di rumore colorato, filtri binaurali, ADSR field shaping, anelli di delay e dinamiche stocastiche (OU, lattice, Hawkes).【F:source/aureo_core/aureo_config.hpp†L1-L35】【F:source/aureo_core/aureo_rng.hpp†L1-L17】

## Flusso audio e caratteristiche DSP

Durante il perform 64-bit il codice:
1. Sincronizza lo stato pinna, ring buffer e modulazioni wow/flutter, con eventuale smoothing dei parametri.【F:source/aureonoise~/aureonoise_tilde.cpp†L205-L249】
2. Aggiorna il processo lattice/OU quando abilitato per modulare pan, ampiezza, ITD e densità di eventi.【F:source/aureonoise~/aureonoise_tilde.cpp†L252-L289】
3. Genera rumore colorato processato via tanh, scrive nel ring buffer e crea nuovi grani quando il contatore raggiunge zero, impostando durata, ampiezza e pan sulla base di sequenze Weyl e modulazioni stocastiche.【F:source/aureonoise~/aureonoise_tilde.cpp†L290-L329】
4. Applica mix binaurale equal-power con filtri notch della pinna, micro ITD/ILD e controlli lo-fi (sample/bit crush, glitch mix) prima di scrivere nei due canali di uscita.【F:source/aureonoise~/aureonoise_tilde.cpp†L330-L392】

## Parametri utente

Gli attributi espongono un set ricco di controlli: densità di eventi (`rate`, ora estesa fino a 120 Hz), durata base e varianza φ (`baselen_ms`, `len_phi`, con durata base riducibile a 1 ms), spazializzazione (`width`, `itd_us`, `ild_db`, `pinna_on`, `pinna_depth`), carattere timbrico (`color`, `color_amt`, `vhs_wow`, `vhs_flutter`, `glitch_mix`, `srcrush_amt`, `bitcrush_amt`), seed, e – se compilato con `AUREO_THERMO_LATTICE` – parametri termodinamici/lattice per modulazioni complesse inclusi `thermo`, `lattice`, `burst`, `T`, `lat_rate`, `lat_eps`, `lat_gamma`, `lat_sigma`, `lat_x/y/z`.【F:source/aureonoise_attributes.cpp†L37-L196】

## Dipendenze esterne

L'inclusione di `max-pretarget.cmake`/`max-posttarget.cmake` implica una dipendenza diretta dal Max SDK. Non sono presenti altre librerie di terze parti: tutte le utilità DSP/stocastiche sono implementate internamente in `aureo_core` e collegate come libreria `INTERFACE` in CMake.【F:CMakeLists.txt†L32-L77】【F:CMakeLists.txt†L78-L108】

## Estendibilità e punti di attenzione

- Il design modulare separa chiaramente la logica Max dallo strato DSP riutilizzabile, rendendo possibile usare `aureo_core` in altri contesti C++17.
- Le opzioni di build prevedono flag aggressivi (`-Ofast`, `-ffast-math`, `-funroll-loops`, LTO) che richiedono test accurati su target differenti, soprattutto se si vuole portare il progetto su piattaforme non-Apple o con ABI differenti.【F:CMakeLists.txt†L44-L77】
- Il codice sfrutta SSE2 e path condizionali ARM (`__aarch64__`), quindi bisogna assicurarsi di gestire fallback adeguati quando si compila su target privi di tali estensioni.【F:source/aureonoise~/aureonoise_tilde.cpp†L10-L18】【F:source/aureonoise~/aureonoise_tilde.cpp†L236-L246】
- La dimensione del ring buffer (`131072` campioni) e il numero massimo di grani (`32`) sono configurati in `aureo_config.hpp`; eventuali modifiche impattano memoria e scheduling e richiedono ricalibrazione delle funzioni di mapping durate/amplitudini.【F:source/aureo_core/aureo_config.hpp†L23-L33】

## Suggerimenti operativi

- Prima di costruire il progetto assicurarsi che `MAX_SDK_PATH` punti alla root del Max SDK e che l'architettura target corrisponda alla piattaforma finale.
- Usare l'attributo `clear` in Max per resettare rapidamente lo stato (ring buffer, grani, modulazioni) durante il debugging di preset complessi.【F:source/aureonoise~/aureonoise_tilde.cpp†L167-L182】
- Il messaggio Max `report` invia sull'outlet di destra un log strutturato di ogni grano (timestamp, ampiezza, pan, valori Weyl/RNG, stato lattice/burst) con statistiche di densità per verificare l'effettiva stochasticità del motore aureo.【F:source/aureonoise~/aureonoise_tilde.cpp†L240-L325】【F:source/aureonoise~/aureonoise_tilde.cpp†L333-L373】
- Per analizzare manualmente il comportamento stocastico dei grani si possono ancora osservare `samples_to_next`, i valori delle sequenze Weyl (`w_phi`, `w_s2`, `w_pl`) e l'output del lattice termico quando `AUREO_THERMO_LATTICE` è attivo.

