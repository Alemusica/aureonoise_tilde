# aureonoise~ CMake Project

Questo repository contiene il progetto CMake per l'external `aureonoise~` destinato a Cycling '74 Max. Il target principale è il bundle `.mxo` prodotto a partire dai sorgenti presenti in `source/`.

## Prerequisiti

1. **Max SDK 8.5+** – Clona o scarica il Max SDK e prendi nota del percorso locale. Il progetto assume che sia disponibile lo script CMake `max-pretarget.cmake`.
2. **CMake 3.21 o successivo** – Il file `CMakeLists.txt` imposta la versione minima richiesta.
3. **Compilatore C++17** – Su macOS è richiesto `clang` (Apple), mentre su altre piattaforme può essere usato `clang`/`gcc` compatibile.
4. **Visual Studio Code** con le estensioni:
   - [CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)
   - [C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)

## Configurazione dell'ambiente

Il progetto cerca automaticamente il Max SDK leggendo la variabile d'ambiente `MAX_SDK_PATH`. Se non è impostata, viene usato il percorso predefinito `~/max-sdk`. Lo stesso valore viene usato anche come fallback per `MAX_SDK_ROOT`, quindi nella maggior parte dei casi non è necessario configurare altre variabili.

Se preferisci impostare esplicitamente il percorso (ad esempio quando lo hai installato in una posizione non standard), definisci `MAX_SDK_PATH` prima di aprire VS Code oppure configuralo tramite le impostazioni di CMake Tools:

```bash
export MAX_SDK_PATH="/percorso/al/max-sdk"
```

Se necessario puoi impostare manualmente `MAX_SDK_ROOT` nello stesso modo (variabile d'ambiente o opzione CMake) per puntare a una copia diversa dell'SDK. Puoi anche definire il percorso di installazione del bundle impostando `MAX_PACKAGE_DIR` (di default `~/Documents/Max 9/Packages/MyDev/externals`).

## Build con Visual Studio Code

1. **Apri la cartella** – Avvia VS Code e scegli `File > Open Folder...`, selezionando la radice del repository.
2. **Seleziona il Kit di compilazione** – CMake Tools mostrerà una notifica "Select a Kit". Scegli il toolchain appropriato (ad es. "Clang 15.0.0" su macOS Apple Silicon). Puoi cambiare kit in seguito dal comando `CMake: Select a Kit`.
3. **Configura la variante di build** – Con il comando `CMake: Select Variant` scegli `Release` (è il tipo predefinito definito dal progetto). Questa variante abilita LTO se supportata e setta le opzioni di ottimizzazione.
4. **Configura il progetto** – CMake Tools lancerà automaticamente `cmake --configure`. Se l'SDK è posizionato correttamente, vedrai la creazione della cartella `build/` con gli script generati per l'external `aureonoise~`. In caso di errore "include could not find ... max-pretarget.cmake" verifica `MAX_SDK_PATH`.
5. **Compila il target** – Usa il comando `CMake: Build` (oppure il pulsante "Build" nella barra di stato). Verrà costruito il bundle `aureonoise~.mxo` usando il target `aureonoise_tilde` definito nel `CMakeLists.txt`.
6. **Installa l'external (opzionale)** – Per copiare automaticamente il bundle nella cartella dei Packages, esegui `CMake: Install`. Il bundle verrà copiato nella directory `MAX_PACKAGE_DIR`.

## Suggerimenti aggiuntivi

- Per forzare un rebuild pulito puoi eseguire `CMake: Delete Cache and Reconfigure`.
- Se stai sviluppando su macOS con architettura differente, modifica `CMAKE_OSX_ARCHITECTURES` nel `CMakeLists.txt` (es. `"arm64;x86_64"` per Universal Binary).
- Assicurati di lanciare Cycling '74 Max dopo l'installazione per verificare che l'external venga caricato senza errori.
- La spazializzazione binaurale sintetica sfrutta ora un modello Woodworth/Gaussian crossfeed con report dettagliato degli azimuth, utile da consultare tramite il messaggio `report` di aureonoise~.

## Inviluppo complesso "vellutato"

L'inviluppo dei grani non è più un semplice ADSR statico: ogni grano evolve in un piano complesso con modulazioni lente e non ripetitive. Il modulo calcola un valore di ampiezza a partire dalla magnitudine del vettore complesso, garantendo un decadimento morbido (*vellutato*) ma sempre leggermente diverso. Gli attributi seguenti permettono di plasmare il comportamento stocastico senza introdurre ciclicità percepibili:

- `env_complex_magvar` (0..1.5) controlla la varianza della magnitudine. Valori più alti rendono il decadimento più mobile.
- `env_complex_phasevar` (0..1.5 rad) regola la diffusione della fase complessa. Incrementandolo aumenta la micro-rotazione tra campioni.
- `env_complex_corr` (0..0.999) imposta la correlazione tra i campioni di rumore usati per la deriva. Valori vicini a 1 producono traiettorie molto levigate.
- `env_complex_tau_ms` (1..2000 ms) definisce la costante di tempo del filtro che insegue il bersaglio ADSR. È espressa in millisecondi e permette di rallentare/accelerare la risposta.
- `env_complex_bias` (0..1) miscela il target istantaneo con il livello di sustain, mantenendo una presenza vellutata anche su code lunghissime.

Le impostazioni predefinite riproducono un decadimento naturale; piccole variazioni dei parametri consentono di ottenere texture differenti mantenendo la non-periodicità aureo.

## Modalità di stimolo

L'attributo `mode` consente di scegliere tra la generazione **texture** classica e la nuova modalità **dichotic** pensata per paradigmi percettivi L/R estremamente lenti (0.1–0.2 Hz). In modalità dichotic:

- I burst vengono alternati automaticamente fra orecchio sinistro e destro rispettando l'intervallo impostato (`dichotic_rate`).
- L'inviluppo è configurabile tramite `dichotic_burst_ms`, `dichotic_attack_ms` e `dichotic_release_ms`, con ampiezza controllata da `dichotic_amp`.
- È possibile programmare trial *match/mismatch* con probabilità `dichotic_match_prob`, scegliendo per ciascun caso se riprodurre **rumore** o **toni** (`dichotic_content_match` / `dichotic_content_mismatch`) e, per i toni, le frequenze di riferimento (`dichotic_tone_match_hz`, `dichotic_tone_mismatch_hz`).
- Ogni burst produce un messaggio diagnostico sull'outlet info (`dichotic index …`) e viene tracciato dal comando `report`, che riporta anche il rapporto match/mismatch e lo stato del buffer.

Quando `mode` è impostato su `texture` le elaborazioni classiche (grani generativi, pinna/crossfeed, ecc.) restano invariate; passando a `dichotic` vengono disattivate automaticamente le sezioni non pertinenti.

## Risoluzione problemi

- **Errore include Max SDK**: indica che il percorso `MAX_SDK_PATH` non punta alla cartella che contiene `source/max-sdk-base/script/max-pretarget.cmake`.
- **Kit mancante**: installa Xcode Command Line Tools o il compilatore clang/gcc adeguato e riavvia VS Code.
- **Firma/notarizzazione**: se distribuisci l'external su macOS, potrebbero essere necessari passaggi aggiuntivi per la firma del bundle.
