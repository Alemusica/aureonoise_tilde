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

Il progetto cerca automaticamente il Max SDK leggendo la variabile d'ambiente `MAX_SDK_PATH`. Se non è impostata, viene usato il percorso predefinito `~/max-sdk`. Imposta la variabile prima di aprire VS Code oppure configuralo tramite le impostazioni di CMake Tools:

```bash
export MAX_SDK_PATH="/percorso/al/max-sdk"
```

Puoi anche definire il percorso di installazione del bundle impostando `MAX_PACKAGE_DIR` (di default `~/Documents/Max 9/Packages/MyDev/externals`).

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

## Risoluzione problemi

- **Errore include Max SDK**: indica che il percorso `MAX_SDK_PATH` non punta alla cartella che contiene `source/max-sdk-base/script/max-pretarget.cmake`.
- **Kit mancante**: installa Xcode Command Line Tools o il compilatore clang/gcc adeguato e riavvia VS Code.
- **Firma/notarizzazione**: se distribuisci l'external su macOS, potrebbero essere necessari passaggi aggiuntivi per la firma del bundle.
