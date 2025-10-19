cmake_minimum_required(VERSION 3.21)

# Nome del progetto = nome cartella/target per gli script C74
project(aureonoise_tilde CXX)

# ----- Parametri richiesti dagli script del Max SDK -----
# Deve esistere un target con questo nome
set(C74_TARGET "aureonoise_tilde")

# (opzionale) nome file bundle in output (quello che vedi in Max)
# Di solito si vuole "aureonoise~.mxo"
set(C74_EXTERNAL_NAME "aureonoise~")

# macOS target
set(CMAKE_OSX_ARCHITECTURES "arm64")            # per Universal: "arm64;x86_64"
set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0")

# C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build-type di default (se non multi-config come Xcode)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Opzioni prestazionali (tuning)
option(AUREONOISE_NATIVE "Use -mcpu/-march native where available" ON)
option(AUREONOISE_LTO    "Enable Link-Time Optimization in Release" ON)
option(AUREONOISE_STRICT "Extra warnings" ON)

# Percorso Max SDK (puoi esportare MAX_SDK_PATH nell'ambiente)
if(NOT DEFINED MAX_SDK_PATH)
  set(MAX_SDK_PATH "$ENV{HOME}/max-sdk")
endif()

# Percorso di installazione dell'external. Alcuni script del Max SDK lo
# utilizzano durante la fase di configurazione, quindi impostiamo un fallback
# *prima* di includerli in modo da evitare percorsi vuoti come "/externals"
# quando la variabile non è definita nell'ambiente (scenario comune nei
# container CI).
if(NOT DEFINED MAX_PACKAGE_DIR OR MAX_PACKAGE_DIR STREQUAL "")
  # Usa una cartella locale (all'interno dell'albero di build) come fallback
  # sicuro. Viene creata esplicitamente per permettere agli script del Max SDK
  # di scrivere i file temporanei (Info.plist ecc.) senza errori.
  set(MAX_PACKAGE_DIR "${CMAKE_BINARY_DIR}/externals")
endif()
file(MAKE_DIRECTORY "${MAX_PACKAGE_DIR}")

# Script del Max SDK: devono essere inclusi PRIMA di creare il target
include("${MAX_SDK_PATH}/source/max-sdk-base/script/max-pretarget.cmake")

# Core condiviso
add_library(aureo_core INTERFACE)
target_include_directories(aureo_core INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/source")
target_compile_options(aureo_core INTERFACE
  -O3
  -ffast-math
  -fno-math-errno
  -funroll-loops
  -fstrict-aliasing
  -fvisibility=hidden
)

# Target modulo bundle .mxo: il nome DEVE essere ${C74_TARGET}
add_library(${C74_TARGET} MODULE
  "${CMAKE_CURRENT_LIST_DIR}/source/aureonoise~/aureonoise_tilde.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/source/aureonoise_attributes.cpp"
)
target_include_directories(${C74_TARGET} PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/source"
  "${CMAKE_CURRENT_LIST_DIR}/source/aureo_core"
  "${CMAKE_CURRENT_LIST_DIR}/source/aureonoise~"
)
set_target_properties(${C74_TARGET} PROPERTIES
  PREFIX ""
  SUFFIX ".mxo"
  BUNDLE TRUE
  BUNDLE_EXTENSION "mxo"
  OUTPUT_NAME "${C74_EXTERNAL_NAME}"   # => aureonoise~
)
target_link_libraries(${C74_TARGET} PRIVATE aureo_core)
target_link_options(${C74_TARGET} PRIVATE "-Wl,-undefined,dynamic_lookup")

# Include dal Max SDK
if(DEFINED C74_INCLUDES)
  target_include_directories(${C74_TARGET} PRIVATE ${C74_INCLUDES})
endif()

# Flag clang (Apple Silicon) "spinti"
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(${C74_TARGET} PRIVATE -fno-trapping-math)
  if(AUREONOISE_NATIVE)
    # In alternativa: -mcpu=native (richiede Xcode recente)
    target_compile_options(${C74_TARGET} PRIVATE -mcpu=apple-m1 -mtune=apple-m1)
  endif()
  if(AUREONOISE_STRICT)
    target_compile_options(${C74_TARGET} PRIVATE
      -Wall -Wextra
      -Wno-unknown-pragmas -Wno-deprecated-declarations
    )
  endif()
endif()

# LTO (se supportata)
if(AUREONOISE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT _lto_ok OUTPUT _lto_msg)
  if(_lto_ok)
    set_property(TARGET ${C74_TARGET} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
  else()
    message(STATUS "LTO not supported: ${_lto_msg}")
  endif()
endif()

# Post-target del Max SDK (link a MaxAPI/framework ecc.)
# DEVE essere incluso DOPO la creazione del target con lo stesso nome di C74_TARGET.
include("${MAX_SDK_PATH}/source/max-sdk-base/script/max-posttarget.cmake")

# --- Firma automatica su macOS ---
if(APPLE)
  find_program(CODESIGN_EXECUTABLE codesign)
  if(CODESIGN_EXECUTABLE)
    # firma il bundle .mxo (cartella al livello superiore di Contents/MacOS)
    add_custom_command(
      TARGET aureonoise_tilde POST_BUILD
      COMMAND "${CODESIGN_EXECUTABLE}" --force --deep -s - "$<TARGET_FILE_DIR:aureonoise_tilde>/../.."
      COMMENT "Ad-hoc signing aureonoise~.mxo"
      VERBATIM)
  endif()
endif()

# Install destino (cambia a piacere)
install(TARGETS ${C74_TARGET}
  BUNDLE  DESTINATION "${MAX_PACKAGE_DIR}"
  LIBRARY DESTINATION "${MAX_PACKAGE_DIR}")
