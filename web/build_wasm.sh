#!/usr/bin/env bash
set -euo pipefail

# Build LAST (lastdb, lastal) to WebAssembly using Emscripten.
# Outputs: web/lastdb.js(.wasm), web/lastal.js(.wasm)

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"
BIN_DIR="$ROOT_DIR/bin"
WEB_DIR="$ROOT_DIR/web"

command -v em++ >/dev/null 2>&1 || {
  echo "error: em++ (Emscripten) not found in PATH." >&2
  echo "Install emsdk and run: source <emsdk>/emsdk_env.sh" >&2
  exit 1
}

mkdir -p "$WEB_DIR"

echo "[1/4] Generating data headers if missing..."
cd "$SRC_DIR"
if [ ! -f CyclicSubsetSeedData.hh ]; then
  ../build/seed-inc.sh ../data/*.seed > CyclicSubsetSeedData.hh
fi
if [ ! -f GeneticCodeData.hh ]; then
  ../build/gc-inc.sh ../data/gc.prt > GeneticCodeData.hh
fi
if [ ! -f ScoreMatrixData.hh ]; then
  ../build/mat-inc.sh ../data/*.mat > ScoreMatrixData.hh
fi

# Provide a simple version.hh if git describe is unavailable under emmake
if [ ! -f version.hh ]; then
  echo 'unknown' > version.hh
fi

echo "[2/4] Cleaning previous build artifacts..."
make -C "$SRC_DIR" clean >/dev/null || true

echo "[3/4] Building with Emscripten (this may take a while)..."
# Notes:
# - Remove -msse4 and -pthread; Emscripten targets web without threads by default.
# - Enable zlib, filesystem, modularized ES6 exports, and memory growth.
# - The src/makefile produces ../bin/lastdb and ../bin/lastal; with em++ these will be JS glue files.

# Tunable debug flags: set DEBUG=1 to enable heavy checks
if [[ "${DEBUG:-0}" == "1" ]]; then
  EM_OPT="-O0 -g3"
  EM_SAFE="-sASSERTIONS=2 -sSAFE_HEAP=1 -sSTACK_OVERFLOW_CHECK=2"
else
  EM_OPT="-O3 -g"
  EM_SAFE=""
fi

EM_CXXFLAGS="$EM_OPT -Wall -std=c++11 -sUSE_ZLIB=1"
EM_LDFLAGS="-sUSE_ZLIB=1 -sENVIRONMENT=web -sFORCE_FILESYSTEM -sMODULARIZE=1 -sEXPORT_ES6=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=134217728 -sSTACK_SIZE=8388608 -sEXPORTED_RUNTIME_METHODS=['FS','callMain'] -sNO_DISABLE_EXCEPTION_CATCHING $EM_SAFE"

make -C "$SRC_DIR" \
  CXX=em++ CC=emcc \
  CPPFLAGS="" \
  CPPF="-DALPHABET_CAPACITY=66" \
  CXXFLAGS="$EM_CXXFLAGS" \
  LDFLAGS="$EM_LDFLAGS" \
  all

echo "[4/4] Placing outputs into web/"
cd "$ROOT_DIR"

# Emscripten names outputs after the -o path in the makefile (../bin/lastdb, ../bin/lastal)
# Copy and give .js extensions for clarity
if [ -f "$BIN_DIR/lastdb" ]; then
  mv -f "$BIN_DIR/lastdb" "$WEB_DIR/lastdb.js"
fi
if [ -f "$BIN_DIR/lastdb.wasm" ]; then
  mv -f "$BIN_DIR/lastdb.wasm" "$WEB_DIR/lastdb.wasm"
fi
if [ -f "$BIN_DIR/lastal" ]; then
  mv -f "$BIN_DIR/lastal" "$WEB_DIR/lastal.js"
fi
if [ -f "$BIN_DIR/lastal.wasm" ]; then
  mv -f "$BIN_DIR/lastal.wasm" "$WEB_DIR/lastal.wasm"
fi

echo "Done. Files generated in: $WEB_DIR"
echo "- lastdb.js / lastdb.wasm"
echo "- lastal.js / lastal.wasm"
