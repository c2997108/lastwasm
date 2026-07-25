#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"
BIN_DIR="$ROOT_DIR/bin"
WEB_DIR="$ROOT_DIR/web"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

if ! command -v em++ >/dev/null 2>&1; then
  echo "em++ was not found. Activate an Emscripten environment first." >&2
  exit 1
fi

echo "[1/3] Cleaning pthread build artifacts..."
make -j"$JOBS" -C "$SRC_DIR" clean >/dev/null || true

echo "[2/3] Building the single-thread asm.js fallback..."
EM_CXXFLAGS="-O3 -g -Wall -std=c++11 -sUSE_ZLIB=1"
EM_LDFLAGS="-sUSE_ZLIB=1 -sWASM=0 -sENVIRONMENT=web,worker -sFORCE_FILESYSTEM -sMODULARIZE=1 -sEXPORT_ES6=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=134217728 -sTOTAL_STACK=8388608 -sINVOKE_RUN=0 -sEXPORTED_FUNCTIONS=['_main'] -sEXPORTED_RUNTIME_METHODS=['FS','callMain'] -sNO_DISABLE_EXCEPTION_CATCHING"

make -j"$JOBS" -C "$SRC_DIR" \
  CXX=em++ CC=emcc \
  CPPFLAGS="" \
  CPPF="-DALPHABET_CAPACITY=66" \
  CXXFLAGS="$EM_CXXFLAGS" \
  LDFLAGS="$EM_LDFLAGS" \
  ../bin/lastal

echo "[3/3] Placing the fallback into web/lastal-asm.js..."
cp "$BIN_DIR/lastal" "$WEB_DIR/lastal-asm.js"
cp "$BIN_DIR/lastal.mem" "$WEB_DIR/lastal.mem"
sed -i 's/[[:blank:]]\+$//' "$WEB_DIR/lastal-asm.js"
echo "Done."
