#!/usr/bin/env bash
set -euo pipefail

# Build Node.js-targeted WASM for lastal (CLI usage in Node)

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"
BIN_DIR="$ROOT_DIR/bin"
OUT_DIR="$ROOT_DIR/node"

command -v em++ >/dev/null 2>&1 || {
  echo "error: em++ (Emscripten) not found in PATH." >&2
  echo "Install emsdk and run: source <emsdk>/emsdk_env.sh" >&2
  exit 1
}

mkdir -p "$OUT_DIR"

# Determine parallel job count (override with JOBS env)
JOBS=${JOBS:-}
if [ -z "$JOBS" ]; then
  if command -v nproc >/dev/null 2>&1; then
    JOBS=$(nproc)
  elif command -v getconf >/dev/null 2>&1; then
    JOBS=$(getconf _NPROCESSORS_ONLN || echo 4)
  elif command -v sysctl >/dev/null 2>&1; then
    JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
  else
    JOBS=4
  fi
fi
echo "Using parallel build jobs: -j${JOBS}"

echo "[node:1/3] Cleaning previous build artifacts..."
make -j"$JOBS" -C "$SRC_DIR" clean >/dev/null || true

echo "[node:2/3] Building lastal (Node target)..."

if [[ "${DEBUG:-0}" == "1" ]]; then
  EM_OPT="-O0 -g3"
  EM_SAFE="-sASSERTIONS=2 -sSAFE_HEAP=1 -sSTACK_OVERFLOW_CHECK=2"
else
  EM_OPT="-O3 -g"
  EM_SAFE=""
fi

EM_CXXFLAGS="$EM_OPT -Wall -std=c++11 -sUSE_ZLIB=1 -sUSE_PTHREADS=1 -pthread"
EM_LDFLAGS="-sUSE_ZLIB=1 -sUSE_PTHREADS=1 -pthread \
  -sENVIRONMENT=node \
  -sFORCE_FILESYSTEM \
  -sMODULARIZE=1 \
  -sEXPORT_ES6=0 \
  -sINVOKE_RUN=0 \
  -sEXIT_RUNTIME=1 \
  -sALLOW_MEMORY_GROWTH=1 \
  -sINITIAL_MEMORY=268435456 \
  -sSTACK_SIZE=8388608 \
  -sEXPORTED_RUNTIME_METHODS=['FS','callMain'] \
  -sNO_DISABLE_EXCEPTION_CATCHING \
  $EM_SAFE"

make -j"$JOBS" -C "$SRC_DIR" \
  CXX=em++ CC=emcc \
  CPPFLAGS="" \
  CPPF="-DALPHABET_CAPACITY=66 -DHAS_CXX_THREADS" \
  CXXFLAGS="$EM_CXXFLAGS" \
  LDFLAGS="$EM_LDFLAGS" \
  ../bin/lastal

echo "[node:3/3] Placing outputs into node/"
mv -f "$BIN_DIR/lastal" "$OUT_DIR/lastal.js"
if [ -f "$BIN_DIR/lastal.wasm" ]; then
  mv -f "$BIN_DIR/lastal.wasm" "$OUT_DIR/lastal.wasm"
fi
if [ -f "$BIN_DIR/lastal.worker.js" ]; then
  mv -f "$BIN_DIR/lastal.worker.js" "$OUT_DIR/lastal.worker.js"
fi

echo "Done. Files generated in: $OUT_DIR"
echo "- lastal.js / lastal.wasm"
