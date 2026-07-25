#!/usr/bin/env bash
set -euo pipefail

# Build LAST lastal with pthread support using Emscripten.
# Outputs: web/lastal.js, web/lastal.wasm, and web/lastal.worker.js

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"
BIN_DIR="$ROOT_DIR/bin"
WEB_DIR="$ROOT_DIR/web"

# Respect emsdk environment if active; otherwise, use a project-local cache.
EMCC_BIN="$(command -v emcc || true)"
if echo "$EMCC_BIN" | grep -q "/emsdk/"; then
  # Use emsdk's default cache/config
  :
else
  if [ -z "${EM_CACHE:-}" ]; then
    export EM_CACHE="$ROOT_DIR/.emscripten_cache"
  fi
  unset EM_FROZEN_CACHE
  mkdir -p "$EM_CACHE"
fi

command -v em++ >/dev/null 2>&1 || {
  echo "error: em++ (Emscripten) not found in PATH." >&2
  echo "Install emsdk and run: source <emsdk>/emsdk_env.sh" >&2
  exit 1
}

mkdir -p "$WEB_DIR"

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

# Ensure a writable EM_CONFIG and cache (only if not already provided by emsdk)
if ! echo "$EMCC_BIN" | grep -q "/emsdk/"; then
  if [ -z "${EM_CONFIG:-}" ]; then
    export EM_CONFIG="$ROOT_DIR/.emscripten"
  fi
fi
if [ -n "${EM_CONFIG:-}" ] && [ "$EM_CONFIG" = "$ROOT_DIR/.emscripten" ] && [ ! -f "$EM_CONFIG" ]; then
  echo "Generating local Emscripten config at $EM_CONFIG"
  emcc --generate-config >/dev/null
  # The above writes to the default location; move it to our desired path if needed
  if [ ! -f "$EM_CONFIG" ]; then
    # Try to copy from env-specified location
    if [ -n "${EM_CONFIG:-}" ] && [ -f "${EM_CONFIG}" ]; then
      cp -f "${EM_CONFIG}" "$ROOT_DIR/.emscripten"
    fi
  fi
fi
# Patch config to use project-local cache and disable frozen cache
if [ -n "${EM_CONFIG:-}" ] && [ "$EM_CONFIG" = "$ROOT_DIR/.emscripten" ] && [ -f "$EM_CONFIG" ]; then
  python3 - "$EM_CONFIG" "$EM_CACHE" << 'PY'
import io,sys,os,re
cfg, cache = sys.argv[1], sys.argv[2]
txt = open(cfg,'r',encoding='utf-8').read()
def set_kv(s, key, value):
    pat = re.compile(r'^(\s*%s\s*=).*$' % re.escape(key), re.M)
    repl = r"\1 %s" % value
    if pat.search(s):
        return pat.sub(repl, s)
    else:
        return s + "\n%s = %s\n" % (key, value)
txt = set_kv(txt, 'FROZEN_CACHE', 'False')
txt = set_kv(txt, 'CACHE', repr(cache))
open(cfg,'w',encoding='utf-8').write(txt)
print('Patched EM_CONFIG:', cfg)
PY
fi

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
make -j"$JOBS" -C "$SRC_DIR" clean >/dev/null || true

echo "[3/4] Building with Emscripten (this may take a while)..."
# Notes:
# - Enable pthreads for multi-CPU (-P > 1) in the browser.
#   This requires cross-origin isolation (COOP/COEP) at runtime.
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

# Emscripten recommends passing -pthread to all compilation units when using threads
EM_CXXFLAGS="$EM_OPT -Wall -std=c++11 -sUSE_ZLIB=1 -sUSE_PTHREADS=1 -pthread"
EM_LDFLAGS="-sUSE_ZLIB=1 -sUSE_PTHREADS=1 -pthread -sENVIRONMENT=web,worker -sFORCE_FILESYSTEM -sMODULARIZE=1 -sEXPORT_ES6=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=134217728 -sTOTAL_STACK=8388608 -sINVOKE_RUN=0 -sEXPORTED_FUNCTIONS=['_main'] -sEXPORTED_RUNTIME_METHODS=['FS','callMain'] -sNO_DISABLE_EXCEPTION_CATCHING $EM_SAFE"

make -j"$JOBS" -C "$SRC_DIR" \
  CXX=em++ CC=emcc \
  CPPFLAGS="" \
  CPPF="-DALPHABET_CAPACITY=66 -DHAS_CXX_THREADS" \
  CXXFLAGS="$EM_CXXFLAGS" \
  LDFLAGS="$EM_LDFLAGS" \
  ../bin/lastal

echo "[4/4] Placing outputs into web/"
cd "$ROOT_DIR"

# Emscripten names the output after the -o path in the makefile (../bin/lastal).
# Copy and give .js extensions for clarity
if [ -f "$BIN_DIR/lastal" ]; then
  mv -f "$BIN_DIR/lastal" "$WEB_DIR/lastal.js"
fi
if [ -f "$BIN_DIR/lastal.wasm" ]; then
  mv -f "$BIN_DIR/lastal.wasm" "$WEB_DIR/lastal.wasm"
fi
if [ -f "$BIN_DIR/lastal.worker.js" ]; then
  mv -f "$BIN_DIR/lastal.worker.js" "$WEB_DIR/lastal.worker.js"
fi
node "$WEB_DIR/patch-pthread-runtime.js" "$WEB_DIR/lastal.js" "$WEB_DIR/lastal.worker.js"
sed -i 's/[[:blank:]]\+$//' "$WEB_DIR/lastal.js"

echo "Done. Files generated in: $WEB_DIR"
echo "- lastal.js / lastal.wasm"
