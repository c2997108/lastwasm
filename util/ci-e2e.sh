#!/usr/bin/env bash
set -euo pipefail

# Optional: provide EMSDK env script via EMSDK_ENV
if [[ -n "${EMSDK_ENV:-}" ]]; then
  if [[ -f "${EMSDK_ENV}" ]]; then
    # shellcheck source=/dev/null
    source "${EMSDK_ENV}"
  else
    echo "warning: EMSDK_ENV points to a non-existent file: ${EMSDK_ENV}" >&2
  fi
fi

if ! command -v emcc >/dev/null 2>&1; then
  echo "error: emcc not found in PATH. Set EMSDK_ENV to emsdk_env.sh or source it before running." >&2
  exit 1
fi

# Install dependencies if needed
if [[ ! -d node_modules ]]; then
  npm ci
fi

# Ensure Playwright Chromium is installed
npx playwright install chromium

# Build WASM (pthreads enabled via build_wasm.sh)
bash web/build_wasm.sh

# Run headless E2E
node web/playwright-test.js

