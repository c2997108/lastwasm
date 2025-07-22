#!/bin/sh
# Simple test that the WebAssembly builds of lastal and lastdb run with Node
set -e
cd "$(dirname "$0")"
PATH="../src:$PATH"
DB=/tmp/last-wasm-test
trap 'rm -f "$DB"* out.maf' EXIT
node ../src/lastdb.js "$DB" ttttt.fa
node ../src/lastal.js "$DB" ttttt.fa > out.maf
grep -q ttttt out.maf
printf 'WebAssembly test succeeded\n'

