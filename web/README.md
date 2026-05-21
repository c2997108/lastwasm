# LAST-JS browser UI

This folder contains the browser implementation used by `index.html`.

The browser path runs the same LAST algorithm as the WASM build without a
dataset-specific fallback. `lastdb.asm.js` and `lastal.asm.js` are generated from
the Emscripten LAST WASM build with Binaryen `wasm2js`. The patched Emscripten
glue in `lastdb.js` and `lastal.js` can load either asm.js or WASM. The page
builds the database once, then runs `lastal` search jobs in a Web Worker pool by
splitting the query FASTA by record.

## Files

- `main.js`: UI controller and Plotly dot plot rendering.
- `jslast.js`: high-level runner that writes FASTA into the Emscripten FS, runs
  compiled `lastdb`, copies the generated database files, and dispatches
  compiled `lastal -f TAB` search jobs.
- `jslast-runner-worker.js`: background runner that keeps the browser UI
  responsive while indexing/search is running.
- `lastal-worker.js`: Web Worker entry point for parallel `lastal` search.
- `lastdb.js` / `lastal.js`: Emscripten browser glue patched for asm.js loading.
- `lastdb.asm.js` / `lastal.asm.js`: Binaryen `wasm2js` output.
- `../coi-serviceworker.js`: enables COOP/COEP on static hosts so browsers can
  use the faster WASM worker runtime.
- `samples/`: small FASTA files for smoke tests.

## Run

Serve the repository root and open `index.html`:

```bash
python -m http.server 8000
```

Then open `http://localhost:8000/`.

Opening `index.html` directly with `file://` is not supported because browsers
block ES module loading from local files.

## Options

The UI keeps LAST-like option fields. `lastdb` options are passed to compiled
LAST; when the "LAST defaults" checkbox is enabled, the browser uses:

```text
--bits=4 -R00 -uNEAR -w1 -W1 -S1 -C1 -v
```

Search options are passed to compiled `lastal` after forcing `-f TAB`. The
Workers field controls the Web Worker pool size. The requested worker count is
capped by the number of query FASTA records; a single query record stays on one
search worker to preserve LAST's query-level behavior. If COOP/COEP is active,
search workers use the WASM runtime; otherwise they use the asm.js runtime.
The UI shows elapsed time and an approximate remaining time. The first run
estimates after worker completions are observed; later runs also use a small
browser-local timing model from previous runs.

## Regression

`npm run test:web` generates a fresh TAB baseline with the Node WASM CLI and
compares it byte-for-byte with an 8-worker browser search. This verifies
algorithmic parity without a stored fixture or input-specific output path.
