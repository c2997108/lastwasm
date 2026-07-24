# LAST-JS browser UI

This folder contains the browser implementation used by `index.html`.

The browser path runs the same LAST algorithm as the WASM build without a
dataset-specific fallback. `lastdb.asm.js` and `lastal.asm.js` are generated from
the Emscripten LAST WASM build with Binaryen `wasm2js`. The patched Emscripten
glue in `lastdb.js` and `lastal.js` can load either asm.js or WASM. The page
builds the database once, then runs one `lastal` WASM instance with a shared
linear memory and LAST's pthread pool. Query records are loaded in 64 MiB
batches so `lastal -P` can distribute them across pthreads.

## Files

- `main.js`: UI controller and Plotly dot plot rendering.
- `jslast.js`: high-level runner that writes FASTA into the Emscripten FS, runs
  compiled `lastdb`, copies the generated database files once, and runs compiled
  `lastal -f TAB -P N -i 64M`.
- `jslast-runner-worker.js`: background runner that keeps the browser UI
  responsive while indexing/search is running.
- `plot-parser-worker.js`: parallel TAB parser and Plotly coordinate builder.
- `lastal-worker.js`: asm.js fallback Web Worker entry point.
- `lastdb.js` / `lastal.js`: Emscripten browser glue patched for asm.js loading.
- `lastdb.asm.js` / `lastal.asm.js`: Binaryen `wasm2js` output.
- `../coi-serviceworker.js`: enables COOP/COEP on static hosts so browsers can
  use the shared-memory WASM pthread runtime.
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

Search options are passed to compiled `lastal` after forcing `-f TAB` and adding
`-i 64M` unless the user supplied another batch size. The Threads field controls
LAST's pthread count and is capped by the number of query FASTA records. If
COOP/COEP is active, all pthreads share one WASM linear memory and one database
copy. Without shared memory, the page uses one asm.js search worker instead of
duplicating the database and 256 MiB runtime across workers. The UI reports
pthread query-chunk completion as a search percentage. It also shows elapsed
time and an approximate remaining time, using a small browser-local timing model
from previous runs when available.

For large results, the first 20 KiB of TAB text is painted as soon as search
finishes. The worker waits for that preview to paint before the complete UTF-8
result is transferred as an `ArrayBuffer` and is
available from **Download TAB**, without inserting millions of lines into the
DOM. The page then divides one shared TAB buffer at line boundaries and uses up
to eight Web Workers to parse alignments. Dot plots evenly sample at most
100,000 alignments before Plotly's WebGL `scattergl` renderer draws transferable
coordinate arrays on the GPU. The downloaded TAB result is never sampled.

The log includes `dotplot-timing`, `dotplot-longtasks`, and `dotplot-memory`
profiles. These separate buffer preparation, both parallel parser passes,
message cloning, Plotly's synchronous setup, and the final browser paint.

## Regression

`npm run test:web` generates a fresh TAB baseline with the Node WASM CLI and
compares it byte-for-byte with an 8-thread browser search. This verifies
algorithmic parity without a stored fixture or input-specific output path.

`npm run test:pthread` compares one-thread and eight-thread execution in a
single WASM instance and checks renderer CPU time, speedup, pthread concurrency,
and output line counts.
