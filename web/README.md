# LAST-JS browser UI

This folder contains the browser implementation used by `index.html`.

The browser path runs the same LAST algorithm as the WASM build without a
dataset-specific fallback. `lastdb.asm.js` is generated with Binaryen `wasm2js`.
`lastal.js` loads the pthread WASM runtime, while `lastal-asm.js` contains the
single-thread asm.js fallback generated from the same LAST source. The page
builds the database once, then runs one `lastal` WASM instance with a shared
linear memory and LAST's pthread pool. Query records are loaded in 100 KiB
batches so progress updates throughout large runs. `lastal -P` distributes each
batch across pthreads, though smaller batches can reduce parallel efficiency.

## Files

- `main.js`: UI controller and dot-plot pipeline coordination.
- `regl-dotplot.js`: full-count regl/WebGL2 renderer with canvas axes,
  pan/zoom, fit, and GPU picking.
- `regl.min.js`: locally bundled regl 2.1.1 runtime.
- `jslast.js`: high-level runner that writes FASTA into the Emscripten FS, runs
  compiled `lastdb`, copies the generated database files once, and runs compiled
  `lastal -f TAB -P N -i 64M`.
- `jslast-runner-worker.js`: background runner that keeps the browser UI
  responsive while indexing/search is running.
- `plot-parser-worker.js`: parallel TAB parser and full-count typed-array
  coordinate builder.
- `lastal-worker.js`: asm.js fallback Web Worker entry point.
- `lastdb.js`: Emscripten browser glue for the database builder.
- `lastal.js`: Emscripten browser glue for the pthread WASM runtime.
- `lastal-asm.js`: Emscripten asm.js fallback runtime.
- `lastdb.asm.js`: Binaryen `wasm2js` database builder.
- `build_asm.sh`: regenerates the single-thread asm.js fallback.
- `../coi-serviceworker.js`: enables COOP/COEP on static hosts so browsers can
  use the shared-memory WASM pthread runtime. The page waits for activation and
  reloads automatically before starting the app when isolation is needed.
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

The browser always runs compiled `lastdb` with these LAST-compatible defaults:

```text
--bits=4 -R00 -uNEAR -w1 -W1 -S1 -C1 -v
```

The **Maximum seed hits per query position** field controls `lastal -m` and
defaults to 100. **Query batch size** controls `lastal -i` and defaults to
`64M`. The browser also forces `-f TAB`. The Threads field controls
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
to eight Web Workers to parse alignments. Alignments at or above the configured
configured cutoff are transferred in compact typed arrays and rendered as line
geometry by the local regl/WebGL2 renderer. The cutoff defaults to 1000 bp, so
alignments of 1000 bp or less are hidden. It only affects the plot, not TAB
output. The plot is not sampled. CPU metadata remains in typed
arrays for hover details, while positions and selection IDs are uploaded to GPU
buffers. To retain the former Plotly view's readability, a deterministic set of
up to 100,000 alignments is emphasized while every remaining alignment is still
drawn at lower opacity. A canvas overlay draws axes and dotted sequence
boundaries without creating per-alignment DOM or JavaScript objects.
By default, reference and query axes follow the sequence order in the original
FASTA files, including sequences without plotted alignments. Clearing **Plot in
FASTA sequence order** restores the alignment-output order.

The log includes `dotplot-timing`, `dotplot-longtasks`, and `dotplot-memory`
profiles. These separate buffer preparation, both parallel parser passes,
message cloning, regl buffer upload, draw submission, GPU completion, and the
final browser paint.

## Regression

`npm run test:web` generates a fresh TAB baseline with the Node WASM CLI and
compares it byte-for-byte with an 8-thread browser search. This verifies
algorithmic parity without a stored fixture or input-specific output path.

`npm run test:pthread` compares one-thread and eight-thread execution in a
single WASM instance and checks renderer CPU time, speedup, pthread concurrency,
and output line counts.
