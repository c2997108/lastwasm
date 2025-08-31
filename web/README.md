# LAST in WebAssembly (Browser)

This folder contains a minimal web UI and build script to compile `lastdb` and `lastal` to WebAssembly using Emscripten, enabling in-browser homology search.

## Prerequisites

- Emscripten SDK installed and activated (`emcc`, `em++` in PATH)
  - https://emscripten.org/docs/getting_started/downloads.html
  - After installing: `source <emsdk>/emsdk_env.sh`

## Build

```bash
cd web
./build_wasm.sh
```

This generates `web/lastdb.js/.wasm` and `web/lastal.js/.wasm`.

## Run

Serve the `web/` directory with a local HTTP server (modules + WASM need HTTP). To enable multi-threading (for `-P > 1`), the page must be cross-origin isolated (COOP/COEP). This repo includes a service worker that sets those headers:

```bash
cd web
python3 -m http.server 8080
# then open http://localhost:8080 in your browser
```

Use the UI to select:
- Target FASTA for `lastdb` (database)
- Query FASTA for `lastal`
- Optional flags for each command

All files run in the browser’s memory (no upload). `lastdb` output files are copied from its virtual filesystem to `lastal`’s virtual filesystem automatically.

## Notes

- Multi-threading: The build enables WebAssembly pthreads. Browsers require `SharedArrayBuffer`, which in turn requires cross-origin isolation (COOP/COEP). The included `coi-serviceworker.js` registers automatically and will reload the page once active; after that `lastal` can use `-P N` (N > 1).
- If opened via `file://` or a server that blocks service workers, threads won’t be available; use `-P 1`.
- zlib is enabled via Emscripten ports (`-sUSE_ZLIB=1`).
- Large datasets may run out of memory in the browser; try smaller inputs.

## Samples

Small example FASTA files are provided in `web/samples/`:

- `web/samples/sample_ref.fasta`
- `web/samples/sample_qry.fasta`

You can use these for a quick sanity check by selecting them for both inputs in the web UI.
