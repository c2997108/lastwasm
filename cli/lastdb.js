#!/usr/bin/env node
/*
 Node.js CLI wrapper for LAST's lastdb compiled to WebAssembly.

 Usage: cli/lastdb.js [lastdb args]
 It will run in the current working directory and read/write files on the host filesystem.
*/

const path = require('path');
const fs = require('fs');

async function main() {
  const argv = process.argv.slice(2);
  if (argv.length === 0 || argv.includes('-h') || argv.includes('--help')) {
    // Let lastdb itself print usage; just run with -h
    argv.length = 0;
    argv.push('-h');
  }

  const runtimePath = path.resolve(__dirname, '..', 'node');
  const wasmPath = path.join(runtimePath, 'lastdb.wasm');
  const jsPath = path.join(runtimePath, 'lastdb.js');

  if (!fs.existsSync(jsPath)) {
    console.error('error: Node-targeted lastdb.js not found.');
    console.error('hint: run: bash node/build_node_wasm.sh');
    process.exit(2);
  }

  const createModule = require(jsPath);

  // Track exit code from the WASM program
  let exitCode = 0;

  // Infer desired pthread pool size from -P option
  let poolSize = 0;
  for (let i = 0; i < argv.length - 1; i++) {
    if (argv[i] === '-P' && /^\d+$/.test(argv[i + 1])) {
      poolSize = Math.max(poolSize, parseInt(argv[i + 1], 10));
    }
  }

  const mod = await createModule({
    locateFile: (p) => {
      if (p.endsWith('.wasm')) return wasmPath;
      if (p.endsWith('.worker.js')) return path.join(runtimePath, p);
      return p;
    },
    // Ensure worker script path is known for pthreads
    mainScriptUrlOrBlob: jsPath,
    // Pre-warm a worker pool if threads are requested
    pthreadPoolSize: poolSize > 1 ? poolSize : 0,
    print: (txt) => process.stdout.write(String(txt) + '\n'),
    printErr: (txt) => process.stderr.write(String(txt) + '\n'),
    onExit: (code) => { exitCode = code | 0; },
  });

  // With -sNODERAWFS=1, WASM syscalls access Node's real filesystem directly.
  // No mounting or path translation is required.
  //
  // Emscripten pthreads in Node sometimes require extra runtime setup.
  // Until we guarantee worker-based pthreads work in all environments,
  // coerce -P > 1 to -P 1 with a warning to avoid silent failures.
  const translatedArgs = (() => {
    const out = [];
    for (let i = 0; i < argv.length; i++) {
      const a = argv[i];
      out.push(a);
      if ((a === '-P' || a === '--threads') && i + 1 < argv.length) {
        const n = parseInt(argv[i + 1], 10);
        const allowThreads = process.env.LASTDB_ALLOW_THREADS === '1' || process.env.LASTDB_ALLOW_THREADS === 'true';
        if (Number.isFinite(n) && n > 1 && !allowThreads) {
          process.stderr.write('[lastdb-wasm] Note: forcing -P to 1 for Node runtime stability.\n');
          out.push('1');
          i++; // skip original value
          continue;
        }
      }
    }
    return out;
  })();

  try {
    mod.callMain(translatedArgs);
  } catch (e) {
    // Emscripten may throw on exit; prefer exitCode if set
    if (typeof e === 'number') {
      exitCode = e | 0;
    } else if (e && typeof e.status === 'number') {
      exitCode = e.status | 0;
    } else {
      console.error(e && e.stack ? e.stack : String(e));
      exitCode = exitCode || 1;
    }
  }

  process.exit(exitCode);
}

main().catch((err) => {
  console.error(err && err.stack ? err.stack : String(err));
  process.exit(1);
});
