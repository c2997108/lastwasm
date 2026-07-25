const fs = require('node:fs');

const [runtimePath, workerPath] = process.argv.slice(2);
if (!runtimePath || !workerPath) {
  throw new Error('usage: node patch-pthread-runtime.js <lastal.js> <lastal.worker.js>');
}

function replaceOnce(source, original, replacement, marker) {
  if (source.includes(marker)) return source;
  if (!source.includes(original)) throw new Error(`generated runtime hook not found: ${marker}`);
  return source.replace(original, replacement);
}

let runtime = fs.readFileSync(runtimePath, 'utf8');
runtime = replaceOnce(runtime,
`  wasmModule = module;
  if (!ENVIRONMENT_IS_PTHREAD) {
   removeRunDependency("wasm-instantiate");
  }`,
`  wasmModule = module;
  if (!ENVIRONMENT_IS_PTHREAD) {
   var pthreadPoolSize = Module["pthreadPoolSize"] || 0;
   if (pthreadPoolSize > 0) {
    addRunDependency("loading-workers");
    var workersToLoad = pthreadPoolSize;
    var workerLoaded = function() {
     if (--workersToLoad === 0) removeRunDependency("loading-workers");
    };
    for (var i = 0; i < pthreadPoolSize; ++i) {
     PThread.allocateUnusedWorker();
     PThread.loadWasmModuleToWorker(PThread.unusedWorkers[PThread.unusedWorkers.length - 1], workerLoaded);
    }
   }
   removeRunDependency("wasm-instantiate");
  }`,
'var pthreadPoolSize = Module["pthreadPoolSize"]');

runtime = replaceOnce(runtime,
`  worker.onerror = (e => {
   var message = "worker sent an error!";
   err(message + " " + e.filename + ":" + e.lineno + ": " + e.message);
   throw e;
  });`,
`  worker.onerror = (e => {
   var detail = "pthread worker failed at " + e.filename + ":" + e.lineno + ": " + e.message;
   err(detail);
   throw new Error(detail);
  });`,
'var detail = "pthread worker failed at "');

runtime = replaceOnce(runtime,
`   "wasmMemory": wasmMemory,
   "wasmModule": wasmModule
  });`,
`   "wasmMemory": wasmMemory,
   "wasmModule": wasmModule,
   "jslastProgressBuffer": Module["jslastProgressBuffer"]
  });`,
'"jslastProgressBuffer": Module["jslastProgressBuffer"]');

runtime = replaceOnce(runtime,
`function ___emscripten_thread_cleanup(thread) {
 if (!ENVIRONMENT_IS_PTHREAD) cleanupThread(thread); else postMessage({`,
`function ___emscripten_thread_cleanup(thread) {
 if (Module["jslastProgressBuffer"]) {
  Atomics.add(new Int32Array(Module["jslastProgressBuffer"]), 0, 1);
 }
 if (!ENVIRONMENT_IS_PTHREAD) cleanupThread(thread); else postMessage({`,
'Atomics.add(new Int32Array(Module["jslastProgressBuffer"]), 0, 1)');

runtime = replaceOnce(runtime,
` PThread.runningWorkers.push(worker);
 var pthread = PThread.pthreads[threadParams.pthread_ptr] = {`,
` PThread.runningWorkers.push(worker);
 Module["pthreadSpawnCount"] = (Module["pthreadSpawnCount"] || 0) + 1;
 Module["pthreadMaxRunning"] = Math.max(Module["pthreadMaxRunning"] || 0, PThread.runningWorkers.length);
 var pthread = PThread.pthreads[threadParams.pthread_ptr] = {`,
'Module["pthreadSpawnCount"] =');

fs.writeFileSync(runtimePath, runtime);

let worker = fs.readFileSync(workerPath, 'utf8');
worker = replaceOnce(worker,
`      Module['buffer'] = Module['wasmMemory'].buffer;

      Module['ENVIRONMENT_IS_PTHREAD'] = true;`,
`      Module['buffer'] = Module['wasmMemory'].buffer;

      Module['jslastProgressBuffer'] = e.data.jslastProgressBuffer;

      Module['ENVIRONMENT_IS_PTHREAD'] = true;`,
"Module['jslastProgressBuffer'] = e.data.jslastProgressBuffer");
fs.writeFileSync(workerPath, worker);
