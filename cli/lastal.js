#!/usr/bin/env node
/*
 Node.js CLI wrapper for LAST's lastal compiled to WebAssembly.

 Usage: cli/lastal.js [lastal args]
*/
const path = require('path');
const fs = require('fs');

async function main() {
  const argv = process.argv.slice(2);
  if (argv.length === 0 || argv.includes('-h') || argv.includes('--help')) {
    argv.length = 0;
    argv.push('-h');
  }

  const runtimePath = path.resolve(__dirname, '..', 'node');
  const wasmPath = path.join(runtimePath, 'lastal.wasm');
  const jsPath = path.join(runtimePath, 'lastal.js');

  if (!fs.existsSync(jsPath)) {
    console.error('error: Node-targeted lastal.js not found.');
    console.error('hint: run: bash node/build_node_lastal.sh');
    process.exit(2);
  }

  const createModule = require(jsPath);
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
    // Ensure worker script path is known and prewarm worker pool
    mainScriptUrlOrBlob: jsPath,
    pthreadPoolSize: poolSize > 1 ? poolSize : 0,
    print: (txt) => process.stdout.write(String(txt) + '\n'),
    printErr: (txt) => process.stderr.write(String(txt) + '\n'),
    onExit: (code) => { exitCode = code | 0; },
  });

  // Prepare in-memory filesystem and copy inputs from host FS.
  const { FS } = mod;
  try { if (!FS.analyzePath('/work').exists) FS.mkdir('/work'); } catch {}
  FS.chdir('/work');

  // Identify database name (first non-option arg) and query files (following non-option args)
  let dbArgIndex = -1;
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a === '--') { dbArgIndex = i + 1; break; }
    if (!a.startsWith('-')) { dbArgIndex = i; break; }
    // skip option parameters for options that take a value (simple heuristic)
    const opt = a.replace(/^--?/, '');
    const takesValue = new Set(['f','p','X','a','b','A','B','F','x','y','z','d','e','m','l','L','k','W','D','H','E','r','q','Q','S','U','w','i','b','B','C','s','u','P']);
    if (a.length === 2 && takesValue.has(opt) && i + 1 < argv.length && !argv[i+1].startsWith('-')) {
      i++;
    }
  }
  if (dbArgIndex < 0 || dbArgIndex >= argv.length) {
    // no db specified; just run as-is (likely -h)
    try { mod.callMain(argv); } catch (e) {}
    process.exit(exitCode);
  }
  const dbPathInput = argv[dbArgIndex];
  const dbDir = path.resolve(path.dirname(dbPathInput));
  const dbBase = path.basename(dbPathInput);

  // Copy db.* files
  for (const name of fs.readdirSync(dbDir)) {
    if (!name.startsWith(dbBase + '.')) continue;
    const buf = fs.readFileSync(path.join(dbDir, name));
    FS.writeFile('/work/' + name, buf);
  }

  // Collect query files (remaining non-option args after db name)
  const qArgs = [];
  for (let i = dbArgIndex + 1; i < argv.length; i++) {
    const a = argv[i];
    if (a.startsWith('-')) break; // stop at next option
    qArgs.push(a);
  }
  const qBasenames = [];
  for (const qa of qArgs) {
    const p = path.resolve(qa);
    const b = path.basename(p);
    FS.writeFile('/work/' + b, fs.readFileSync(p));
    qBasenames.push(b);
  }

  // Translate args: replace db path with dbBase, queries with basenames
  const translatedArgs = argv.slice(0, dbArgIndex)
    .concat([dbBase])
    .concat(qBasenames)
    .concat(argv.slice(dbArgIndex + 1 + qArgs.length).map(String));

  try {
    mod.callMain(translatedArgs);
  } catch (e) {
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
