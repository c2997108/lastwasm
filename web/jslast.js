import { JSLAST_VERSION } from './version.js?v=20260724-wasm-pthreads';
import createLastdbModule from './lastdb.js?v=20260724-wasm-pthreads';
import createLastalModule from './lastal.js?v=20260724-wasm-pthreads';

export { JSLAST_VERSION };

const DEFAULT_LASTDB_ARGS = ['--bits=4', '-R00', '-uNEAR', '-w1', '-W1', '-S1', '-C1', '-v'];
const DEFAULT_QUERY_BATCH_SIZE = '64M';
const LASTDB_MEMORY = 134217728;
const LASTAL_MEMORY = 268435456;

export async function runJsLast({
  refText,
  qryText,
  lastdbArgs = '',
  lastalArgs = '',
  requestedThreads = null,
  useDefaultDbArgs = true,
  onProgress = null,
} = {}) {
  const runStartedAt = performance.now();
  const progress = (event) => {
    if (onProgress) onProgress(event);
  };

  progress({ stage: 'parse', message: 'Reading FASTA...' });
  const refs = parseFasta(refText, 'ref');
  const queryRecords = parseFastaRecords(qryText, 'query');
  const queries = queryRecords.map(record => ({ name: record.name, seq: record.seq.toUpperCase() }));
  if (refs.length === 0) throw new Error('Reference FASTA has no sequences.');
  if (queries.length === 0) throw new Error('Query FASTA has no sequences.');
  const refLetters = totalSequenceLetters(refs);
  const queryLetters = totalSequenceLetters(queries);

  const requested = Math.max(1, Number(requestedThreads || 1));
  const threads = Math.max(1, Math.min(requested, queryRecords.length));
  const canUseWasmThreads = hasWasmThreadSupport();
  const runtime = canUseWasmThreads
    ? 'wasm-pthreads'
    : (threads > 1 ? 'asmjs-worker-pool' : 'asmjs');
  const warnings = [];
  if (requested > queryRecords.length) {
    warnings.push(`Search threads were limited to ${queryRecords.length} query FASTA record(s). LAST parallelism is split by query record.`);
  }
  if (!canUseWasmThreads && threads > 1) {
    warnings.push('Shared-memory WASM is unavailable; using the higher-memory asm.js worker-pool fallback.');
  }

  progress({
    stage: 'configured',
    message: `Reference ${refs.length} / query ${queries.length} / LAST ${runtime} runtime`,
    refs,
    queries,
    runtime,
    threads,
    warnings,
    metrics: {
      refLetters,
      queryLetters,
      queryRecords: queryRecords.length,
    },
  });

  const dbName = 'refdb';
  const dbStdout = [];
  const dbStderr = [];
  const dbArgs = normalizeDbArgs(splitArgs(lastdbArgs), useDefaultDbArgs);

  const dbStartedAt = performance.now();
  progress({ stage: 'lastdb', status: 'start', message: 'Running lastdb...' });
  const lastdbModule = await createLastdbModule(moduleOptions({
    moduleName: 'lastdb',
    memory: LASTDB_MEMORY,
    runtime: 'asmjs',
    threads: 1,
    stdout: line => dbStdout.push(line),
    stderr: line => {
      dbStderr.push(line);
      progress({ stage: 'lastdb-log', status: 'running', message: line });
    },
  }));

  const dbFS = lastdbModule.FS;
  ensureWorkdir(dbFS);
  dbFS.writeFile('ref.fa', refText || '');
  runMain(lastdbModule, ['-P', '1', ...dbArgs, dbName, 'ref.fa'], 'lastdb', dbStderr);
  const dbElapsedMs = performance.now() - dbStartedAt;
  progress({
    stage: 'lastdb',
    status: 'complete',
    message: 'lastdb complete',
    elapsedMs: dbElapsedMs,
  });
  const dbThreadStats = threadStats(lastdbModule);

  const dbFiles = dbFS.readdir('.')
    .filter(name => name !== '.' && name !== '..' && (name === dbName || name.startsWith(`${dbName}.`)));
  if (dbFiles.length === 0) throw new Error('lastdb did not create database files.');
  const dbFilePayloads = dbFiles.map(name => ({
    name,
    data: dbFS.readFile(name),
  }));

  const alStderr = [];
  const alArgs = ensureQueryBatchArgs(normalizeAlArgs(splitArgs(lastalArgs)));

  const searchStartedAt = performance.now();
  progress({
    stage: 'lastal',
    status: 'start',
    message: runtime === 'wasm-pthreads'
      ? `Running lastal with ${threads} shared-memory thread(s)...`
      : (threads > 1 ? `Running lastal on ${threads} search workers...` : 'Running lastal...'),
    workers: runtime.endsWith('worker-pool') ? threads : 0,
    threads,
    queryLetters,
  });

  let tabText = '';
  let alThreadStats = { spawned: 0, maxRunning: 0 };
  let searchWorkerStats = { started: 0, completed: 0 };
  let searchProgressStats = null;
  if (runtime.endsWith('worker-pool')) {
    const chunks = splitRecords(queryRecords, threads);
    const workerResult = await runLastalWorkerPool({
      chunks,
      dbFiles: dbFilePayloads,
      dbName,
      alArgs,
      useAsmJs: true,
      onProgress: progress,
    });
    tabText = workerResult.tabText;
    alStderr.push(...workerResult.alLog);
    searchWorkerStats = workerResult.workerStats;
    searchProgressStats = workerResult.progressStats;
  } else {
    const tabLines = [];
    const lastalModule = await createLastalModule(moduleOptions({
      moduleName: 'lastal',
      memory: LASTAL_MEMORY,
      runtime,
      threads,
      stdout: line => tabLines.push(line),
      stderr: line => alStderr.push(line),
    }));

    const alFS = lastalModule.FS;
    ensureWorkdir(alFS);
    for (const file of dbFilePayloads) {
      alFS.writeFile(file.name, file.data);
    }
    alFS.writeFile('query.fa', qryText || '');
    runMain(lastalModule, ['-f', 'TAB', '-P', String(threads), ...alArgs, dbName, 'query.fa'], 'lastal', alStderr);
    alThreadStats = threadStats(lastalModule);
    tabText = tabLines.length > 0 ? `${tabLines.join('\n')}\n` : '';
  }
  const searchElapsedMs = performance.now() - searchStartedAt;
  progress({
    stage: 'lastal',
    status: 'complete',
    message: 'lastal complete',
    elapsedMs: searchElapsedMs,
    workers: runtime.endsWith('worker-pool') ? threads : 0,
    threads,
  });

  const alignments = parseTabAlignments(tabText);
  const totalElapsedMs = performance.now() - runStartedAt;
  progress({ stage: 'done', message: `Done: ${alignments.length} alignments`, elapsedMs: totalElapsedMs });

  return {
    tabText,
    alignments,
    refs,
    queries,
    threads,
    runtime,
    warnings,
    dbThreadStats,
    alThreadStats,
    searchWorkerStats,
    searchProgressStats,
    timing: {
      totalMs: totalElapsedMs,
      dbMs: dbElapsedMs,
      searchMs: searchElapsedMs,
      refLetters,
      queryLetters,
      searchWorkers: threads,
    },
    dbLog: dbStderr.concat(dbStdout),
    alLog: alStderr,
  };
}

function threadStats(module) {
  return {
    spawned: Number(module.pthreadSpawnCount || 0),
    maxRunning: Number(module.pthreadMaxRunning || 0),
  };
}

function hasWasmThreadSupport() {
  return typeof WebAssembly === 'object'
    && typeof SharedArrayBuffer !== 'undefined'
    && typeof crossOriginIsolated !== 'undefined'
    && crossOriginIsolated;
}

function moduleOptions({ moduleName, memory, runtime, threads, stdout, stderr }) {
  const useAsmJs = runtime !== 'wasm-pthreads';
  return {
    noInitialRun: true,
    useAsmJs,
    INITIAL_MEMORY: memory,
    locateFile: path => new URL(path, import.meta.url).href,
    mainScriptUrlOrBlob: new URL(`./${moduleName}.js?v=${JSLAST_VERSION}`, import.meta.url).href,
    pthreadPoolSize: useAsmJs ? 0 : Math.max(0, threads - 1),
    print: stdout,
    printErr: stderr,
  };
}

async function runLastalWorkerPool({ chunks, dbFiles, dbName, alArgs, useAsmJs, onProgress }) {
  const workerUrl = new URL(`./lastal-worker.js?v=${JSLAST_VERSION}`, import.meta.url);
  const startedAt = performance.now();
  const progressState = {
    completed: 0,
    completedLetters: 0,
    completedWorkerMs: 0,
    totalLetters: chunks.reduce((sum, chunk) => sum + chunk.letters, 0),
    chunkLetters: chunks.map(chunk => chunk.letters),
    done: new Array(chunks.length).fill(false),
  };
  onProgress?.({
    stage: 'lastal-progress',
    status: 'start',
    message: `Search workers started: 0/${chunks.length}`,
    completed: 0,
    total: chunks.length,
    completedLetters: 0,
    totalLetters: progressState.totalLetters,
    searchElapsedMs: 0,
    estimatedRemainingMs: null,
  });

  const workerResults = await Promise.all(chunks.map((chunk, index) => runLastalWorker({
    workerUrl,
    index,
    queryText: recordsToFasta(chunk.records),
    queryLetters: chunk.letters,
    dbFiles,
    dbName,
    alArgs,
    useAsmJs,
    onProgress: event => {
      if (event.status === 'complete') {
        progressState.completed += 1;
        progressState.completedLetters += event.queryLetters || 0;
        progressState.completedWorkerMs += event.elapsedMs || 0;
        progressState.done[event.workerIndex] = true;
      }
      onProgress?.({
        ...event,
        completed: progressState.completed,
        total: chunks.length,
        completedLetters: progressState.completedLetters,
        totalLetters: progressState.totalLetters,
        searchElapsedMs: performance.now() - startedAt,
        estimatedRemainingMs: estimateSearchRemaining(progressState, performance.now() - startedAt),
      });
    },
  })));
  workerResults.sort((a, b) => a.index - b.index);
  return {
    tabText: mergeTabOutputs(workerResults.map(result => result.tabText)),
    alLog: workerResults.flatMap(result => result.alLog || []),
    workerStats: {
      started: chunks.length,
      completed: workerResults.length,
    },
    progressStats: {
      completedLetters: progressState.completedLetters,
      totalLetters: progressState.totalLetters,
      searchElapsedMs: performance.now() - startedAt,
    },
  };
}

function runLastalWorker({ workerUrl, index, queryText, queryLetters, dbFiles, dbName, alArgs, useAsmJs, onProgress }) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(workerUrl, {
      type: 'module',
      name: `lastal-search-${index + 1}`,
    });

    worker.onmessage = event => {
      const data = event.data || {};
      if (data.type === 'progress') {
        onProgress?.({
          stage: 'lastal-worker',
          status: 'start',
          message: `Running lastal worker ${data.index + 1}...`,
          workerIndex: data.index,
          queryLetters,
        });
        return;
      }
      worker.terminate();
      if (data.ok) {
        onProgress?.({
          stage: 'lastal-worker',
          status: 'complete',
          message: `lastal worker ${data.index + 1} complete`,
          workerIndex: data.index,
          queryLetters,
          elapsedMs: data.elapsedMs,
        });
        resolve(data);
      } else {
        reject(new Error(data.error || `lastal worker ${index + 1} failed`));
      }
    };

    worker.onerror = error => {
      worker.terminate();
      reject(new Error(error.message || `lastal worker ${index + 1} failed`));
    };

    worker.postMessage({
      index,
      queryText,
      queryLetters,
      dbFiles,
      dbName,
      alArgs,
      useAsmJs,
      memory: LASTAL_MEMORY,
    });
  });
}

function estimateSearchRemaining(state, searchElapsedMs) {
  if (state.completedLetters <= 0 || state.completedWorkerMs <= 0) return null;
  const msPerLetter = state.completedWorkerMs / state.completedLetters;
  let longestRemainingMs = 0;
  for (let i = 0; i < state.chunkLetters.length; i++) {
    if (state.done[i]) continue;
    const estimatedChunkMs = state.chunkLetters[i] * msPerLetter;
    longestRemainingMs = Math.max(longestRemainingMs, Math.max(0, estimatedChunkMs - searchElapsedMs));
  }
  return longestRemainingMs;
}

function splitRecords(records, requestedParts) {
  const parts = Math.max(1, Math.min(requestedParts, records.length));
  if (parts === 1) return [{ records: records.slice(), letters: totalRecordLetters(records) }];

  const chunks = [];
  const totalLetters = Math.max(1, totalRecordLetters(records));
  let start = 0;
  let currentLetters = 0;
  for (let i = 0; i < records.length; i++) {
    const remainingRecords = records.length - i;
    const remainingChunks = parts - chunks.length - 1;
    currentLetters += records[i].seq.length;
    const target = totalLetters * (chunks.length + 1) / parts;
    if (remainingChunks > 0 && currentLetters >= target && remainingRecords > remainingChunks) {
      const chunkRecords = records.slice(start, i + 1);
      chunks.push({ records: chunkRecords, letters: totalRecordLetters(chunkRecords) });
      start = i + 1;
    }
  }
  const finalRecords = records.slice(start);
  if (finalRecords.length > 0) chunks.push({ records: finalRecords, letters: totalRecordLetters(finalRecords) });
  return chunks;
}

function totalRecordLetters(records) {
  return records.reduce((sum, record) => sum + record.seq.length, 0);
}

function totalSequenceLetters(records) {
  return records.reduce((sum, record) => sum + record.seq.length, 0);
}

function recordsToFasta(records) {
  return records.map(record => record.text).join('');
}

function mergeTabOutputs(outputs) {
  const header = [];
  const alignments = [];
  let querySequences = 0;
  let queryLetters = 0;

  for (let outputIndex = 0; outputIndex < outputs.length; outputIndex++) {
    let inHeader = true;
    for (const line of String(outputs[outputIndex] || '').split(/\r?\n/)) {
      if (!line) continue;
      const summary = line.match(/^# Query sequences=(\d+) normal letters=(\d+)/);
      if (summary) {
        querySequences += Number(summary[1]);
        queryLetters += Number(summary[2]);
        continue;
      }
      if (line.startsWith('#') && inHeader) {
        if (outputIndex === 0) header.push(line);
        continue;
      }
      inHeader = false;
      if (!line.startsWith('#')) alignments.push(line);
    }
  }

  const lines = header.concat(alignments);
  if (querySequences > 0) {
    lines.push(`# Query sequences=${querySequences} normal letters=${queryLetters}`);
  }
  return lines.length > 0 ? `${lines.join('\n')}\n` : '';
}

function ensureWorkdir(fs) {
  if (!fs.analyzePath('/work').exists) fs.mkdir('/work');
  fs.chdir('/work');
}

function runMain(module, argv, program, stderrLines) {
  try {
    module.callMain(argv);
  } catch (error) {
    const status = typeof error?.status === 'number' ? error.status : undefined;
    if (error?.name === 'ExitStatus' && status === 0) return;
    const rawError = error?.stack || error?.message || String(error);
    const detail = stderrLines.length > 0 ? `\n${stderrLines.join('\n')}` : `\n${rawError}`;
    const exitInfo = typeof status === 'number' ? ` exited with code ${status}` : ' failed';
    throw new Error(`${program}${exitInfo}.${detail}`);
  }
}

function normalizeDbArgs(args, useDefaults) {
  const normalized = stripThreadArgs(args);
  return normalized.length > 0 || !useDefaults ? normalized : DEFAULT_LASTDB_ARGS.slice();
}

function normalizeAlArgs(args) {
  const normalized = [];
  for (let i = 0; i < args.length; i++) {
    const arg = args[i];
    const upper = arg.toUpperCase();
    if (arg === '-f' && i + 1 < args.length) {
      i++;
      continue;
    }
    if (upper.startsWith('-F')) continue;
    if (arg === '-P' && i + 1 < args.length) {
      i++;
      continue;
    }
    if (/^-P\d+$/.test(arg)) continue;
    normalized.push(arg);
  }
  return normalized;
}

function ensureQueryBatchArgs(args) {
  if (args.some(arg => arg === '-i' || /^-i.+/.test(arg))) return args;
  return ['-i', DEFAULT_QUERY_BATCH_SIZE, ...args];
}

function stripThreadArgs(args) {
  const normalized = [];
  for (let i = 0; i < args.length; i++) {
    const arg = args[i];
    if (arg === '-P' && i + 1 < args.length) {
      i++;
      continue;
    }
    if (/^-P\d+$/.test(arg)) continue;
    normalized.push(arg);
  }
  return normalized;
}

function splitArgs(text) {
  const args = [];
  const input = String(text || '');
  let current = '';
  let quote = '';
  let escaped = false;

  for (const char of input) {
    if (escaped) {
      current += char;
      escaped = false;
      continue;
    }
    if (char === '\\') {
      escaped = true;
      continue;
    }
    if (quote) {
      if (char === quote) quote = '';
      else current += char;
      continue;
    }
    if (char === '"' || char === "'") {
      quote = char;
      continue;
    }
    if (/\s/.test(char)) {
      if (current) {
        args.push(current);
        current = '';
      }
      continue;
    }
    current += char;
  }

  if (escaped) current += '\\';
  if (current) args.push(current);
  return args;
}

function parseFasta(text, fallbackPrefix) {
  const sequences = [];
  let name = '';
  let seq = [];

  const flush = () => {
    if (!name && seq.length === 0) return;
    const rawName = name || `${fallbackPrefix}_${sequences.length + 1}`;
    sequences.push({
      name: rawName.trim().split(/\s+/)[0] || `${fallbackPrefix}_${sequences.length + 1}`,
      seq: seq.join('').replace(/\s+/g, '').toUpperCase(),
    });
  };

  for (const rawLine of String(text || '').split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      flush();
      name = line.slice(1).trim();
      seq = [];
    } else {
      seq.push(line);
    }
  }
  flush();
  return sequences;
}

function parseFastaRecords(text, fallbackPrefix) {
  const records = [];
  let header = '';
  let seqLines = [];

  const flush = () => {
    if (!header && seqLines.length === 0) return;
    const fallbackName = `${fallbackPrefix}_${records.length + 1}`;
    const rawHeader = header || fallbackName;
    const name = rawHeader.trim().split(/\s+/)[0] || fallbackName;
    const seq = seqLines.join('').replace(/\s+/g, '');
    records.push({
      name,
      seq,
      text: `>${rawHeader}\n${seqLines.join('\n')}\n`,
    });
  };

  for (const rawLine of String(text || '').split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      flush();
      header = line.slice(1).trim();
      seqLines = [];
    } else {
      seqLines.push(line);
    }
  }

  flush();
  return records;
}

function parseTabAlignments(tabText) {
  const alignments = [];
  for (const line of String(tabText || '').split(/\r?\n/)) {
    if (!/^\d+\t/.test(line)) continue;
    const cols = line.split(/\t/);
    if (cols.length < 12) continue;
    alignments.push({
      score: Number(cols[0]),
      refName: cols[1],
      refStart: Number(cols[2]),
      refLen: Number(cols[3]),
      refStrand: cols[4],
      refTotal: Number(cols[5]),
      qryName: cols[6],
      qryStart: Number(cols[7]),
      qryLen: Number(cols[8]),
      qryStrand: cols[9],
      qryTotal: Number(cols[10]),
      blocks: cols[11],
    });
  }
  return alignments;
}
