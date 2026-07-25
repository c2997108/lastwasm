import { JSLAST_VERSION } from './version.js?v=20260725-axis-labels';
import { renderReglDotPlot } from './regl-dotplot.js?v=20260725-axis-labels';

window.__JSLAST_MODULE_READY = true;

const $ = (sel) => document.querySelector(sel);
const statusEl = $('#status');
const logEl = $('#log');
const tabEl = $('#tab');
const tabSummaryEl = $('#tabSummary');
const downloadTabBtn = $('#downloadTabBtn');
const timeEl = $('#timeEstimate');
const runBtn = $('#runBtn');
const workerInput = $('#workerCount');
const seedHitLimitInput = $('#seedHitLimit');
const queryBatchSizeInput = $('#queryBatchSize');
const plotLengthCutoffInput = $('#plotLengthCutoff');
const useFastaOrderInput = $('#useFastaOrder');
const TIMING_MODEL_KEY = `lastjsTimingModel:${JSLAST_VERSION}`;
const THREAD_PROGRESS_COMPLETED = 0;
const THREAD_PROGRESS_STARTED_BATCHES = 1;
const THREAD_PROGRESS_THREADS = 2;
const THREAD_PROGRESS_TOTAL_BATCHES = 3;
const THREAD_PROGRESS_DONE = 4;
const MAX_PLOT_WORKERS = 8;
const MAX_UNSHARED_PLOT_BYTES = 32 * 1024 * 1024;

let runClock = null;
let estimateState = null;
let threadProgressView = null;
let threadProgressTimer = null;
let latestTabBuffer = null;
let activeDotPlot = null;

function setStatus(message) {
  if (statusEl) statusEl.textContent = message;
}

function appendLog(message) {
  if (!logEl) return;
  logEl.textContent += `${message}\n`;
  logEl.scrollTop = logEl.scrollHeight;
}

async function readFileInput(input) {
  const file = input?.files?.[0];
  if (!file) return null;
  return file.text();
}

function configureDefaults() {
  const verEl = $('#appVer');
  if (verEl) verEl.textContent = `(v${JSLAST_VERSION})`;

  if (workerInput && !workerInput.value) {
    const hardware = navigator.hardwareConcurrency || 1;
    workerInput.value = String(Math.max(1, Math.min(hardware, 8)));
    workerInput.max = String(Math.max(1, hardware || 64));
  }

  if (seedHitLimitInput && !seedHitLimitInput.value) seedHitLimitInput.value = '100';
  if (queryBatchSizeInput && !queryBatchSizeInput.value) queryBatchSizeInput.value = '64M';
}

async function run() {
  if (runBtn) runBtn.disabled = true;
  if (logEl) logEl.textContent = '';
  if (tabEl) tabEl.textContent = '';
  if (timeEl) timeEl.textContent = '';
  if (tabSummaryEl) tabSummaryEl.textContent = '';
  latestTabBuffer = null;
  if (downloadTabBtn) downloadTabBtn.disabled = true;
  clearPlot();

  try {
    setStatus('Preparing...');
    const refText = await readFileInput($('#refFasta'));
    const qryText = await readFileInput($('#qryFasta'));
    if (!refText || !qryText) {
      setStatus('Select both reference and query FASTA files.');
      return;
    }

    const requestedThreads = Math.max(1, Number(workerInput?.value || navigator.hardwareConcurrency || 1));
    const seedHitLimit = parseNonNegativeInteger(seedHitLimitInput?.value, 100);
    if (seedHitLimitInput) seedHitLimitInput.value = String(seedHitLimit);
    const queryBatchSize = normalizeQueryBatchSize(queryBatchSizeInput?.value);
    if (queryBatchSizeInput) queryBatchSizeInput.value = queryBatchSize;
    const lastalArgs = `-m${seedHitLimit} -i${queryBatchSize}`;
    const plotLengthCutoff = parsePlotLengthCutoff(plotLengthCutoffInput?.value);
    if (plotLengthCutoffInput) plotLengthCutoffInput.value = String(plotLengthCutoff);
    const useFastaOrder = Boolean(useFastaOrderInput?.checked);
    appendLog(`$ lastdb + lastal ${lastalArgs}`);
    appendLog(`dotplot hide-length<=${plotLengthCutoff}bp`);
    appendLog(`dotplot sequence-order=${useFastaOrder ? 'fasta' : 'tab'}`);
    startEstimateClock({
      refText,
      qryText,
      requestedThreads,
    });

    const result = await runJsLastInWorker({
      refText,
      qryText,
      lastdbArgs: '',
      lastalArgs,
      requestedThreads,
      useDefaultDbArgs: true,
      onProgress: (event) => {
        updateEstimateFromProgress(event);
        if (event.message) setStatus(event.message);
        if (event.stage === 'configured') {
          appendLog(`reference=${event.refs.length}`);
          appendLog(`query=${event.queries.length}`);
          appendLog(`runtime=compiled LAST ${event.runtime} threads=${event.threads}`);
          for (const warning of event.warnings || []) appendLog(`[WARN] ${warning}`);
        }
      },
      onTabPreview: showTabPreview,
    });

    latestTabBuffer = result.tabBuffer;
    if (downloadTabBtn) downloadTabBtn.disabled = !latestTabBuffer;
    if (tabSummaryEl) tabSummaryEl.textContent = formatTabSummary(result.tabByteLength, true);
    updateTimingModel(result.timing);
    appendLog(formatResultPipelineProfile(result));
    if (result.runtime.endsWith('worker-pool')) {
      appendLog(`lastal-workers started=${result.searchWorkerStats.started} completed=${result.searchWorkerStats.completed}`);
    } else if (result.runtime === 'wasm-pthreads') {
      appendLog(`lastal-pthreads spawned=${result.alThreadStats.spawned} maxRunning=${result.alThreadStats.maxRunning}`);
    }
    appendLog(`threads=${result.threads}, alignments=${result.alignmentCount}`);
    setStatus(`Result ready: ${result.alignmentCount} alignments`);
    finishEstimateClock(result.timing);
    await afterNextPaint();
    setStatus('Rendering dot plot...');
    await afterNextPaint();
    const plotStats = await renderReglFromTab(
      result.tabBuffer,
      requestedThreads,
      plotLengthCutoff,
      useFastaOrder ? { refs: result.refs, queries: result.queries } : null,
    );
    for (const line of formatPlotProfile(plotStats)) appendLog(line);
    if (plotStats.workers > 0) {
      appendLog(`dotplot-parser workers=${plotStats.workers} alignments=${plotStats.alignments} rendered=${plotStats.rendered} hide-length<=${plotStats.lengthCutoff}bp`);
    }
    setStatus('Done');
  } catch (error) {
    setStatus('Error');
    appendLog(error?.stack || String(error));
  } finally {
    stopEstimateClock();
    if (runBtn) runBtn.disabled = false;
  }
}

function runJsLastInWorker(payload) {
  return new Promise((resolve, reject) => {
    let transferStartedAt = null;
    const worker = new Worker(new URL(`./jslast-runner-worker.js?v=${JSLAST_VERSION}`, import.meta.url), {
      type: 'module',
      name: 'jslast-runner',
    });

    worker.onmessage = event => {
      const data = event.data || {};
      if (data.type === 'progress') {
        if (data.event?.stage === 'tab-encode' && data.event.status === 'complete') {
          transferStartedAt = performance.now();
        }
        payload.onProgress?.(data.event);
        return;
      }
      if (data.type === 'tab-preview') {
        payload.onTabPreview?.(data);
        afterNextPaint().then(() => worker.postMessage({ type: 'tab-preview-painted' }));
        return;
      }
      worker.terminate();
      stopThreadProgress();
      if (data.type === 'done') {
        data.result.timing.tabTransferMs = transferStartedAt === null ? 0 : performance.now() - transferStartedAt;
        resolve({ ...data.result, tabBuffer: data.tabBuffer });
      } else {
        reject(new Error(data.error || 'LAST worker failed'));
      }
    };

    worker.onerror = error => {
      worker.terminate();
      stopThreadProgress();
      reject(new Error(error.message || 'LAST worker failed'));
    };

    const { onProgress, onTabPreview, ...runPayload } = payload;
    if (hasSharedMemory()) {
      const threadProgressBuffer = new SharedArrayBuffer(Int32Array.BYTES_PER_ELEMENT * 5);
      runPayload.threadProgressBuffer = threadProgressBuffer;
      startThreadProgress(threadProgressBuffer);
    }
    worker.postMessage({ type: 'run', payload: runPayload });
  });
}

function showTabPreview({ text, truncated, totalChars }) {
  if (tabEl) {
    tabEl.textContent = truncated
      ? `${text}\n# Preview limited to the first 20 KiB. Download TAB contains the complete result.\n`
      : text;
  }
  if (tabSummaryEl) {
    tabSummaryEl.textContent = truncated
      ? `Preview ready · ${formatNumber(totalChars)} characters total`
      : 'Complete result ready';
  }
  setStatus(truncated ? 'Result preview ready; preparing full TAB...' : 'Result ready');
}

function formatTabSummary(bytes, downloadable) {
  if (!Number.isFinite(bytes)) return '';
  const size = bytes >= 1024 * 1024
    ? `${(bytes / (1024 * 1024)).toFixed(1)} MiB`
    : `${Math.ceil(bytes / 1024)} KiB`;
  return downloadable ? `Complete TAB · ${size}` : size;
}

function formatResultPipelineProfile(result) {
  const timing = result?.timing || {};
  const tabMiB = Number(result?.tabByteLength || 0) / (1024 * 1024);
  return [
    'result-timing',
    `lastdb=${formatMs(timing.dbMs)}`,
    `lastal-compute=${formatMs(timing.lastalComputeMs)}`,
    `tab-preview=${formatMs(timing.tabPreviewMs)}`,
    `tab-join=${formatMs(timing.tabFinalizeMs)}`,
    `utf8-encode=${formatMs(timing.tabEncodeMs)}`,
    `transfer=${formatMs(timing.tabTransferMs)}`,
    `tab=${tabMiB.toFixed(1)}MiB`,
  ].join(' ');
}

function formatNumber(value) {
  return Number(value || 0).toLocaleString();
}

function downloadLatestTab() {
  if (!latestTabBuffer) return;
  const url = URL.createObjectURL(new Blob([latestTabBuffer], { type: 'text/tab-separated-values' }));
  const link = document.createElement('a');
  link.href = url;
  link.download = 'lastal-result.tab';
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function startThreadProgress(buffer) {
  stopThreadProgress();
  threadProgressView = new Int32Array(buffer);
  threadProgressTimer = setInterval(updateThreadProgress, 50);
}

function stopThreadProgress() {
  if (threadProgressTimer) clearInterval(threadProgressTimer);
  threadProgressTimer = null;
  threadProgressView = null;
}

function updateThreadProgress() {
  const percent = readThreadProgressPercent();
  if (percent === null || !estimateState || estimateState.phase !== 'Searching') return;
  estimateState.searchPercent = percent;
  setStatus(`Searching... ${percent}%`);
  updateTimeEstimate();
}

function readThreadProgressPercent() {
  if (!threadProgressView) return null;
  if (Atomics.load(threadProgressView, THREAD_PROGRESS_DONE) > 0) return 100;
  const threads = Atomics.load(threadProgressView, THREAD_PROGRESS_THREADS);
  const totalBatches = Atomics.load(threadProgressView, THREAD_PROGRESS_TOTAL_BATCHES);
  if (threads < 2 || totalBatches < 1) return null;

  const startedBatches = Atomics.load(threadProgressView, THREAD_PROGRESS_STARTED_BATCHES);
  const completedChildren = Atomics.load(threadProgressView, THREAD_PROGRESS_COMPLETED);
  const completedBatches = Math.max(0, startedBatches - 1);
  const childrenPerBatch = threads - 1;
  const currentBatchChildren = Math.max(0, completedChildren - completedBatches * childrenPerBatch);
  const completedUnits = completedBatches * threads + Math.min(childrenPerBatch, currentBatchChildren);
  return Math.max(0, Math.min(99, Math.floor(completedUnits * 100 / (totalBatches * threads))));
}

function startEstimateClock({ refText, qryText, requestedThreads }) {
  stopEstimateClock();
  const summary = summarizeFastaWork(refText, qryText, requestedThreads);
  const model = loadTimingModel();
  estimateState = {
    startedAt: performance.now(),
    phase: 'Preparing',
    estimatedTotalMs: estimateTotalMs(summary, model),
    searchEstimatedRemainingMs: null,
    searchStartedAt: null,
    searchPercent: null,
    completedWorkers: 0,
    totalWorkers: 0,
  };
  updateTimeEstimate();
  runClock = setInterval(updateTimeEstimate, 1000);
}

function stopEstimateClock() {
  if (runClock) clearInterval(runClock);
  runClock = null;
}

function finishEstimateClock(timing) {
  if (!timeEl || !timing) return;
  timeEl.textContent = `Elapsed ${formatDuration(timing.totalMs)}`;
}

function updateEstimateFromProgress(event) {
  if (!estimateState) return;
  if (event.stage === 'lastdb') {
    estimateState.phase = event.status === 'complete' ? 'Index ready' : 'Building index';
  } else if (event.stage === 'lastal') {
    estimateState.phase = event.status === 'complete' ? 'Search complete' : 'Searching';
    if (event.status !== 'complete' && estimateState.searchStartedAt === null) {
      estimateState.searchStartedAt = performance.now();
    }
    if (event.status === 'complete') estimateState.searchEstimatedRemainingMs = 0;
    estimateState.totalWorkers = event.workers || estimateState.totalWorkers;
  } else if (event.stage === 'lastal-progress' || event.stage === 'lastal-worker') {
    estimateState.phase = 'Searching';
    estimateState.completedWorkers = event.completed ?? estimateState.completedWorkers;
    estimateState.totalWorkers = event.total ?? estimateState.totalWorkers;
    if (Number.isFinite(event.estimatedRemainingMs)) {
      estimateState.searchEstimatedRemainingMs = Math.max(0, event.estimatedRemainingMs);
    }
    if (Number.isFinite(event.percent)) estimateState.searchPercent = event.percent;
  } else if (event.stage === 'done') {
    estimateState.phase = 'Done';
  }
  updateTimeEstimate();
}

function updateTimeEstimate() {
  if (!timeEl || !estimateState) return;
  const elapsedMs = performance.now() - estimateState.startedAt;
  let remainingMs = null;
  if (Number.isFinite(estimateState.searchEstimatedRemainingMs)) {
    remainingMs = estimateState.searchEstimatedRemainingMs;
  } else if (estimateState.phase === 'Searching'
      && Number.isFinite(estimateState.searchPercent)
      && estimateState.searchPercent > 0
      && Number.isFinite(estimateState.searchStartedAt)) {
    const searchElapsedMs = performance.now() - estimateState.searchStartedAt;
    remainingMs = searchElapsedMs * (100 - estimateState.searchPercent) / estimateState.searchPercent;
  } else if (Number.isFinite(estimateState.estimatedTotalMs)) {
    remainingMs = Math.max(0, estimateState.estimatedTotalMs - elapsedMs);
  }

  const parts = [`Elapsed ${formatDuration(elapsedMs)}`];
  if (remainingMs !== null) {
    parts.push(`remaining ~${formatDuration(remainingMs)}`);
  } else {
    parts.push('estimating...');
  }
  if (estimateState.totalWorkers > 0) {
    parts.push(`${estimateState.completedWorkers}/${estimateState.totalWorkers} workers`);
  } else if (Number.isFinite(estimateState.searchPercent)) {
    parts.push(`${estimateState.searchPercent}%`);
  } else {
    parts.push(estimateState.phase);
  }
  timeEl.textContent = parts.join(' · ');
}

function summarizeFastaWork(refText, qryText, requestedThreads) {
  const ref = fastaSummary(refText);
  const qry = fastaSummary(qryText);
  return {
    refLetters: ref.letters,
    queryLetters: qry.letters,
    queryRecords: qry.records,
    workers: Math.max(1, Math.min(requestedThreads, Math.max(1, qry.records))),
  };
}

function fastaSummary(text) {
  let records = 0;
  let letters = 0;
  let hasSequence = false;
  for (const rawLine of String(text || '').split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      records += 1;
    } else {
      letters += line.replace(/\s+/g, '').length;
      hasSequence = true;
    }
  }
  return {
    records: records || (hasSequence ? 1 : 0),
    letters,
  };
}

function estimateTotalMs(summary, model) {
  if (!model) return null;
  const dbMs = model.dbMsPerRefLetter * summary.refLetters;
  const searchMs = model.searchWorkerMsPerQueryLetter * summary.queryLetters / summary.workers;
  const estimate = dbMs + searchMs;
  return Number.isFinite(estimate) && estimate > 0 ? estimate : null;
}

function loadTimingModel() {
  try {
    const raw = localStorage.getItem(TIMING_MODEL_KEY);
    if (!raw) return null;
    const model = JSON.parse(raw);
    if (!Number.isFinite(model.dbMsPerRefLetter) || !Number.isFinite(model.searchWorkerMsPerQueryLetter)) return null;
    return model;
  } catch {
    return null;
  }
}

function updateTimingModel(timing) {
  if (!timing || timing.refLetters <= 0 || timing.queryLetters <= 0) return;
  const next = {
    dbMsPerRefLetter: timing.dbMs / timing.refLetters,
    searchWorkerMsPerQueryLetter: timing.searchMs * Math.max(1, timing.searchWorkers) / timing.queryLetters,
  };
  if (!Number.isFinite(next.dbMsPerRefLetter) || !Number.isFinite(next.searchWorkerMsPerQueryLetter)) return;
  const previous = loadTimingModel();
  const alpha = previous ? 0.25 : 1;
  const model = previous ? {
    dbMsPerRefLetter: previous.dbMsPerRefLetter * (1 - alpha) + next.dbMsPerRefLetter * alpha,
    searchWorkerMsPerQueryLetter: previous.searchWorkerMsPerQueryLetter * (1 - alpha) + next.searchWorkerMsPerQueryLetter * alpha,
  } : next;
  try {
    localStorage.setItem(TIMING_MODEL_KEY, JSON.stringify(model));
  } catch {}
}

function formatDuration(ms) {
  const totalSeconds = Math.max(0, Math.round(ms / 1000));
  const hours = Math.floor(totalSeconds / 3600);
  const minutes = Math.floor((totalSeconds % 3600) / 60);
  const seconds = totalSeconds % 60;
  if (hours > 0) return `${hours}:${String(minutes).padStart(2, '0')}:${String(seconds).padStart(2, '0')}`;
  return `${minutes}:${String(seconds).padStart(2, '0')}`;
}

function afterNextPaint() {
  return new Promise(resolve => {
    let settled = false;
    const finish = () => {
      if (settled) return;
      settled = true;
      resolve();
    };
    setTimeout(finish, 250);
    requestAnimationFrame(() => requestAnimationFrame(finish));
  });
}

function clearPlot() {
  const plotDiv = $('#dotplot');
  if (!plotDiv) return;
  activeDotPlot?.destroy();
  activeDotPlot = null;
  plotDiv.textContent = '';
}

async function renderReglFromTab(tabBuffer, requestedWorkers, lengthCutoff, sequenceOrder) {
  const totalStartedAt = performance.now();
  const memoryStartMiB = usedHeapMiB();
  const plotDiv = $('#dotplot');
  if (!plotDiv || !tabBuffer) return { workers: 0, alignments: 0, rendered: 0 };
  if (typeof window.createREGL !== 'function') {
    plotDiv.textContent = 'The WebGL renderer is unavailable; TAB output is shown above.';
    return { workers: 0, alignments: 0, rendered: 0 };
  }
  const canShare = hasSharedMemory();
  if (!canShare && tabBuffer.byteLength > MAX_UNSHARED_PLOT_BYTES) {
    plotDiv.textContent = 'Dot plot skipped because shared memory is unavailable for this large result.';
    return { workers: 0, alignments: 0, rendered: 0 };
  }

  const plotData = await parsePlotDataInParallel(
    tabBuffer,
    requestedWorkers,
    canShare,
    lengthCutoff,
    sequenceOrder,
  );
  if (plotData.alignmentCount === 0) {
    plotDiv.textContent = 'No plottable alignments.';
    return { workers: plotData.workerCount, alignments: 0, rendered: 0, lengthCutoff };
  }
  if (plotData.renderedCount === 0) {
    plotDiv.textContent = `No alignments longer than ${lengthCutoff.toLocaleString()} bp.`;
    return {
      workers: plotData.workerCount,
      alignments: plotData.alignmentCount,
      rendered: 0,
      lengthCutoff,
    };
  }
  const longTaskProfile = startLongTaskProfile();
  const memoryBeforeRendererMiB = usedHeapMiB();
  const rendered = await renderReglDotPlot(plotDiv, plotData);
  activeDotPlot = rendered.renderer;
  const longTasks = longTaskProfile.stop();
  return {
    workers: plotData.workerCount,
    alignments: plotData.alignmentCount,
    rendered: plotData.renderedCount,
    lengthCutoff,
    timing: {
      ...plotData.timing,
      rendererUploadMs: rendered.timing.uploadMs,
      rendererDrawSubmitMs: rendered.timing.drawSubmitMs,
      rendererGpuMs: rendered.timing.gpuCompleteMs,
      rendererPaintMs: rendered.timing.firstPaintMs,
      rendererTotalMs: rendered.timing.totalMs,
      longTasks,
      totalMs: performance.now() - totalStartedAt,
    },
    memory: {
      startMiB: memoryStartMiB,
      ...plotData.memory,
      beforeRendererMiB: memoryBeforeRendererMiB,
      afterRendererMiB: usedHeapMiB(),
    },
  };
}

async function parsePlotDataInParallel(tabBuffer, requestedWorkers, canShare, lengthCutoff, sequenceOrder) {
  const totalStartedAt = performance.now();
  const workerLimit = canShare
    ? Math.max(1, Math.min(MAX_PLOT_WORKERS, Number(requestedWorkers) || 1))
    : 1;
  let buffer;
  const bufferStartedAt = performance.now();
  if (canShare) {
    buffer = new SharedArrayBuffer(tabBuffer.byteLength);
    new Uint8Array(buffer).set(new Uint8Array(tabBuffer));
  } else {
    buffer = tabBuffer.slice(0);
  }
  const bufferCopyMs = performance.now() - bufferStartedAt;
  const afterBufferMiB = usedHeapMiB();
  const rangesStartedAt = performance.now();
  const ranges = splitPlotRanges(new Uint8Array(buffer), workerLimit);
  const rangeSplitMs = performance.now() - rangesStartedAt;

  const workerUrl = new URL(`./plot-parser-worker.js?v=${JSLAST_VERSION}`, import.meta.url);
  const workersStartedAt = performance.now();
  const sessions = ranges.map((range, index) => ({
    index,
    range,
    worker: new Worker(workerUrl, { type: 'module', name: `jslast-plot-${index + 1}` }),
  }));
  const workerCreateMs = performance.now() - workersStartedAt;

  try {
    const parseStartedAt = performance.now();
    const metadata = await Promise.all(sessions.map((session, index) => requestPlotWorker(
      session.worker,
      {
        type: 'parse',
        index: session.index,
        buffer,
        start: session.range.start,
        end: session.range.end,
        lengthCutoff,
      },
      canShare || index > 0 ? [] : [buffer],
    )));
    const parsePassMs = performance.now() - parseStartedAt;
    const parseWorkers = summarizeWorkerTimings(metadata);
    const afterParseMiB = usedHeapMiB();
    const mergeStartedAt = performance.now();
    const merged = mergePlotMetadata(metadata, sequenceOrder);
    const metadataMergeMs = performance.now() - mergeStartedAt;
    const buildStartedAt = performance.now();
    const chunks = await Promise.all(sessions.map(session => requestPlotWorker(session.worker, {
      type: 'build',
      index: session.index,
      refOffsets: Array.from(merged.refOffsets),
      qryOffsets: Array.from(merged.qryOffsets),
      alignmentStart: merged.alignmentStartByIndex.get(session.index) || 0,
      lengthCutoff,
    })));
    const buildPassMs = performance.now() - buildStartedAt;
    const buildWorkers = summarizeWorkerTimings(chunks);
    const afterBuildMiB = usedHeapMiB();
    const finalizeStartedAt = performance.now();
    chunks.sort((left, right) => left.index - right.index);
    const result = {
      ...merged,
      chunks,
      workerCount: sessions.length,
      renderedCount: chunks.reduce((sum, chunk) => sum + (chunk.alignmentCount || 0), 0),
      timing: {
        bufferCopyMs,
        rangeSplitMs,
        workerCreateMs,
        parsePassMs,
        parseWorkers,
        parseOverheadMs: Math.max(0, parsePassMs - parseWorkers.maxTotalMs),
        metadataMergeMs,
        buildPassMs,
        buildWorkers,
        buildOverheadMs: Math.max(0, buildPassMs - buildWorkers.maxTotalMs),
        finalizeMs: performance.now() - finalizeStartedAt,
        parserTotalMs: performance.now() - totalStartedAt,
      },
      memory: {
        afterBufferMiB,
        afterParseMiB,
        afterBuildMiB,
      },
    };
    return result;
  } finally {
    for (const session of sessions) session.worker.terminate();
  }
}

function summarizeWorkerTimings(results) {
  const timings = results.map(result => result.timing || {});
  const max = key => timings.reduce((value, timing) => Math.max(value, Number(timing[key]) || 0), 0);
  const sum = key => timings.reduce((value, timing) => value + (Number(timing[key]) || 0), 0);
  return {
    maxTotalMs: max('totalMs'),
    sumTotalMs: sum('totalMs'),
    maxDecodeMs: max('decodeMs'),
    sumDecodeMs: sum('decodeMs'),
    maxScanMs: max('scanMs'),
    sumScanMs: sum('scanMs'),
    maxArraysMs: max('arraysMs'),
    sumArraysMs: sum('arraysMs'),
  };
}

function startLongTaskProfile() {
  const durations = [];
  let observer = null;
  if (typeof PerformanceObserver === 'function'
    && PerformanceObserver.supportedEntryTypes?.includes('longtask')) {
    observer = new PerformanceObserver(list => {
      for (const entry of list.getEntries()) durations.push(entry.duration);
    });
    observer.observe({ type: 'longtask', buffered: false });
  }
  return {
    stop() {
      if (observer) {
        for (const entry of observer.takeRecords()) durations.push(entry.duration);
        observer.disconnect();
      }
      return {
        count: durations.length,
        totalMs: durations.reduce((sum, duration) => sum + duration, 0),
        maxMs: durations.reduce((max, duration) => Math.max(max, duration), 0),
      };
    },
  };
}

function formatPlotProfile(plotStats) {
  const timing = plotStats?.timing;
  if (!timing) return [];
  const parse = timing.parseWorkers || {};
  const build = timing.buildWorkers || {};
  const lines = [
    `dotplot-timing total=${formatMs(timing.totalMs)} parser-total=${formatMs(timing.parserTotalMs)} renderer-total=${formatMs(timing.rendererTotalMs)}`,
    `dotplot-timing setup buffer-copy=${formatMs(timing.bufferCopyMs)} range-split=${formatMs(timing.rangeSplitMs)} worker-create=${formatMs(timing.workerCreateMs)} metadata-merge=${formatMs(timing.metadataMergeMs)} finalize=${formatMs(timing.finalizeMs)}`,
    `dotplot-timing parse-pass=${formatMs(timing.parsePassMs)} worker-max=${formatMs(parse.maxTotalMs)} decode-max=${formatMs(parse.maxDecodeMs)} scan-max=${formatMs(parse.maxScanMs)} startup-message=${formatMs(timing.parseOverheadMs)}`,
    `dotplot-timing build-pass=${formatMs(timing.buildPassMs)} worker-max=${formatMs(build.maxTotalMs)} decode-max=${formatMs(build.maxDecodeMs)} scan-max=${formatMs(build.maxScanMs)} arrays-max=${formatMs(build.maxArraysMs)} message-clone=${formatMs(timing.buildOverheadMs)}`,
    `dotplot-timing regl upload=${formatMs(timing.rendererUploadMs)} draw-submit=${formatMs(timing.rendererDrawSubmitMs)} gpu=${formatMs(timing.rendererGpuMs)} paint=${formatMs(timing.rendererPaintMs)}`,
    `dotplot-longtasks count=${timing.longTasks?.count || 0} total=${formatMs(timing.longTasks?.totalMs)} max=${formatMs(timing.longTasks?.maxMs)}`,
  ];
  const memory = plotStats.memory || {};
  if (Number.isFinite(memory.startMiB)) {
    lines.push(`dotplot-memory heap start=${formatMiB(memory.startMiB)} after-buffer=${formatMiB(memory.afterBufferMiB)} after-parse=${formatMiB(memory.afterParseMiB)} after-build=${formatMiB(memory.afterBuildMiB)} before-renderer=${formatMiB(memory.beforeRendererMiB)} after-renderer=${formatMiB(memory.afterRendererMiB)}`);
  }
  return lines;
}

function usedHeapMiB() {
  const bytes = performance.memory?.usedJSHeapSize;
  return Number.isFinite(bytes) ? bytes / (1024 * 1024) : null;
}

function formatMs(value) {
  return `${(Number(value) || 0).toFixed(1)}ms`;
}

function formatMiB(value) {
  return `${(Number(value) || 0).toFixed(1)}MiB`;
}

function splitPlotRanges(bytes, workerCount) {
  if (bytes.length === 0) return [{ start: 0, end: 0 }];
  const boundaries = [0];
  for (let index = 1; index < workerCount; index++) {
    let boundary = Math.floor(bytes.length * index / workerCount);
    while (boundary < bytes.length && bytes[boundary - 1] !== 10) boundary += 1;
    if (boundary > boundaries[boundaries.length - 1] && boundary < bytes.length) {
      boundaries.push(boundary);
    }
  }
  boundaries.push(bytes.length);
  return boundaries.slice(0, -1).map((start, index) => ({ start, end: boundaries[index + 1] }));
}

function requestPlotWorker(worker, message, transfer = []) {
  return new Promise((resolve, reject) => {
    const cleanup = () => {
      worker.removeEventListener('message', onMessage);
      worker.removeEventListener('error', onError);
    };
    const onMessage = event => {
      cleanup();
      if (event.data?.type === 'error') {
        reject(new Error(event.data.error || 'Dot-plot worker failed'));
      } else {
        resolve(event.data);
      }
    };
    const onError = error => {
      cleanup();
      reject(new Error(error.message || 'Dot-plot worker failed'));
    };
    worker.addEventListener('message', onMessage);
    worker.addEventListener('error', onError);
    worker.postMessage(message, transfer);
  });
}

function mergePlotMetadata(metadata, sequenceOrder) {
  const refOrders = [];
  const qryOrders = [];
  const refLenByName = new Map();
  const qryLenByName = new Map();
  mergeSequenceMetadata(sequenceOrder?.refs, refOrders, refLenByName);
  mergeSequenceMetadata(sequenceOrder?.queries, qryOrders, qryLenByName);
  const alignmentStartByIndex = new Map();
  let alignmentCount = 0;
  let renderedCount = 0;
  metadata.sort((left, right) => left.index - right.index);
  for (const chunk of metadata) {
    alignmentStartByIndex.set(chunk.index, renderedCount);
    alignmentCount += chunk.alignmentCount || 0;
    renderedCount += chunk.renderedCount || 0;
    mergeSequenceMetadata(chunk.refs, refOrders, refLenByName);
    mergeSequenceMetadata(chunk.queries, qryOrders, qryLenByName);
  }
  const refOffsets = offsetsFor(refOrders, refLenByName);
  const qryOffsets = offsetsFor(qryOrders, qryLenByName);
  return {
    alignmentCount,
    renderedCount,
    alignmentStartByIndex,
    refOrders,
    qryOrders,
    refLenByName,
    qryLenByName,
    refOffsets,
    qryOffsets,
    maxX: totalLength(refOrders, refLenByName),
    maxY: totalLength(qryOrders, qryLenByName),
  };
}

function parsePlotLengthCutoff(value) {
  return parseNonNegativeInteger(value, 1000);
}

function parseNonNegativeInteger(value, fallback) {
  if (String(value ?? '').trim() === '') return fallback;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? Math.max(0, Math.floor(parsed)) : fallback;
}

function normalizeQueryBatchSize(value) {
  const normalized = String(value ?? '').trim().toUpperCase();
  const match = normalized.match(/^(\d+)([KMGT]?)$/);
  return match && Number(match[1]) > 0 ? normalized : '64M';
}

function hasSharedMemory() {
  return typeof SharedArrayBuffer !== 'undefined'
    && typeof crossOriginIsolated !== 'undefined'
    && crossOriginIsolated;
}

function mergeSequenceMetadata(records, order, lengths) {
  for (const record of records || []) {
    if (lengths.has(record.name)) continue;
    lengths.set(record.name, record.length);
    order.push(record.name);
  }
}

function offsetsFor(names, lenByName) {
  const offsets = new Map();
  let offset = 0;
  for (const name of names) {
    offsets.set(name, offset);
    offset += lenByName.get(name) || 0;
  }
  return offsets;
}

function totalLength(names, lenByName) {
  return Math.max(1, names.reduce((sum, name) => sum + (lenByName.get(name) || 0), 0));
}

configureDefaults();
downloadTabBtn?.addEventListener('click', downloadLatestTab);
runBtn?.addEventListener('click', () => {
  run();
});
