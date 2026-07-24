import { JSLAST_VERSION } from './version.js?v=20260724-plot-profile';

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
const lastdbArgsInput = $('#lastdbArgs');
const lastalArgsInput = $('#lastalArgs');
const safeMode = $('#safeMode');
const TIMING_MODEL_KEY = `lastjsTimingModel:${JSLAST_VERSION}`;
const THREAD_PROGRESS_COMPLETED = 0;
const THREAD_PROGRESS_STARTED_BATCHES = 1;
const THREAD_PROGRESS_THREADS = 2;
const THREAD_PROGRESS_TOTAL_BATCHES = 3;
const THREAD_PROGRESS_DONE = 4;
const MAX_PLOT_WORKERS = 8;
const MAX_PLOT_ALIGNMENTS = 100000;
const MAX_UNSHARED_PLOT_BYTES = 32 * 1024 * 1024;

let runClock = null;
let estimateState = null;
let threadProgressView = null;
let threadProgressTimer = null;
let latestTabBuffer = null;

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

  if (lastalArgsInput && !lastalArgsInput.value) {
    lastalArgsInput.value = '-m100';
  }

  const updateSafeMode = () => {
    if (!lastdbArgsInput || !safeMode) return;
    if (safeMode.checked) {
      lastdbArgsInput.value = '';
      lastdbArgsInput.placeholder = 'Uses LAST-compatible database defaults';
      lastdbArgsInput.setAttribute('disabled', 'disabled');
    } else {
      lastdbArgsInput.removeAttribute('disabled');
      lastdbArgsInput.placeholder = '--bits=4 -R00 -uNEAR -w1 -W1 -S1 -C1 -v';
    }
  };
  safeMode?.addEventListener('change', updateSafeMode);
  updateSafeMode();
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
    appendLog(`$ lastdb + lastal ${lastalArgsInput?.value || ''}`);
    startEstimateClock({
      refText,
      qryText,
      requestedThreads,
    });

    const result = await runJsLastInWorker({
      refText,
      qryText,
      lastdbArgs: lastdbArgsInput?.value || '',
      lastalArgs: lastalArgsInput?.value || '',
      requestedThreads,
      useDefaultDbArgs: Boolean(safeMode?.checked),
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
    const plotStats = await renderPlotlyFromTab(result.tabBuffer, requestedThreads);
    for (const line of formatPlotProfile(plotStats)) appendLog(line);
    if (plotStats.workers > 0) {
      appendLog(`dotplot-parser workers=${plotStats.workers} alignments=${plotStats.alignments} plotted=${plotStats.plotted}`);
      if (plotStats.plotted < plotStats.alignments) {
        appendLog(`[WARN] Dot plot sampled ${plotStats.plotted} of ${plotStats.alignments} alignments to limit memory use.`);
      }
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
    const worker = new Worker(new URL(`./jslast-runner-worker.js?v=${JSLAST_VERSION}`, import.meta.url), {
      type: 'module',
      name: 'jslast-runner',
    });

    worker.onmessage = event => {
      const data = event.data || {};
      if (data.type === 'progress') {
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
  if (window.Plotly) window.Plotly.purge(plotDiv);
  plotDiv.textContent = '';
}

async function renderPlotlyFromTab(tabBuffer, requestedWorkers) {
  const totalStartedAt = performance.now();
  const memoryStartMiB = usedHeapMiB();
  const plotDiv = $('#dotplot');
  if (!plotDiv || !tabBuffer) return { workers: 0, alignments: 0, plotted: 0 };
  if (!window.Plotly) {
    plotDiv.textContent = 'Plotly is unavailable; TAB output is shown above.';
    return { workers: 0, alignments: 0, plotted: 0 };
  }
  const canShare = hasSharedMemory();
  if (!canShare && tabBuffer.byteLength > MAX_UNSHARED_PLOT_BYTES) {
    plotDiv.textContent = 'Dot plot skipped because shared memory is unavailable for this large result.';
    return { workers: 0, alignments: 0, plotted: 0 };
  }

  const plotData = await parsePlotDataInParallel(tabBuffer, requestedWorkers, canShare);
  if (plotData.alignmentCount === 0) {
    plotDiv.textContent = 'No plottable alignments.';
    return { workers: plotData.workerCount, alignments: 0, plotted: 0 };
  }

  const hovertemplate = [
    'Score: %{customdata.score}',
    'Ref: %{customdata.refName}:%{customdata.refStart}-%{customdata.refEnd}',
    'Qry: %{customdata.qryName}:%{customdata.qryStart}-%{customdata.qryEnd}',
    'Length: %{customdata.len}',
    '<extra></extra>',
  ].join('<br>');
  const tracesStartedAt = performance.now();
  const traces = plotTraces(plotData.chunks, hovertemplate);
  const tracesMs = performance.now() - tracesStartedAt;

  const layoutStartedAt = performance.now();
  const layout = {
    xaxis: { title: 'Reference', range: [0, plotData.maxX], zeroline: false, showgrid: true },
    yaxis: { title: 'Query', range: [0, plotData.maxY], scaleanchor: 'x', scaleratio: 1, zeroline: false, showgrid: true },
    shapes: boundaryShapes(
      plotData.refOrders,
      plotData.refOffsets,
      plotData.qryOrders,
      plotData.qryOffsets,
      plotData.maxX,
      plotData.maxY,
    ),
    annotations: sequenceAnnotations(
      plotData.refOrders,
      plotData.refOffsets,
      plotData.refLenByName,
      plotData.qryOrders,
      plotData.qryOffsets,
      plotData.qryLenByName,
    ),
    showlegend: true,
    dragmode: 'pan',
    hovermode: 'closest',
    margin: { l: 70, r: 20, t: 20, b: 60 },
  };
  const config = { responsive: true, scrollZoom: true, displaylogo: false, modeBarButtonsToRemove: ['select2d', 'lasso2d'] };
  const layoutMs = performance.now() - layoutStartedAt;
  const memoryBeforePlotlyMiB = usedHeapMiB();
  const longTaskProfile = startLongTaskProfile();
  const plotlyStartedAt = performance.now();
  const plotPromise = window.Plotly.newPlot(plotDiv, traces, layout, config);
  const plotlyCallMs = performance.now() - plotlyStartedAt;
  const plotlyPromiseStartedAt = performance.now();
  await plotPromise;
  const plotlyPromiseWaitMs = performance.now() - plotlyPromiseStartedAt;
  const plotlyNewPlotMs = performance.now() - plotlyStartedAt;
  const paintTiming = await measureFinalPlotPaint();
  const longTasks = longTaskProfile.stop();
  return {
    workers: plotData.workerCount,
    alignments: plotData.alignmentCount,
    plotted: plotData.plottedCount,
    timing: {
      ...plotData.timing,
      tracesMs,
      layoutMs,
      plotlyCallMs,
      plotlyPromiseWaitMs,
      plotlyNewPlotMs,
      finalPaintMs: paintTiming.totalMs,
      eventLoopYieldMs: paintTiming.eventLoopYieldMs,
      firstFrameMs: paintTiming.firstFrameMs,
      secondFrameMs: paintTiming.secondFrameMs,
      longTasks,
      totalMs: performance.now() - totalStartedAt,
    },
    memory: {
      startMiB: memoryStartMiB,
      ...plotData.memory,
      beforePlotlyMiB: memoryBeforePlotlyMiB,
      afterPlotlyMiB: usedHeapMiB(),
    },
  };
}

async function parsePlotDataInParallel(tabBuffer, requestedWorkers, canShare) {
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
      },
      canShare || index > 0 ? [] : [buffer],
    )));
    const parsePassMs = performance.now() - parseStartedAt;
    const parseWorkers = summarizeWorkerTimings(metadata);
    const afterParseMiB = usedHeapMiB();
    const mergeStartedAt = performance.now();
    const merged = mergePlotMetadata(metadata);
    const sampleEvery = Math.max(1, Math.ceil(merged.alignmentCount / MAX_PLOT_ALIGNMENTS));
    const metadataMergeMs = performance.now() - mergeStartedAt;
    const buildStartedAt = performance.now();
    const chunks = await Promise.all(sessions.map(session => requestPlotWorker(session.worker, {
      type: 'build',
      index: session.index,
      refOffsets: Array.from(merged.refOffsets),
      qryOffsets: Array.from(merged.qryOffsets),
      alignmentStart: merged.alignmentStartByIndex.get(session.index) || 0,
      sampleEvery,
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
      plottedCount: chunks.reduce((sum, chunk) => sum + (chunk.plottedCount || 0), 0),
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

async function measureFinalPlotPaint() {
  const startedAt = performance.now();
  const eventLoopStartedAt = performance.now();
  await new Promise(resolve => setTimeout(resolve, 0));
  const eventLoopYieldMs = performance.now() - eventLoopStartedAt;
  const firstFrameStartedAt = performance.now();
  await nextAnimationFrame();
  const firstFrameMs = performance.now() - firstFrameStartedAt;
  const secondFrameStartedAt = performance.now();
  await nextAnimationFrame();
  const secondFrameMs = performance.now() - secondFrameStartedAt;
  return {
    totalMs: performance.now() - startedAt,
    eventLoopYieldMs,
    firstFrameMs,
    secondFrameMs,
  };
}

function nextAnimationFrame() {
  return new Promise(resolve => {
    let settled = false;
    const finish = () => {
      if (settled) return;
      settled = true;
      resolve();
    };
    setTimeout(finish, 250);
    requestAnimationFrame(finish);
  });
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
    `dotplot-timing total=${formatMs(timing.totalMs)} parser-total=${formatMs(timing.parserTotalMs)} plotly-newPlot=${formatMs(timing.plotlyNewPlotMs)} final-paint=${formatMs(timing.finalPaintMs)}`,
    `dotplot-timing setup buffer-copy=${formatMs(timing.bufferCopyMs)} range-split=${formatMs(timing.rangeSplitMs)} worker-create=${formatMs(timing.workerCreateMs)} metadata-merge=${formatMs(timing.metadataMergeMs)} finalize=${formatMs(timing.finalizeMs)}`,
    `dotplot-timing parse-pass=${formatMs(timing.parsePassMs)} worker-max=${formatMs(parse.maxTotalMs)} decode-max=${formatMs(parse.maxDecodeMs)} scan-max=${formatMs(parse.maxScanMs)} startup-message=${formatMs(timing.parseOverheadMs)}`,
    `dotplot-timing build-pass=${formatMs(timing.buildPassMs)} worker-max=${formatMs(build.maxTotalMs)} decode-max=${formatMs(build.maxDecodeMs)} scan-max=${formatMs(build.maxScanMs)} arrays-max=${formatMs(build.maxArraysMs)} message-clone=${formatMs(timing.buildOverheadMs)}`,
    `dotplot-timing main traces=${formatMs(timing.tracesMs)} layout=${formatMs(timing.layoutMs)}`,
    `dotplot-timing plotly call=${formatMs(timing.plotlyCallMs)} promise-wait=${formatMs(timing.plotlyPromiseWaitMs)} event-loop=${formatMs(timing.eventLoopYieldMs)} frame-1=${formatMs(timing.firstFrameMs)} frame-2=${formatMs(timing.secondFrameMs)}`,
    `dotplot-longtasks count=${timing.longTasks?.count || 0} total=${formatMs(timing.longTasks?.totalMs)} max=${formatMs(timing.longTasks?.maxMs)}`,
  ];
  const memory = plotStats.memory || {};
  if (Number.isFinite(memory.startMiB)) {
    lines.push(`dotplot-memory heap start=${formatMiB(memory.startMiB)} after-buffer=${formatMiB(memory.afterBufferMiB)} after-parse=${formatMiB(memory.afterParseMiB)} after-build=${formatMiB(memory.afterBuildMiB)} before-plotly=${formatMiB(memory.beforePlotlyMiB)} after-plotly=${formatMiB(memory.afterPlotlyMiB)}`);
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

function mergePlotMetadata(metadata) {
  const refOrders = [];
  const qryOrders = [];
  const refLenByName = new Map();
  const qryLenByName = new Map();
  const alignmentStartByIndex = new Map();
  let alignmentCount = 0;
  metadata.sort((left, right) => left.index - right.index);
  for (const chunk of metadata) {
    alignmentStartByIndex.set(chunk.index, alignmentCount);
    alignmentCount += chunk.alignmentCount || 0;
    mergeSequenceMetadata(chunk.refs, refOrders, refLenByName);
    mergeSequenceMetadata(chunk.queries, qryOrders, qryLenByName);
  }
  const refOffsets = offsetsFor(refOrders, refLenByName);
  const qryOffsets = offsetsFor(qryOrders, qryLenByName);
  return {
    alignmentCount,
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

function plotTraces(chunks, hovertemplate) {
  const traces = [];
  let hasForwardLegend = false;
  let hasReverseLegend = false;
  for (const chunk of chunks) {
    if (chunk.forward.count > 0) {
      traces.push(plotTrace('Forward', '#1664d9', hovertemplate, chunk.forward, !hasForwardLegend));
      hasForwardLegend = true;
    }
    if (chunk.reverse.count > 0) {
      traces.push(plotTrace('Reverse', '#c92a2a', hovertemplate, chunk.reverse, !hasReverseLegend));
      hasReverseLegend = true;
    }
  }
  return traces;
}

function plotTrace(name, color, hovertemplate, chunk, showlegend) {
  return {
    x: new Float64Array(chunk.x),
    y: new Float64Array(chunk.y),
    customdata: chunk.customdata,
    mode: 'lines+markers',
    type: 'scattergl',
    name,
    legendgroup: name,
    showlegend,
    line: { color, width: 1.8 },
    marker: { size: 5, opacity: 0 },
    hovertemplate,
  };
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

function boundaryShapes(refOrders, refOffsets, qryOrders, qryOffsets, maxX, maxY) {
  const shapes = [];
  for (let i = 1; i < refOrders.length; i++) {
    const off = refOffsets.get(refOrders[i]) || 0;
    shapes.push({ type: 'line', x0: off, x1: off, y0: 0, y1: maxY, line: { color: 'rgba(0,0,0,0.2)', width: 1, dash: 'dot' } });
  }
  for (let i = 1; i < qryOrders.length; i++) {
    const off = qryOffsets.get(qryOrders[i]) || 0;
    shapes.push({ type: 'line', x0: 0, x1: maxX, y0: off, y1: off, line: { color: 'rgba(0,0,0,0.2)', width: 1, dash: 'dot' } });
  }
  return shapes;
}

function sequenceAnnotations(refOrders, refOffsets, refLenByName, qryOrders, qryOffsets, qryLenByName) {
  const annotations = [];
  for (const name of refOrders) {
    const off = refOffsets.get(name) || 0;
    const len = refLenByName.get(name) || 0;
    annotations.push({ x: off + len / 2, y: 1, xref: 'x', yref: 'paper', yanchor: 'bottom', yshift: 8, text: name, showarrow: false, font: { size: 10 }, textangle: -35 });
  }
  for (const name of qryOrders) {
    const off = qryOffsets.get(name) || 0;
    const len = qryLenByName.get(name) || 0;
    annotations.push({ x: 0, y: off + len / 2, xref: 'paper', yref: 'y', xanchor: 'right', xshift: -8, text: name, showarrow: false, font: { size: 10 } });
  }
  return annotations;
}

configureDefaults();
downloadTabBtn?.addEventListener('click', downloadLatestTab);
runBtn?.addEventListener('click', () => {
  run();
});
