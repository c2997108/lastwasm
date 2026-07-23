import { JSLAST_VERSION } from './version.js?v=20260724-wasm-pthreads';

window.__JSLAST_MODULE_READY = true;

const $ = (sel) => document.querySelector(sel);
const statusEl = $('#status');
const logEl = $('#log');
const tabEl = $('#tab');
const timeEl = $('#timeEstimate');
const runBtn = $('#runBtn');
const workerInput = $('#workerCount');
const lastdbArgsInput = $('#lastdbArgs');
const lastalArgsInput = $('#lastalArgs');
const safeMode = $('#safeMode');
const TIMING_MODEL_KEY = `lastjsTimingModel:${JSLAST_VERSION}`;

let runClock = null;
let estimateState = null;

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
    lastalArgsInput.value = '-m1000';
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
    });

    if (tabEl) tabEl.textContent = result.tabText;
    updateTimingModel(result.timing);
    if (result.runtime.endsWith('worker-pool')) {
      appendLog(`lastal-workers started=${result.searchWorkerStats.started} completed=${result.searchWorkerStats.completed}`);
    } else if (result.runtime === 'wasm-pthreads') {
      appendLog(`lastal-pthreads spawned=${result.alThreadStats.spawned} maxRunning=${result.alThreadStats.maxRunning}`);
    }
    appendLog(`threads=${result.threads}, alignments=${result.alignments.length}`);
    renderPlotlyFromTab(result.tabText);
    setStatus('Done');
    finishEstimateClock(result.timing);
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
      worker.terminate();
      if (data.type === 'done') {
        resolve(data.result);
      } else {
        reject(new Error(data.error || 'LAST worker failed'));
      }
    };

    worker.onerror = error => {
      worker.terminate();
      reject(new Error(error.message || 'LAST worker failed'));
    };

    const { onProgress, ...runPayload } = payload;
    worker.postMessage({ type: 'run', payload: runPayload });
  });
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

function clearPlot() {
  const plotDiv = $('#dotplot');
  if (!plotDiv) return;
  if (window.Plotly) window.Plotly.purge(plotDiv);
  plotDiv.textContent = '';
}

function renderPlotlyFromTab(tabText) {
  const plotDiv = $('#dotplot');
  if (!plotDiv || !tabText) return;
  if (!window.Plotly) {
    plotDiv.textContent = 'Plotly is unavailable; TAB output is shown above.';
    return;
  }

  const segments = [];
  const refOrders = [];
  const qryOrders = [];
  const refLenByName = new Map();
  const qryLenByName = new Map();

  for (const rawLine of tabText.split(/\n/)) {
    const line = rawLine.trim();
    if (!line || line.startsWith('#') || line.startsWith('$') || line.startsWith('[')) continue;
    const cols = line.split(/\s+/);
    if (cols.length < 12) continue;

    const score = Number(cols[0]);
    const refName = cols[1];
    const refStart = Number(cols[2]);
    const refLen = Number(cols[3]);
    const refStrand = cols[4];
    const refTotal = Number(cols[5]);
    const qryName = cols[6];
    const qryStart = Number(cols[7]);
    const qryLen = Number(cols[8]);
    const qryStrand = cols[9];
    const qryTotal = Number(cols[10]);
    if ([score, refStart, refLen, refTotal, qryStart, qryLen, qryTotal].some(n => !Number.isFinite(n))) continue;

    if (!refLenByName.has(refName)) {
      refLenByName.set(refName, refTotal);
      refOrders.push(refName);
    }
    if (!qryLenByName.has(qryName)) {
      qryLenByName.set(qryName, qryTotal);
      qryOrders.push(qryName);
    }
    segments.push({ score, refName, qryName, refStart, refLen, refTotal, refStrand, qryStart, qryLen, qryTotal, qryStrand });
  }

  if (segments.length === 0) {
    plotDiv.textContent = 'No plottable alignments.';
    return;
  }

  const refOffsets = offsetsFor(refOrders, refLenByName);
  const qryOffsets = offsetsFor(qryOrders, qryLenByName);
  const maxX = totalLength(refOrders, refLenByName);
  const maxY = totalLength(qryOrders, qryLenByName);

  const hovertemplate = [
    'Score: %{customdata.score}',
    'Ref: %{customdata.refName}:%{customdata.refStart}-%{customdata.refEnd}',
    'Qry: %{customdata.qryName}:%{customdata.qryStart}-%{customdata.qryEnd}',
    'Length: %{customdata.len}',
    '<extra></extra>',
  ].join('<br>');
  const forward = trace('Forward', '#1664d9', hovertemplate);
  const reverse = trace('Reverse', '#c92a2a', hovertemplate);

  for (const s of segments) {
    const roff = refOffsets.get(s.refName) || 0;
    const qoff = qryOffsets.get(s.qryName) || 0;
    const reverseQuery = s.qryStrand === '-';
    const x1 = roff + s.refStart + 1;
    const x2 = roff + s.refStart + s.refLen;
    const y1 = qoff + (reverseQuery ? (s.qryTotal - s.qryStart) : (s.qryStart + 1));
    const y2 = qoff + (reverseQuery ? (s.qryTotal - s.qryStart - s.qryLen + 1) : (s.qryStart + s.qryLen));
    const target = reverseQuery ? reverse : forward;
    const customdata = {
      score: s.score,
      refName: s.refName,
      refStart: s.refStart + 1,
      refEnd: s.refStart + s.refLen,
      qryName: s.qryName,
      qryStart: reverseQuery ? (s.qryTotal - s.qryStart) : (s.qryStart + 1),
      qryEnd: reverseQuery ? (s.qryTotal - s.qryStart - s.qryLen + 1) : (s.qryStart + s.qryLen),
      len: Math.min(s.refLen, s.qryLen),
    };
    target.x.push(x1, x2, null);
    target.y.push(y1, y2, null);
    target.customdata.push(customdata, customdata, null);
  }

  const layout = {
    xaxis: { title: 'Reference', range: [0, maxX], zeroline: false, showgrid: true },
    yaxis: { title: 'Query', range: [0, maxY], scaleanchor: 'x', scaleratio: 1, zeroline: false, showgrid: true },
    shapes: boundaryShapes(refOrders, refOffsets, qryOrders, qryOffsets, maxX, maxY),
    annotations: sequenceAnnotations(refOrders, refOffsets, refLenByName, qryOrders, qryOffsets, qryLenByName),
    showlegend: true,
    dragmode: 'pan',
    hovermode: 'closest',
    margin: { l: 70, r: 20, t: 20, b: 60 },
  };
  const config = { responsive: true, scrollZoom: true, displaylogo: false, modeBarButtonsToRemove: ['select2d', 'lasso2d'] };
  window.Plotly.newPlot(plotDiv, [forward, reverse], layout, config);
}

function trace(name, color, hovertemplate) {
  return {
    x: [],
    y: [],
    customdata: [],
    mode: 'lines+markers',
    type: 'scattergl',
    name,
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
runBtn?.addEventListener('click', () => {
  run();
});
