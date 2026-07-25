import { runJsLast } from './jslast.js?v=20260725-axis-labels';

const TAB_PREVIEW_CHARS = 20 * 1024;
const PREVIEW_PAINT_TIMEOUT_MS = 5000;

let previewPaintResolver = null;

self.onmessage = async event => {
  const data = event.data || {};
  if (data.type === 'tab-preview-painted') {
    previewPaintResolver?.();
    return;
  }
  if (data.type !== 'run') return;

  try {
    let previewSent = false;
    const sendTabPreview = async ({ text = '', totalChars = text.length } = {}) => {
      if (previewSent) return;
      previewSent = true;
      const previewEnd = previewBoundary(text, TAB_PREVIEW_CHARS, totalChars);
      const paintPromise = waitForPreviewPaint();
      self.postMessage({
        type: 'tab-preview',
        text: text.slice(0, previewEnd),
        truncated: totalChars > previewEnd,
        totalChars,
      });
      await paintPromise;
    };

    const result = await runJsLast({
      ...data.payload,
      tabPreviewChars: TAB_PREVIEW_CHARS,
      onTabPreview: sendTabPreview,
      onProgress: progressEvent => {
        self.postMessage({
          type: 'progress',
          event: sanitizeProgress(progressEvent),
        });
      },
    });
    let tabText = result.tabText || '';
    await sendTabPreview({ text: tabText, totalChars: tabText.length });

    self.postMessage({
      type: 'progress',
      event: { stage: 'tab-encode', status: 'start', message: 'Encoding complete TAB for transfer...' },
    });
    const encodeStartedAt = performance.now();
    const tabBytes = new TextEncoder().encode(tabText);
    result.timing.tabEncodeMs = performance.now() - encodeStartedAt;
    self.postMessage({
      type: 'progress',
      event: {
        stage: 'tab-encode',
        status: 'complete',
        message: 'TAB encoded; transferring result to the page...',
        elapsedMs: result.timing.tabEncodeMs,
        tabBytes: tabBytes.byteLength,
      },
    });
    delete result.tabText;
    tabText = '';
    self.postMessage({
      type: 'done',
      result: {
        ...sanitizeResult(result),
        tabByteLength: tabBytes.byteLength,
      },
      tabBuffer: tabBytes.buffer,
    }, [tabBytes.buffer]);
  } catch (error) {
    self.postMessage({
      type: 'error',
      error: error?.stack || error?.message || String(error),
    });
  }
};

function waitForPreviewPaint() {
  return new Promise(resolve => {
    let timeoutId = null;
    const finish = () => {
      if (previewPaintResolver !== finish) return;
      previewPaintResolver = null;
      if (timeoutId !== null) clearTimeout(timeoutId);
      resolve();
    };
    previewPaintResolver = finish;
    timeoutId = setTimeout(finish, PREVIEW_PAINT_TIMEOUT_MS);
  });
}

function previewBoundary(text, limit, totalChars = text.length) {
  if (totalChars <= limit) return text.length;
  const newline = text.lastIndexOf('\n', limit);
  return newline > 0 ? newline + 1 : limit;
}

function sanitizeProgress(event) {
  const clean = { ...event };
  if (event.refs) clean.refs = event.refs.map(sequenceSummary);
  if (event.queries) clean.queries = event.queries.map(sequenceSummary);
  return clean;
}

function sanitizeResult(result) {
  return {
    ...result,
    refs: (result.refs || []).map(sequenceSummary),
    queries: (result.queries || []).map(sequenceSummary),
  };
}

function sequenceSummary(record) {
  return {
    name: record.name,
    length: record.seq?.length || record.length || 0,
  };
}
