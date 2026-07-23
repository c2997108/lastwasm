import { runJsLast } from './jslast.js?v=20260724-large-results';

const TAB_PREVIEW_CHARS = 2 * 1024 * 1024;

self.onmessage = async event => {
  const data = event.data || {};
  if (data.type !== 'run') return;

  try {
    const result = await runJsLast({
      ...data.payload,
      onProgress: progressEvent => {
        self.postMessage({
          type: 'progress',
          event: sanitizeProgress(progressEvent),
        });
      },
    });
    let tabText = result.tabText || '';
    const previewEnd = previewBoundary(tabText, TAB_PREVIEW_CHARS);
    self.postMessage({
      type: 'tab-preview',
      text: tabText.slice(0, previewEnd),
      truncated: previewEnd < tabText.length,
      totalChars: tabText.length,
    });

    const tabBytes = new TextEncoder().encode(tabText);
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

function previewBoundary(text, limit) {
  if (text.length <= limit) return text.length;
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
