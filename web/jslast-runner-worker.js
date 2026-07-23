import { runJsLast } from './jslast.js?v=20260724-wasm-pthreads';

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
    self.postMessage({
      type: 'done',
      result: sanitizeResult(result),
    });
  } catch (error) {
    self.postMessage({
      type: 'error',
      error: error?.stack || error?.message || String(error),
    });
  }
};

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
