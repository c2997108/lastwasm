let chunkBuffer = null;
let chunkStart = 0;
let chunkEnd = 0;

self.onmessage = event => {
  const data = event.data || {};
  try {
    if (data.type === 'parse') {
      parseChunk(data);
    } else if (data.type === 'build') {
      buildPlotChunks(data);
    }
  } catch (error) {
    self.postMessage({
      type: 'error',
      error: error?.stack || error?.message || String(error),
    });
  }
};

function parseChunk({ index, buffer, start, end }) {
  chunkBuffer = buffer;
  chunkStart = start;
  chunkEnd = end;
  const refLengths = new Map();
  const qryLengths = new Map();
  let alignmentCount = 0;

  forEachAlignment(segment => {
    if (!refLengths.has(segment.refName)) refLengths.set(segment.refName, segment.refTotal);
    if (!qryLengths.has(segment.qryName)) qryLengths.set(segment.qryName, segment.qryTotal);
    alignmentCount += 1;
  });

  self.postMessage({
    type: 'parsed',
    index,
    alignmentCount,
    refs: Array.from(refLengths, ([name, length]) => ({ name, length })),
    queries: Array.from(qryLengths, ([name, length]) => ({ name, length })),
  });
}

function buildPlotChunks({ index, refOffsets, qryOffsets, alignmentStart = 0, sampleEvery = 1 }) {
  const refOffsetByName = new Map(refOffsets);
  const qryOffsetByName = new Map(qryOffsets);
  const selected = [];
  let localIndex = 0;
  forEachAlignment(segment => {
    if ((alignmentStart + localIndex) % sampleEvery === 0) selected.push(segment);
    localIndex += 1;
  });

  const forwardCount = selected.reduce((count, segment) => count + (segment.qryStrand === '-' ? 0 : 1), 0);
  const reverseCount = selected.length - forwardCount;
  const forward = makePlotChunk(selected, false, forwardCount, refOffsetByName, qryOffsetByName);
  const reverse = makePlotChunk(selected, true, reverseCount, refOffsetByName, qryOffsetByName);
  chunkBuffer = null;

  self.postMessage({
    type: 'built',
    index,
    plottedCount: selected.length,
    forward,
    reverse,
  }, [forward.x, forward.y, reverse.x, reverse.y]);
}

function forEachAlignment(callback) {
  const text = decodeChunk();
  for (let start = 0; start < text.length;) {
    const newline = text.indexOf('\n', start);
    const end = newline >= 0 ? newline : text.length;
    const segment = parseAlignment(text.slice(start, end).trim());
    if (segment) callback(segment);
    start = newline >= 0 ? newline + 1 : text.length;
  }
}

function decodeChunk() {
  const source = new Uint8Array(chunkBuffer, chunkStart, chunkEnd - chunkStart);
  const bytes = typeof SharedArrayBuffer !== 'undefined' && chunkBuffer instanceof SharedArrayBuffer
    ? source.slice()
    : source;
  return new TextDecoder().decode(bytes);
}

function parseAlignment(line) {
  if (!/^\d+\t/.test(line)) return null;
  const cols = line.split(/\s+/);
  if (cols.length < 12) return null;
  const segment = {
    score: Number(cols[0]),
    refName: cols[1],
    refStart: Number(cols[2]),
    refLen: Number(cols[3]),
    refTotal: Number(cols[5]),
    qryName: cols[6],
    qryStart: Number(cols[7]),
    qryLen: Number(cols[8]),
    qryStrand: cols[9],
    qryTotal: Number(cols[10]),
  };
  return [
    segment.score,
    segment.refStart,
    segment.refLen,
    segment.refTotal,
    segment.qryStart,
    segment.qryLen,
    segment.qryTotal,
  ].some(value => !Number.isFinite(value)) ? null : segment;
}

function makePlotChunk(chunkSegments, reverseTarget, segmentCount, refOffsets, qryOffsets) {
  const pointCount = segmentCount * 3;
  const x = new Float64Array(pointCount);
  const y = new Float64Array(pointCount);
  const customdata = new Array(pointCount);
  let pointIndex = 0;

  for (const segment of chunkSegments) {
    if ((segment.qryStrand === '-') !== reverseTarget) continue;
    const refOffset = refOffsets.get(segment.refName) || 0;
    const qryOffset = qryOffsets.get(segment.qryName) || 0;
    const reverseQuery = segment.qryStrand === '-';
    const x1 = refOffset + segment.refStart + 1;
    const x2 = refOffset + segment.refStart + segment.refLen;
    const y1 = qryOffset + (reverseQuery
      ? segment.qryTotal - segment.qryStart
      : segment.qryStart + 1);
    const y2 = qryOffset + (reverseQuery
      ? segment.qryTotal - segment.qryStart - segment.qryLen + 1
      : segment.qryStart + segment.qryLen);
    const hoverData = {
      score: segment.score,
      refName: segment.refName,
      refStart: segment.refStart + 1,
      refEnd: segment.refStart + segment.refLen,
      qryName: segment.qryName,
      qryStart: reverseQuery ? segment.qryTotal - segment.qryStart : segment.qryStart + 1,
      qryEnd: reverseQuery
        ? segment.qryTotal - segment.qryStart - segment.qryLen + 1
        : segment.qryStart + segment.qryLen,
      len: Math.min(segment.refLen, segment.qryLen),
    };

    x[pointIndex] = x1;
    x[pointIndex + 1] = x2;
    x[pointIndex + 2] = NaN;
    y[pointIndex] = y1;
    y[pointIndex + 1] = y2;
    y[pointIndex + 2] = NaN;
    customdata[pointIndex] = hoverData;
    customdata[pointIndex + 1] = hoverData;
    customdata[pointIndex + 2] = null;
    pointIndex += 3;
  }

  return {
    count: segmentCount,
    x: x.buffer,
    y: y.buffer,
    customdata,
  };
}
