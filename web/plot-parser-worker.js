let chunkBuffer = null;
let chunkStart = 0;
let chunkEnd = 0;
let chunkRenderedCount = 0;

self.onmessage = event => {
  const data = event.data || {};
  try {
    if (data.type === 'parse') {
      parseChunk(data);
    } else if (data.type === 'build') {
      buildPlotChunk(data);
    }
  } catch (error) {
    self.postMessage({
      type: 'error',
      error: error?.stack || error?.message || String(error),
    });
  }
};

function parseChunk({ index, buffer, start, end, lengthCutoff = 0 }) {
  const startedAt = performance.now();
  chunkBuffer = buffer;
  chunkStart = start;
  chunkEnd = end;
  const refLengths = new Map();
  const qryLengths = new Map();
  let alignmentCount = 0;
  let renderedCount = 0;

  const scanTiming = forEachAlignment(segment => {
    if (!refLengths.has(segment.refName)) refLengths.set(segment.refName, segment.refTotal);
    if (!qryLengths.has(segment.qryName)) qryLengths.set(segment.qryName, segment.qryTotal);
    alignmentCount += 1;
    if (alignmentLength(segment) > lengthCutoff) renderedCount += 1;
  });
  chunkRenderedCount = renderedCount;

  self.postMessage({
    type: 'parsed',
    index,
    alignmentCount,
    renderedCount,
    refs: Array.from(refLengths, ([name, length]) => ({ name, length })),
    queries: Array.from(qryLengths, ([name, length]) => ({ name, length })),
    timing: {
      totalMs: performance.now() - startedAt,
      decodeMs: scanTiming.decodeMs,
      scanMs: scanTiming.scanMs,
    },
  });
}

function buildPlotChunk({ index, refOffsets, qryOffsets, alignmentStart = 0, lengthCutoff = 0 }) {
  const startedAt = performance.now();
  const refOffsetByName = new Map(refOffsets);
  const qryOffsetByName = new Map(qryOffsets);
  const refIdByName = new Map(refOffsets.map(([name], id) => [name, id]));
  const qryIdByName = new Map(qryOffsets.map(([name], id) => [name, id]));
  const NameIdArray = Math.max(refOffsets.length, qryOffsets.length) <= 0xffff ? Uint16Array : Uint32Array;
  const positions = new Float32Array(chunkRenderedCount * 4);
  const alignmentIds = new Float32Array(chunkRenderedCount * 2);
  const scores = new Float64Array(chunkRenderedCount);
  const refNameIds = new NameIdArray(chunkRenderedCount);
  const qryNameIds = new NameIdArray(chunkRenderedCount);
  const refStarts = new Uint32Array(chunkRenderedCount);
  const refLengths = new Uint32Array(chunkRenderedCount);
  const qryStarts = new Uint32Array(chunkRenderedCount);
  const qryLengths = new Uint32Array(chunkRenderedCount);
  const strands = new Uint8Array(chunkRenderedCount * 2);
  let localIndex = 0;

  const scanTiming = forEachAlignment(segment => {
    if (alignmentLength(segment) <= lengthCutoff) return;
    const refOffset = refOffsetByName.get(segment.refName) || 0;
    const qryOffset = qryOffsetByName.get(segment.qryName) || 0;
    const reverseQuery = segment.qryStrand === '-';
    const positionIndex = localIndex * 4;
    positions[positionIndex] = refOffset + segment.refStart + 1;
    positions[positionIndex + 1] = qryOffset + (reverseQuery
      ? segment.qryTotal - segment.qryStart
      : segment.qryStart + 1);
    positions[positionIndex + 2] = refOffset + segment.refStart + segment.refLen;
    positions[positionIndex + 3] = qryOffset + (reverseQuery
      ? segment.qryTotal - segment.qryStart - segment.qryLen + 1
      : segment.qryStart + segment.qryLen);
    alignmentIds[localIndex * 2] = alignmentStart + localIndex;
    alignmentIds[localIndex * 2 + 1] = alignmentStart + localIndex;
    scores[localIndex] = segment.score;
    refNameIds[localIndex] = refIdByName.get(segment.refName) || 0;
    qryNameIds[localIndex] = qryIdByName.get(segment.qryName) || 0;
    refStarts[localIndex] = segment.refStart;
    refLengths[localIndex] = segment.refLen;
    qryStarts[localIndex] = segment.qryStart;
    qryLengths[localIndex] = segment.qryLen;
    strands[localIndex * 2] = reverseQuery ? 1 : 0;
    strands[localIndex * 2 + 1] = reverseQuery ? 1 : 0;
    localIndex += 1;
  });
  chunkBuffer = null;

  const arrays = {
    positions: positions.buffer,
    alignmentIds: alignmentIds.buffer,
    scores: scores.buffer,
    refNameIds: refNameIds.buffer,
    qryNameIds: qryNameIds.buffer,
    refStarts: refStarts.buffer,
    refLengths: refLengths.buffer,
    qryStarts: qryStarts.buffer,
    qryLengths: qryLengths.buffer,
    strands: strands.buffer,
  };
  self.postMessage({
    type: 'built',
    index,
    alignmentStart,
    alignmentCount: localIndex,
    nameIdBytes: NameIdArray.BYTES_PER_ELEMENT,
    arrays,
    timing: {
      totalMs: performance.now() - startedAt,
      decodeMs: scanTiming.decodeMs,
      scanMs: scanTiming.scanMs,
      arraysMs: Math.max(0, performance.now() - startedAt - scanTiming.decodeMs - scanTiming.scanMs),
    },
  }, Object.values(arrays));
}

function forEachAlignment(callback) {
  const decodeStartedAt = performance.now();
  const text = decodeChunk();
  const decodeMs = performance.now() - decodeStartedAt;
  const scanStartedAt = performance.now();
  for (let start = 0; start < text.length;) {
    const newline = text.indexOf('\n', start);
    const end = newline >= 0 ? newline : text.length;
    const segment = parseAlignment(text.slice(start, end).trim());
    if (segment) callback(segment);
    start = newline >= 0 ? newline + 1 : text.length;
  }
  return {
    decodeMs,
    scanMs: performance.now() - scanStartedAt,
  };
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

function alignmentLength(segment) {
  return Math.min(segment.refLen, segment.qryLen);
}
