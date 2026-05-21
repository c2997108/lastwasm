#!/usr/bin/env node

const fs = require('fs');

function usage() {
  console.error('usage: node util/compare-tab.js <baseline.tab> <candidate.tab> [--threshold N] [--reciprocal] [--offdiag]');
  process.exit(2);
}

const args = process.argv.slice(2);
if (args.length < 2) usage();

const baselineFile = args[0];
const candidateFile = args[1];
let threshold = 0.5;
let reciprocal = false;
let offdiag = false;

for (let i = 2; i < args.length; i++) {
  const arg = args[i];
  if (arg === '--threshold' && args[i + 1]) threshold = Number(args[++i]);
  else if (arg.startsWith('--threshold=')) threshold = Number(arg.split('=').slice(1).join('='));
  else if (arg === '--reciprocal') reciprocal = true;
  else if (arg === '--offdiag') offdiag = true;
  else usage();
}

if (!Number.isFinite(threshold) || threshold < 0 || threshold > 1) usage();

const baseline = maybeFilter(parseTab(baselineFile));
const candidate = maybeFilter(parseTab(candidateFile));

const baselineCovered = coverage(baseline, candidate, threshold, reciprocal);
const candidateCovered = coverage(candidate, baseline, threshold, reciprocal);

console.log(JSON.stringify({
  mode: reciprocal ? 'reciprocal-overlap' : 'smaller-interval-overlap',
  threshold,
  offdiag,
  baseline: summarize(baseline),
  candidate: summarize(candidate),
  baselineCoveredByCandidate: baselineCovered.summary,
  candidateCoveredByBaseline: candidateCovered.summary,
  baselineTopMissing: baselineCovered.missing.slice(0, 20).map(formatMissing),
  candidateTopMissing: candidateCovered.missing.slice(0, 20).map(formatMissing),
}, null, 2));

function maybeFilter(rows) {
  if (!offdiag) return rows;
  return rows.filter(row => !(row.qryStrand === '+' && Math.abs(row.refStart - row.qryStart) < 100));
}

function parseTab(file) {
  const text = readText(file);
  return text.split(/\r?\n/)
    .filter(line => /^\s*\d+\s+/.test(line))
    .map((line) => {
      const cols = line.trim().split(/\s+/);
      return {
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
        blocks: cols[11] || '',
        line: line.trim(),
      };
    })
    .filter(row => Number.isFinite(row.score)
      && Number.isFinite(row.refStart)
      && Number.isFinite(row.refLen)
      && Number.isFinite(row.qryStart)
      && Number.isFinite(row.qryLen));
}

function readText(file) {
  const bytes = fs.readFileSync(file);
  if (bytes[0] === 0xff && bytes[1] === 0xfe) return bytes.toString('utf16le');
  if (bytes[0] === 0xfe && bytes[1] === 0xff) return bytes.swap16().toString('utf16le');
  return bytes.toString('utf8');
}

function summarize(rows) {
  return {
    count: rows.length,
    scores: quantiles(rows.map(row => row.score)),
    lengths: quantiles(rows.map(row => Math.min(row.refLen, row.qryLen))),
    strands: rows.reduce((acc, row) => {
      acc[row.qryStrand] = (acc[row.qryStrand] || 0) + 1;
      return acc;
    }, {}),
  };
}

function quantiles(values) {
  if (values.length === 0) return { min: 0, p25: 0, median: 0, p75: 0, p95: 0, max: 0 };
  const sorted = [...values].sort((a, b) => a - b);
  const q = p => sorted[Math.min(sorted.length - 1, Math.floor((sorted.length - 1) * p))];
  return {
    min: q(0),
    p25: q(0.25),
    median: q(0.5),
    p75: q(0.75),
    p95: q(0.95),
    max: q(1),
  };
}

function coverage(base, compared, minOverlap, useReciprocal) {
  let covered = 0;
  const missing = [];

  for (const row of base) {
    let best = 0;
    for (const other of compared) {
      best = Math.max(best, overlap(row, other, useReciprocal));
    }
    if (best >= minOverlap) covered++;
    else missing.push({ row, best });
  }

  missing.sort((a, b) => b.row.score - a.row.score);
  return {
    summary: {
      covered,
      total: base.length,
      rate: base.length === 0 ? 0 : covered / base.length,
    },
    missing,
  };
}

function overlap(a, b, useReciprocal) {
  if (a.refName !== b.refName || a.qryName !== b.qryName || a.qryStrand !== b.qryStrand) return 0;
  const refOverlap = intervalOverlap(a.refStart, a.refStart + a.refLen, b.refStart, b.refStart + b.refLen);
  const qryOverlap = intervalOverlap(a.qryStart, a.qryStart + a.qryLen, b.qryStart, b.qryStart + b.qryLen);
  const overlapLen = Math.min(refOverlap, qryOverlap);
  if (overlapLen <= 0) return 0;

  if (useReciprocal) {
    const aLen = Math.max(1, Math.max(a.refLen, a.qryLen));
    const bLen = Math.max(1, Math.max(b.refLen, b.qryLen));
    return Math.min(overlapLen / aLen, overlapLen / bLen);
  }

  const smaller = Math.max(1, Math.min(a.refLen, a.qryLen, b.refLen, b.qryLen));
  return overlapLen / smaller;
}

function intervalOverlap(a0, a1, b0, b1) {
  return Math.max(0, Math.min(a1, b1) - Math.max(a0, b0));
}

function formatMissing(item) {
  return {
    bestOverlap: item.best,
    score: item.row.score,
    ref: `${item.row.refName}:${item.row.refStart}-${item.row.refStart + item.row.refLen}`,
    query: `${item.row.qryName}:${item.row.qryStart}-${item.row.qryStart + item.row.qryLen}${item.row.qryStrand}`,
    blocks: item.row.blocks,
  };
}
