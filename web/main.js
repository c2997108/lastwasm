// Minimal browser UI to run lastdb and lastal compiled to WASM.
// Requires: lastdb.js/.wasm and lastal.js/.wasm built via build_wasm.sh

const $ = (sel) => document.querySelector(sel);
const APP_VER = '20250823-3';
const statusEl = $('#status');
const logEl = $('#log');
const outEl = $('#out');
const tabEl = $('#tab');
const SAFE_DB_ARGS = ['--bits=4','-R00','-uNEAR','-w1','-W1','-S1','-C1','-v'];

// Plotly-based dotplot

function setStatus(msg) {
  statusEl.textContent = msg;
}
function append(el, line) {
  el.textContent += (line.endsWith('\n') ? line : line + '\n');
}
function splitArgs(str) {
  // naive whitespace split; good enough for simple flags
  if (!str) return [];
  return str.trim().split(/\s+/);
}

function tick() {
  return new Promise((resolve) => requestAnimationFrame(() => resolve()));
}

function parseTabAndPlotly(tabText) {
  const plotDiv = document.getElementById('dotplot');
  if (!plotDiv || !tabText) return;

  // Parse TAB lines and collect alignment segments and contig lengths
  const segments = [];
  const refOrders = [];
  const qryOrders = [];
  const refLenByName = new Map();
  const qryLenByName = new Map();
  const lines = tabText.split(/\n+/).map(s => s.trim()).filter(Boolean);
  for (const line of lines) {
    if (!line || line.startsWith('$') || line.startsWith('[')) continue;
    const cols = line.split(/\s+/);
    if (cols.length < 12) continue;
    // Expected columns (1-based LAST -f TAB):
    // 1:score 2:refName 3:refStart 4:refLen 5:refStrand 6:refTotalLen
    // 7:qryName 8:qryStart 9:qryLen 10:qryStrand 11:qryTotalLen 12:blocks
    const refName = cols[1];
    const refStart = parseInt(cols[2], 10);
    const refLen = parseInt(cols[3], 10);
    const refStrand = cols[4];
    const refTotal = parseInt(cols[5], 10);
    const qryName = cols[6];
    const qryStart = parseInt(cols[7], 10);
    const qryLen = parseInt(cols[8], 10);
    const qryStrand = cols[9];
    const qryTotal = parseInt(cols[10], 10);
    const blocksSpec = cols[11];
    if ([refStart, refLen, refTotal, qryStart, qryLen, qryTotal].some(n => !Number.isFinite(n))) continue;
    if (!refLenByName.has(refName)) {
      refLenByName.set(refName, refTotal);
      refOrders.push(refName);
    }
    if (!qryLenByName.has(qryName)) {
      qryLenByName.set(qryName, qryTotal);
      qryOrders.push(qryName);
    }
    // Store one segment per TAB line; plotting handles +/- orientation directly
    segments.push({ refName, qryName, refStart, refLen, qryStart, qryLen, refTotal, qryTotal, refStrand, qryStrand });
  }

  if (segments.length === 0) return;

  // Build offsets and extents
  const refOffsets = new Map();
  const qryOffsets = new Map();
  let acc = 0;
  for (const n of refOrders) { refOffsets.set(n, acc); acc += refLenByName.get(n) || 0; }
  const maxX = Math.max(acc, 1);
  acc = 0;
  for (const n of qryOrders) { qryOffsets.set(n, acc); acc += qryLenByName.get(n) || 0; }
  const maxY = Math.max(acc, 1);

  // Build Plotly traces (forward=blue, reverse=red)
  const hoverT = 'Ref: %{customdata.refName}:%{customdata.refStart}-%{customdata.refEnd}' +
                 '<br>Qry: %{customdata.qryName}:%{customdata.qryStart}-%{customdata.qryEnd}' +
                 '<br>Len: %{customdata.len} bp<extra></extra>';
  const fwd = { x: [], y: [], customdata: [], mode: 'lines+markers', type: 'scattergl', name: 'Forward', line: { color: 'blue', width: 1.5 }, marker: { size: 6, opacity: 0 }, hovertemplate: hoverT };
  const rev = { x: [], y: [], customdata: [], mode: 'lines+markers', type: 'scattergl', name: 'Reverse', line: { color: 'red', width: 1.5 }, marker: { size: 6, opacity: 0 }, hovertemplate: hoverT };

  for (const s of segments) {
    const roff = refOffsets.get(s.refName) || 0;
    const qoff = qryOffsets.get(s.qryName) || 0;
    const isRev = (s.qryStrand === '-');
    const x1 = roff + s.refStart;
    const x2 = roff + s.refStart + s.refLen;
    const y1 = qoff + (isRev ? (s.qryTotal - s.qryStart) : s.qryStart);
    const y2 = qoff + (isRev ? (s.qryTotal - s.qryStart - s.qryLen) : (s.qryStart + s.qryLen));
    const T = isRev ? rev : fwd;
    const cd = { refName: s.refName, refStart: s.refStart, refEnd: s.refStart + s.refLen, qryName: s.qryName, qryStart: isRev ? (s.qryTotal - s.qryStart) : s.qryStart, qryEnd: isRev ? (s.qryTotal - s.qryStart - s.qryLen) : (s.qryStart + s.qryLen), len: Math.min(s.refLen, s.qryLen) };
    T.x.push(x1, x2, null);
    T.y.push(y1, y2, null);
    T.customdata.push(cd, cd, null);
  }

  // Shapes for contig boundaries (dashed grid)
  const shapes = [];
  for (let i = 1; i < refOrders.length; i++) {
    const off = refOffsets.get(refOrders[i]) || 0;
    shapes.push({ type: 'line', x0: off, x1: off, y0: 0, y1: maxY, line: { color: 'rgba(0,0,0,0.2)', width: 1, dash: 'dot' } });
  }
  for (let i = 1; i < qryOrders.length; i++) {
    const off = qryOffsets.get(qryOrders[i]) || 0;
    shapes.push({ type: 'line', x0: 0, x1: maxX, y0: off, y1: off, line: { color: 'rgba(0,0,0,0.2)', width: 1, dash: 'dot' } });
  }

  // Annotations for contig names
  const annotations = [];
  for (const n of refOrders) {
    const off = refOffsets.get(n) || 0;
    const len = refLenByName.get(n) || 0;
    const mid = off + len / 2;
    annotations.push({ x: mid, y: 1, xref: 'x', yref: 'paper', yanchor: 'bottom', yshift: 5, text: n, showarrow: false, font: { size: 10 }, textangle: -45 });
  }
  for (const n of qryOrders) {
    const off = qryOffsets.get(n) || 0;
    const len = qryLenByName.get(n) || 0;
    const mid = off + len / 2;
    annotations.push({ x: 0, y: mid, xref: 'paper', yref: 'y', xanchor: 'right', xshift: -5, text: n, showarrow: false, font: { size: 10 } });
  }

  const layout = {
    xaxis: { title: 'Reference', range: [0, maxX], zeroline: false, showgrid: true },
    yaxis: { title: 'Query', range: [0, maxY], scaleanchor: 'x', scaleratio: 1, zeroline: false, showgrid: true },
    shapes,
    annotations,
    showlegend: true,
    dragmode: 'pan',
    hovermode: 'closest',
    margin: { l: 60, r: 10, t: 10, b: 40 }
  };

  const config = { responsive: true, scrollZoom: true, displaylogo: false, modeBarButtonsToRemove: ['select2d','lasso2d'] };
  // Render
  // eslint-disable-next-line no-undef
  Plotly.newPlot(plotDiv, [fwd, rev], layout, config);
}

function worldToScreenX(x) {
  const { l, r } = DOT.margin;
  const w = DOT.canvas.width - l - r;
  return l + (x - DOT.view.x0) / (DOT.view.x1 - DOT.view.x0) * w;
}
function worldToScreenY(y) {
  const { t, b } = DOT.margin;
  const h = DOT.canvas.height - t - b;
  return t + h - (y - DOT.view.y0) / (DOT.view.y1 - DOT.view.y0) * h;
}
function screenToWorldX(px) {
  const { l, r } = DOT.margin;
  const w = DOT.canvas.width - l - r;
  return DOT.view.x0 + (px - l) / w * (DOT.view.x1 - DOT.view.x0);
}
function screenToWorldY(py) {
  const { t, b } = DOT.margin;
  const h = DOT.canvas.height - t - b;
  return DOT.view.y0 + (t + h - py) / h * (DOT.view.y1 - DOT.view.y0);
}

function niceTickStep(span, pixels, targetPx = 80) {
  if (span <= 0) return 1;
  const approxTicks = Math.max(1, Math.floor(pixels / targetPx));
  const raw = span / approxTicks;
  const pow10 = Math.pow(10, Math.floor(Math.log10(raw)));
  const candidates = [1, 2, 5, 10].map(m => m * pow10);
  let best = candidates[0];
  for (const c of candidates) if (Math.abs(raw - c) < Math.abs(raw - best)) best = c;
  return best;
}

function drawDotplot() {
  if (!DOT.canvas) return;
  const ctx = DOT.ctx;
  const { l, r, t, b } = DOT.margin;
  const W = DOT.canvas.width, H = DOT.canvas.height;
  const plotW = W - l - r;
  const plotH = H - t - b;
  ctx.clearRect(0, 0, W, H);

  // Axes
  ctx.strokeStyle = '#444';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(l, t);
  ctx.lineTo(l, t + plotH);
  ctx.lineTo(l + plotW, t + plotH);
  ctx.stroke();

  // Ticks
  ctx.fillStyle = '#666';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  const xStep = niceTickStep(DOT.view.x1 - DOT.view.x0, plotW);
  const yStep = niceTickStep(DOT.view.y1 - DOT.view.y0, plotH);
  const xStart = Math.ceil(DOT.view.x0 / xStep) * xStep;
  const yStart = Math.ceil(DOT.view.y0 / yStep) * yStep;
  ctx.strokeStyle = '#ddd';
  ctx.lineWidth = 1;
  // vertical grid + x labels
  for (let xv = xStart; xv <= DOT.view.x1; xv += xStep) {
    const X = worldToScreenX(xv);
    ctx.beginPath();
    ctx.moveTo(X, t);
    ctx.lineTo(X, t + plotH);
    ctx.stroke();
    ctx.fillText(String(Math.round(xv)), X, t + plotH + 14);
  }
  // horizontal grid + y labels
  ctx.textAlign = 'right';
  for (let yv = yStart; yv <= DOT.view.y1; yv += yStep) {
    const Y = worldToScreenY(yv);
    ctx.beginPath();
    ctx.moveTo(l, Y);
    ctx.lineTo(l + plotW, Y);
    ctx.stroke();
    ctx.fillText(String(Math.round(yv)), l - 4, Y + 4);
  }

  // Sequence boundary grid lines (stronger, dashed)
  ctx.strokeStyle = '#eee';
  ctx.setLineDash([4, 3]);
  for (let i = 1; i < DOT.refOrders.length; i++) {
    const off = DOT.refOffsets.get(DOT.refOrders[i]) || 0;
    const X = worldToScreenX(off);
    ctx.beginPath(); ctx.moveTo(X, t); ctx.lineTo(X, t + plotH); ctx.stroke();
  }
  for (let i = 1; i < DOT.qryOrders.length; i++) {
    const off = DOT.qryOffsets.get(DOT.qryOrders[i]) || 0;
    const Y = worldToScreenY(off);
    ctx.beginPath(); ctx.moveTo(l, Y); ctx.lineTo(l + plotW, Y); ctx.stroke();
  }
  ctx.setLineDash([]);

  // Labels per sequence (centered)
  ctx.fillStyle = '#555';
  ctx.textAlign = 'center';
  for (const n of DOT.refOrders) {
    const off = DOT.refOffsets.get(n) || 0;
    const len = DOT.refLenByName.get(n) || 0;
    const cx = worldToScreenX(off + len / 2);
    if (cx >= l && cx <= l + plotW) ctx.fillText(n, cx, t + plotH + 28);
  }
  ctx.save();
  ctx.translate(12, 0);
  ctx.rotate(-Math.PI / 2);
  for (const n of DOT.qryOrders) {
    const off = DOT.qryOffsets.get(n) || 0;
    const len = DOT.qryLenByName.get(n) || 0;
    const cy = -(worldToScreenY(off + len / 2) - t);
    if (cy >= 0 && cy <= plotH) ctx.fillText(n, cy, 0);
  }
  ctx.restore();

  // Draw segments
  ctx.globalAlpha = 0.9;
  for (const s of DOT.segments) {
    const roff = DOT.refOffsets.get(s.refName) || 0;
    const qoff = DOT.qryOffsets.get(s.qryName) || 0;
    const isRev = (s.qryStrand === '-');
    const x1 = worldToScreenX(roff + s.refStart);
    const x2 = worldToScreenX(roff + s.refStart + s.refLen);
    const y1 = worldToScreenY(qoff + (isRev ? (s.qryTotal - s.qryStart) : s.qryStart));
    const y2 = worldToScreenY(qoff + (isRev ? (s.qryTotal - s.qryStart - s.qryLen) : (s.qryStart + s.qryLen)));
    // Skip if completely outside view
    if ((x1 < DOT.margin.l && x2 < DOT.margin.l) || (x1 > W - DOT.margin.r && x2 > W - DOT.margin.r)) continue;
    if ((y1 < DOT.margin.t && y2 < DOT.margin.t) || (y1 > H - DOT.margin.b && y2 > H - DOT.margin.b)) continue;
    // Color by orientation: forward=blue, reverse=red
    ctx.strokeStyle = isRev ? '#d12e2e' : '#2e7dd1';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();
  }
  ctx.globalAlpha = 1;
}

function attachDotplotEvents() {
  const c = DOT.canvas;
  // Wheel zoom
  c.addEventListener('wheel', (ev) => {
    if (!DOT.canvas) return;
    ev.preventDefault();
    const factor = ev.deltaY < 0 ? 0.9 : 1.1;
    const mx = ev.offsetX, my = ev.offsetY;
    const wx = screenToWorldX(mx);
    const wy = screenToWorldY(my);
    const { x0, x1, y0, y1 } = DOT.view;
    const nx0 = wx - (wx - x0) * factor;
    const nx1 = wx + (x1 - wx) * factor;
    const ny0 = wy - (wy - y0) * factor;
    const ny1 = wy + (y1 - wy) * factor;
    // Clamp
    const minSpan = 10; // bp
    DOT.view.x0 = Math.max(0, Math.min(nx0, nx1 - minSpan));
    DOT.view.x1 = Math.min(DOT.maxX, Math.max(nx1, DOT.view.x0 + minSpan));
    DOT.view.y0 = Math.max(0, Math.min(ny0, ny1 - minSpan));
    DOT.view.y1 = Math.min(DOT.maxY, Math.max(ny1, DOT.view.y0 + minSpan));
    drawDotplot();
  }, { passive: false });

  // Drag pan
  c.addEventListener('mousedown', (ev) => {
    DOT.dragging = true;
    DOT.lastMouse = { x: ev.offsetX, y: ev.offsetY };
  });
  window.addEventListener('mouseup', () => { DOT.dragging = false; });
  window.addEventListener('mousemove', (ev) => {
    if (!DOT.dragging) return;
    const dx = ev.offsetX - DOT.lastMouse.x;
    const dy = ev.offsetY - DOT.lastMouse.y;
    DOT.lastMouse = { x: ev.offsetX, y: ev.offsetY };
    const { l, r, t, b } = DOT.margin;
    const plotW = DOT.canvas.width - l - r;
    const plotH = DOT.canvas.height - t - b;
    const spanX = DOT.view.x1 - DOT.view.x0;
    const spanY = DOT.view.y1 - DOT.view.y0;
    const dWx = -dx / plotW * spanX;
    const dWy = dy / plotH * spanY; // invert y
    DOT.view.x0 = Math.max(0, Math.min(DOT.view.x0 + dWx, DOT.maxX - spanX));
    DOT.view.x1 = DOT.view.x0 + spanX;
    DOT.view.y0 = Math.max(0, Math.min(DOT.view.y0 + dWy, DOT.maxY - spanY));
    DOT.view.y1 = DOT.view.y0 + spanY;
    drawDotplot();
  });

  // Double-click to reset
  c.addEventListener('dblclick', () => {
    DOT.view = { x0: 0, x1: DOT.maxX, y0: 0, y1: DOT.maxY };
    drawDotplot();
  });

  // Tooltip showing sequence names at cursor position
  const tip = document.getElementById('dotTip');
  c.addEventListener('mousemove', (ev) => {
    if (!tip) return;
    const xw = screenToWorldX(ev.offsetX);
    const yw = screenToWorldY(ev.offsetY);
    const refName = findSeqAtCoord(DOT.refOrders, DOT.refOffsets, DOT.refLenByName, xw);
    const qryName = findSeqAtCoord(DOT.qryOrders, DOT.qryOffsets, DOT.qryLenByName, yw);
    // Only show inside plot area
    const { l, r, t, b } = DOT.margin;
    const inside = ev.offsetX >= l && ev.offsetX <= (c.width - r) && ev.offsetY >= t && ev.offsetY <= (c.height - b);
    if (inside && (refName || qryName)) {
      tip.style.display = 'block';
      tip.style.left = (ev.pageX + 12) + 'px';
      tip.style.top = (ev.pageY + 12) + 'px';
      tip.textContent = `${refName || '-'}  |  ${qryName || '-'}`;
    } else {
      tip.style.display = 'none';
    }
  });
  c.addEventListener('mouseleave', () => {
    const tip = document.getElementById('dotTip');
    if (tip) tip.style.display = 'none';
  });
}

function findSeqAtCoord(orders, offsets, lenByName, coord) {
  if (!orders || orders.length === 0) return null;
  for (const n of orders) {
    const off = offsets.get(n) || 0;
    const len = lenByName.get(n) || 0;
    if (coord >= off && coord < off + len) return n;
  }
  return null;
}

async function readFileInput(fileInput) {
  const f = fileInput.files && fileInput.files[0];
  if (!f) return null;
  return await f.text();
}

async function run() {
  logEl.textContent = '';
  outEl.textContent = '';
  if (tabEl) tabEl.textContent = '';
  setStatus('準備中...');
  const verEl = document.getElementById('appVer');
  if (verEl) verEl.textContent = `(v${APP_VER})`;

  const refText = await readFileInput($('#refFasta'));
  const qryText = await readFileInput($('#qryFasta'));
  const dbName = ($('#dbName').value || 'refdb').trim();

  if (!refText || !qryText) {
    setStatus('対象/クエリのFASTAを選択してください');
    return;
  }

  // Dynamically import the Emscripten ES-module factories
  let lastdbFactory, lastalFactory;
  try {
    [lastdbFactory, lastalFactory] = await Promise.all([
      import(`./lastdb.js?v=${APP_VER}`).then(m => m.default || m),
      import(`./lastal.js?v=${APP_VER}`).then(m => m.default || m),
    ]);
  } catch (e) {
    setStatus('エラー');
    append(logEl, 'モジュール読み込みエラー: ' + (e?.message || e));
    return;
  }

  // lastdb: build database
  setStatus('lastdb 実行中...');
  const lastdbModule = await lastdbFactory({
    noInitialRun: true, // prevent auto main() run on module load
    print: (line) => append(logEl, line),
    printErr: (line) => append(logEl, '[ERR] ' + line),
  });

  const dbFS = lastdbModule.FS;
  if (!dbFS.analyzePath('/work').exists) dbFS.mkdir('/work');
  dbFS.chdir('/work');
  dbFS.writeFile('ref.fa', refText);

  const safe = $('#safeMode')?.checked;
  const dbArgs = safe ? SAFE_DB_ARGS.slice() : splitArgs($('#lastdbArgs').value);
  // Enforce single-thread: use 1 thread in browser builds
  const lastdbArgv = ['-P', '1', ...dbArgs, dbName, 'ref.fa'];
  append(logEl, '$ lastdb ' + lastdbArgv.map(a => /\s/.test(a)? ('"'+a+'"'):a).join(' '));
  await tick();
  try {
    lastdbModule.callMain(lastdbArgv);
  } catch (e) {
    // Emscripten throws ExitStatus on exits; treat 0 as success
    const status = (typeof e?.status === 'number') ? e.status : undefined;
    if (!(e?.name === 'ExitStatus' && status === 0)) {
      append(logEl, `lastdb 実行エラー: ${e?.message || e}`);
      if (typeof status === 'number') append(logEl, `exit code: ${status}`);
      setStatus('エラー');
      return;
    }
  }

  // Collect generated DB files
  const entries = dbFS.readdir('.')
    .filter(n => n !== '.' && n !== '..' && (n === dbName || n.startsWith(dbName + '.')));
  append(logEl, '生成ファイル: ' + (entries.join(', ') || '(なし)'));

  // lastal: run alignment (TAB only)
  setStatus('lastal 実行中...');
  const alArgs = splitArgs($('#lastalArgs').value);
  if (tabEl) {
    const lastalModule2 = await lastalFactory({
      noInitialRun: true,
      print: (line) => append(tabEl, line),
      printErr: (line) => append(tabEl, '[ERR] ' + line),
    });
    const alFS2 = lastalModule2.FS;
    if (!alFS2.analyzePath('/work').exists) alFS2.mkdir('/work');
    alFS2.chdir('/work');
    // copy DB files into FS
    for (const name of entries) {
      const data = dbFS.readFile('/work/' + name);
      alFS2.writeFile('/work/' + name, data);
    }
    // write query
    alFS2.writeFile('query.fa', qryText);

    const lastalTabArgv = ['-P', '1', '-f', 'TAB', ...alArgs, dbName, 'query.fa'];
    append(tabEl, '$ lastal ' + lastalTabArgv.map(a => /\s/.test(a)? ('"'+a+'"'):a).join(' '));
    await tick();
    try {
      lastalModule2.callMain(lastalTabArgv);
    } catch (e) {
      const status = (typeof e?.status === 'number') ? e.status : undefined;
      if (!(e?.name === 'ExitStatus' && status === 0)) {
        append(tabEl, `lastal(TAB) 実行エラー: ${e?.message || e}`);
        if (typeof status === 'number') append(tabEl, `exit code: ${status}`);
        setStatus('エラー');
        return;
      }
    }
  }

  // Render dot plot from TAB text
  if (tabEl) parseTabAndPlotly(tabEl.textContent);

  setStatus('完了');
}

$('#runBtn').addEventListener('click', () => {
  run().catch(err => {
    setStatus('エラー');
    append(outEl, String(err?.stack || err));
  });
});

// lastdb only runner
async function runDbOnly() {
  logEl.textContent = '';
  outEl.textContent = '';
  setStatus('lastdb準備中...');

  const refText = await readFileInput($('#refFasta'));
  const dbName = ($('#dbName').value || 'refdb').trim();
  if (!refText) {
    setStatus('対象FASTAを選択してください');
    return;
  }

  let lastdbFactory;
  try {
    lastdbFactory = await import(`./lastdb.js?v=${APP_VER}`).then(m => m.default || m);
  } catch (e) {
    setStatus('エラー');
    append(logEl, 'モジュール読み込みエラー: ' + (e?.message || e));
    return;
  }

  const lastdbModule = await lastdbFactory({
    noInitialRun: true, // prevent auto main() run on module load
    print: (line) => append(logEl, line),
    printErr: (line) => append(logEl, '[ERR] ' + line),
  });
  const dbFS = lastdbModule.FS;
  if (!dbFS.analyzePath('/work').exists) dbFS.mkdir('/work');
  dbFS.chdir('/work');
  dbFS.writeFile('ref.fa', refText);

  const safe = $('#safeMode')?.checked;
  const dbArgs = safe ? SAFE_DB_ARGS.slice() : splitArgs($('#lastdbArgs').value);
  const lastdbArgv = ['-v', '-P', '1', ...dbArgs, dbName, 'ref.fa'];
  append(logEl, '$ lastdb ' + lastdbArgv.map(a => /\s/.test(a)? ('"'+a+'"'):a).join(' '));
  await tick();
  try {
    lastdbModule.callMain(lastdbArgv);
  } catch (e) {
    const status = (typeof e?.status === 'number') ? e.status : undefined;
    if (!(e?.name === 'ExitStatus' && status === 0)) {
      append(logEl, `lastdb 実行エラー: ${e?.message || e}`);
      if (typeof status === 'number') append(logEl, `exit code: ${status}`);
      setStatus('エラー');
      return;
    }
  }

  const list = dbFS.readdir('.')
    .filter(n => n !== '.' && n !== '..')
    .map(n => `${n} (${dbFS.isDir(dbFS.lookupPath(n).node.mode) ? 'dir' : 'file'})`);
  append(logEl, '現在のディレクトリ: /work');
  append(logEl, '一覧: ' + (list.join(', ') || '(なし)'));
  const entries = dbFS.readdir('.')
    .filter(n => n !== '.' && n !== '..' && (n === dbName || n.startsWith(dbName + '.')));
  append(logEl, '生成ファイル: ' + (entries.join(', ') || '(なし)'));
  setStatus('lastdb完了');
}

$('#runDbBtn').addEventListener('click', () => {
  runDbOnly().catch(err => {
    setStatus('エラー');
    append(logEl, String(err?.stack || err));
  });
});

// Toggle: disable manual args in safe mode
const safeCb = document.getElementById('safeMode');
if (safeCb) {
  const lastdbArgsInput = document.getElementById('lastdbArgs');
  safeCb.addEventListener('change', () => {
    if (!lastdbArgsInput) return;
    if (safeCb.checked) {
      lastdbArgsInput.value = SAFE_DB_ARGS.join(' ');
      lastdbArgsInput.setAttribute('disabled', 'disabled');
    } else {
      lastdbArgsInput.removeAttribute('disabled');
    }
  });
  // Apply initial state on load
  if (safeCb.checked && lastdbArgsInput) {
    lastdbArgsInput.value = SAFE_DB_ARGS.join(' ');
    lastdbArgsInput.setAttribute('disabled', 'disabled');
  }
}
