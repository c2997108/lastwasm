// Headless smoke test for the browser page, including the threaded runtime.

const { chromium } = require('playwright');
const path = require('path');
const os = require('os');
const fs = require('fs');
const { spawnSync } = require('child_process');
const TAB_PREVIEW_CHARS = 20 * 1024;

async function main() {
  const rootDir = path.resolve(__dirname, '..');
  const http = require('http');
  const mime = {
    '.html': 'text/html',
    '.js': 'application/javascript',
    '.mjs': 'application/javascript',
    '.css': 'text/css',
    '.wasm': 'application/wasm',
    '.json': 'application/json',
    '.txt': 'text/plain',
    '.tab': 'text/plain',
    '.fa': 'text/plain',
    '.fasta': 'text/plain',
  };

  const server = http.createServer((req, res) => {
    const requestUrl = new URL(req.url, 'http://127.0.0.1');
    let pathname = requestUrl.pathname;
    const unisolated = requestUrl.searchParams.has('unisolated');
    if (pathname === '/') pathname = '/index.html';
    const file = path.join(rootDir, pathname);
    if (!file.startsWith(rootDir)) {
      res.writeHead(403);
      res.end('Forbidden');
      return;
    }
    fs.readFile(file, (err, data) => {
      if (err) {
        res.writeHead(404);
        res.end('Not found');
        return;
      }
      if (!unisolated) {
        res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
        res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
      }
      res.setHeader('Content-Type', mime[path.extname(file)] || 'application/octet-stream');
      res.end(data);
    });
  });

  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const port = server.address().port;

  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext();
  const page = await context.newPage();
  const wasmRequests = [];

  page.on('request', request => {
    if (request.url().includes('.wasm')) wasmRequests.push(request.url());
  });
  page.on('console', message => console.log('[browser]', message.type(), message.text()));
  page.on('pageerror', error => console.error('[browser]', error.stack || error));

  let queryFasta = '';
  try {
    const refFasta = path.join(rootDir, 'test', 'galGal3-M-32.fa');
    queryFasta = makeQueryFasta(rootDir);
    const expectedTab = buildWasmBaseline(rootDir, refFasta, queryFasta);

    await page.goto(`http://127.0.0.1:${port}/`, { waitUntil: 'domcontentloaded' });
    const isolated = await page.evaluate(() => Boolean(crossOriginIsolated));
    if (!isolated) throw new Error('test page is not cross-origin isolated');
    await page.evaluate(() => {
      window.__statusHistory = [];
      window.__resultReadyTabLength = 0;
      const status = document.querySelector('#status');
      const tab = document.querySelector('#tab');
      new MutationObserver(() => {
        const value = status?.textContent || '';
        window.__statusHistory.push(value);
        if (value.startsWith('Result ready:')) {
          window.__resultReadyTabLength = tab?.textContent.length || 0;
        }
      })
        .observe(status, { childList: true, characterData: true, subtree: true });
    });
    await page.setInputFiles('#refFasta', refFasta);
    await page.setInputFiles('#qryFasta', queryFasta);
    await page.fill('#workerCount', '8');
    await page.click('#runBtn');

    await page.waitForFunction(() => ['Done', 'Error'].includes(document.querySelector('#status')?.textContent), null, { timeout: 300000 });
    if (await page.textContent('#status') === 'Error') {
      throw new Error(`browser run failed:\n${await page.textContent('#log')}`);
    }
    const tab = await page.textContent('#tab');
    const log = await page.textContent('#log');
    const timeEstimate = await page.textContent('#timeEstimate');
    const statusHistory = await page.evaluate(() => window.__statusHistory || []);
    const resultReadyTabLength = await page.evaluate(() => window.__resultReadyTabLength || 0);
    const plotCanvasCount = await page.locator('#dotplot canvas').count();
    const downloadPromise = page.waitForEvent('download');
    await page.click('#downloadTabBtn');
    const download = await downloadPromise;
    const completeTab = fs.readFileSync(await download.path(), 'utf8');

    if (comparableTab(completeTab) !== comparableTab(expectedTab)) {
      const mismatch = firstMismatch(comparableTab(completeTab), comparableTab(expectedTab));
      throw new Error(`downloaded TAB differs from WASM baseline at ${mismatch}: expected ${expectedTab.length} chars, got ${completeTab.length}\n` +
        `expected=${JSON.stringify(expectedTab.slice(mismatch, mismatch + 160))}\n` +
        `actual=${JSON.stringify(completeTab.slice(mismatch, mismatch + 160))}`);
    }
    const expectedPreview = tabPreview(completeTab);
    if (tab !== expectedPreview) {
      const mismatch = firstMismatch(tab, expectedPreview);
      throw new Error(`TAB preview differs from complete download at ${mismatch}: expected ${expectedPreview.length} chars, got ${tab.length}\nexpected=${JSON.stringify(expectedPreview.slice(mismatch, mismatch + 120))}\nactual=${JSON.stringify(tab.slice(mismatch, mismatch + 120))}`);
    }
    if (!/Elapsed \d+:\d{2}/.test(timeEstimate || '')) {
      throw new Error(`elapsed time display not found: ${timeEstimate}`);
    }
    if (!log || !/runtime=compiled LAST wasm-pthreads threads=8/.test(log)) {
      throw new Error(`pthread runtime log not found:\n${log}`);
    }
    const pthreadStats = log.match(/lastal-pthreads spawned=(\d+) maxRunning=(\d+)/);
    if (!pthreadStats || Number(pthreadStats[1]) < 7 || Number(pthreadStats[1]) % 7 !== 0 || Number(pthreadStats[2]) !== 7) {
      throw new Error(`pthread completion log not found:\n${log}`);
    }
    if (!/dotplot-parser workers=8 alignments=\d+ rendered=\d+/.test(log)) {
      throw new Error(`parallel dot-plot parser log not found:\n${log}`);
    }
    if (!/dotplot-timing total=\d+\.\dms parser-total=\d+\.\dms renderer-total=\d+\.\dms/.test(log)) {
      throw new Error(`detailed dot-plot timing log not found:\n${log}`);
    }
    if (!/result-timing lastdb=\d+\.\dms lastal-compute=\d+\.\dms tab-preview=\d+\.\dms tab-join=\d+\.\dms utf8-encode=\d+\.\dms transfer=\d+\.\dms tab=\d+\.\dMiB/.test(log)) {
      throw new Error(`result-pipeline timing log not found:\n${log}`);
    }
    if (!statusHistory.includes('Searching... 0%') || !statusHistory.includes('Searching... 100%')) {
      throw new Error(`pthread progress start/completion was not displayed: ${statusHistory.join(' | ')}`);
    }
    const resultReadyIndex = statusHistory.findIndex(status => status.startsWith('Result ready:'));
    const renderingIndex = statusHistory.indexOf('Rendering dot plot...');
    if (resultReadyIndex < 0 || renderingIndex <= resultReadyIndex || resultReadyTabLength < 1) {
      throw new Error(`TAB was not made ready before dot plot rendering: ${statusHistory.join(' | ')}`);
    }
    if (plotCanvasCount < 2 || await page.locator('.dotplot-gl').count() !== 1 || await page.locator('.dotplot-axis').count() !== 1) {
      throw new Error('regl WebGL and axis canvases were not created');
    }
    if (wasmRequests.length < 1 || wasmRequests.some(url => !url.includes('v=20260725-axis-labels'))) {
      throw new Error(`expected versioned worker WASM requests, got: ${wasmRequests.join(', ')}`);
    }

    const singleFasta = path.join(rootDir, 'test', 'huma.fa');
    const singleExpectedTab = buildWasmBaseline(rootDir, singleFasta, singleFasta);
    await page.setInputFiles('#refFasta', singleFasta);
    await page.setInputFiles('#qryFasta', singleFasta);
    await page.fill('#workerCount', '8');
    await page.click('#runBtn');
    await page.waitForFunction(() => ['Done', 'Error'].includes(document.querySelector('#status')?.textContent), null, { timeout: 300000 });
    if (await page.textContent('#status') === 'Error') {
      throw new Error(`single-record browser run failed:\n${await page.textContent('#log')}`);
    }

    const singleTab = await page.textContent('#tab');
    const singleLog = await page.textContent('#log');
    const singleExpectedPreview = tabPreview(singleExpectedTab);
    if (comparableTab(singleTab) !== comparableTab(singleExpectedPreview)) {
      throw new Error(`single-record TAB preview differs from WASM baseline: expected ${singleExpectedPreview.length} chars, got ${singleTab.length}`);
    }
    if (!/runtime=compiled LAST wasm-pthreads threads=1/.test(singleLog)) {
      throw new Error(`single-record WASM runtime log not found:\n${singleLog}`);
    }
    if (!/Search threads were limited to 1 query FASTA record/.test(singleLog)) {
      throw new Error(`single-record thread limit warning not found:\n${singleLog}`);
    }

    await testServiceWorkerIsolation(browser, rootDir, port, singleExpectedTab);
    await testAsmFallback(browser, rootDir, port, singleExpectedTab);
    await testLargeResultPipeline(context, rootDir, port);

    const largeFasta = process.env.LAST_WEB_LARGE_FASTA;
    if (largeFasta) {
      const largePath = path.resolve(largeFasta);
      await page.setInputFiles('#refFasta', largePath);
      await page.setInputFiles('#qryFasta', largePath);
      await page.fill('#workerCount', '8');
      await page.click('#runBtn');
      const profileIntervalMs = Math.max(0, Number(process.env.LAST_WEB_PROFILE_INTERVAL_MS) || 0);
      const profileStartedAt = Date.now();
      let profileBusy = false;
      const profileTimer = profileIntervalMs > 0 ? setInterval(async () => {
        if (profileBusy || page.isClosed()) return;
        profileBusy = true;
        try {
          const snapshot = await page.evaluate(() => {
            const logLines = (document.querySelector('#log')?.textContent || '').trim().split('\n');
            return {
              status: document.querySelector('#status')?.textContent || '',
              estimate: document.querySelector('#timeEstimate')?.textContent || '',
              tabSummary: document.querySelector('#tabSummary')?.textContent || '',
              tabChars: document.querySelector('#tab')?.textContent.length || 0,
              canvases: document.querySelectorAll('#dotplot canvas').length,
              heapMiB: performance.memory ? performance.memory.usedJSHeapSize / (1024 * 1024) : null,
              logTail: logLines.slice(-2),
            };
          });
          const profileLine = `Large FASTA monitor t=${Math.round((Date.now() - profileStartedAt) / 1000)}s ${JSON.stringify(snapshot)}`;
          console.log(profileLine);
          if (process.env.LAST_WEB_PROFILE_LOG) {
            fs.appendFileSync(process.env.LAST_WEB_PROFILE_LOG, `${profileLine}\n`);
          }
        } catch (error) {
          console.log(`Large FASTA monitor unavailable: ${error.message}`);
        } finally {
          profileBusy = false;
        }
      }, profileIntervalMs) : null;
      try {
        await page.waitForFunction(() => {
          const status = document.querySelector('#status')?.textContent;
          return status === 'Done' || status === 'Error';
        }, null, { timeout: 900000 });
      } finally {
        if (profileTimer) clearInterval(profileTimer);
      }

      const largeStatus = await page.textContent('#status');
      const largeLog = await page.textContent('#log');
      const largeElapsed = await page.textContent('#timeEstimate');
      const expectedThreads = Math.max(1, Math.min(8, fastaRecordCount(largePath)));
      const largeProfile = String(largeLog || '').split('\n')
        .filter(line => line.startsWith('dotplot-timing') || line.startsWith('dotplot-longtasks') || line.startsWith('dotplot-memory'));
      console.log(`Large FASTA profile:\n${largeProfile.join('\n')}`);
      if (largeStatus !== 'Done') {
        throw new Error(`large FASTA browser run failed:\n${largeLog}`);
      }
      if (!String(largeLog || '').includes(`runtime=compiled LAST wasm-pthreads threads=${expectedThreads}`)) {
        throw new Error(`large FASTA pthread runtime log not found:\n${largeLog}`);
      }
      const childThreads = expectedThreads - 1;
      const pthreadStats = String(largeLog || '').match(/lastal-pthreads spawned=(\d+) maxRunning=(\d+)/);
      const spawned = Number(pthreadStats?.[1]);
      const maxRunning = Number(pthreadStats?.[2]);
      if (!pthreadStats || spawned < childThreads || spawned % childThreads !== 0 || maxRunning !== childThreads) {
        throw new Error(`large FASTA pthread completion log not found:\n${largeLog}`);
      }
      console.log(`Large FASTA E2E OK: ${largeElapsed}`);
    }
  } finally {
    await browser.close();
    server.close();
    if (queryFasta) fs.rmSync(path.dirname(queryFasta), { recursive: true, force: true });
  }

  console.log('E2E OK');
}

async function testServiceWorkerIsolation(browser, rootDir, port, expectedTab) {
  const context = await browser.newContext();
  const page = await context.newPage();
  page.on('pageerror', error => console.error('[service-worker browser]', error.stack || error));

  try {
    await page.goto(`http://127.0.0.1:${port}/?unisolated=1`, { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => (
      crossOriginIsolated
      && typeof SharedArrayBuffer !== 'undefined'
      && Boolean(navigator.serviceWorker.controller)
    ), null, { timeout: 30000 });
    await page.waitForSelector('#runBtn');

    const fasta = path.join(rootDir, 'test', 'huma.fa');
    await page.setInputFiles('#refFasta', fasta);
    await page.setInputFiles('#qryFasta', fasta);
    await page.fill('#workerCount', '8');
    await page.click('#runBtn');
    await page.waitForFunction(() => ['Done', 'Error'].includes(document.querySelector('#status')?.textContent), null, { timeout: 300000 });

    const status = await page.textContent('#status');
    const log = await page.textContent('#log');
    const tab = await page.textContent('#tab');
    if (status !== 'Done') throw new Error(`Service Worker pthread run failed:\n${log}`);
    if (!/runtime=compiled LAST wasm-pthreads threads=1/.test(log || '')) {
      throw new Error(`Service Worker did not enable the pthread runtime:\n${log}`);
    }
    if (comparableTab(tab) !== comparableTab(tabPreview(expectedTab))) {
      throw new Error(`Service Worker pthread TAB differs from baseline: expected ${expectedTab.length} chars, got ${tab.length}`);
    }
  } finally {
    await context.close();
  }
}

async function testAsmFallback(browser, rootDir, port, expectedTab) {
  const context = await browser.newContext({ serviceWorkers: 'block' });
  const page = await context.newPage();
  const runtimeRequests = [];
  page.on('request', request => {
    if (/lastal(?:[-.]asm)?\.(?:js|wasm|mem)/.test(request.url())) runtimeRequests.push(request.url());
  });
  page.on('pageerror', error => console.error('[asm.js browser]', error.stack || error));

  try {
    await page.goto(`http://127.0.0.1:${port}/?unisolated=1`, { waitUntil: 'domcontentloaded' });
    if (await page.evaluate(() => Boolean(crossOriginIsolated) || typeof SharedArrayBuffer !== 'undefined')) {
      throw new Error('asm.js fallback test unexpectedly has shared-memory support');
    }

    const fasta = path.join(rootDir, 'test', 'huma.fa');
    await page.setInputFiles('#refFasta', fasta);
    await page.setInputFiles('#qryFasta', fasta);
    await page.fill('#workerCount', '8');
    await page.click('#runBtn');
    await page.waitForFunction(() => ['Done', 'Error'].includes(document.querySelector('#status')?.textContent), null, { timeout: 300000 });

    const status = await page.textContent('#status');
    const log = await page.textContent('#log');
    const tab = await page.textContent('#tab');
    if (status !== 'Done') throw new Error(`asm.js fallback run failed:\n${log}`);
    if (!/runtime=compiled LAST asmjs threads=1/.test(log || '')) {
      throw new Error(`asm.js fallback runtime log not found:\n${log}`);
    }
    if (comparableTab(tab) !== comparableTab(tabPreview(expectedTab))) {
      throw new Error(`asm.js fallback TAB differs from baseline: expected ${expectedTab.length} chars, got ${tab.length}`);
    }
    if (!runtimeRequests.some(url => url.includes('/web/lastal-asm.js') && url.includes('v=20260725-axis-labels'))) {
      throw new Error(`versioned asm.js runtime was not loaded: ${runtimeRequests.join(', ')}`);
    }
    if (runtimeRequests.some(url => url.includes('/web/lastal.wasm'))) {
      throw new Error(`pthread WASM was loaded by the asm.js fallback: ${runtimeRequests.join(', ')}`);
    }
  } finally {
    await context.close();
  }
}

async function testLargeResultPipeline(context, rootDir, port) {
  const alignmentCount = Math.max(1, Number(process.env.LAST_WEB_SYNTHETIC_ALIGNMENTS) || 100001);
  const sequenceCount = Math.max(1, Number(process.env.LAST_WEB_SYNTHETIC_SEQUENCES) || 3);
  const page = await context.newPage();
  await page.addInitScript(({ alignmentCount, previewChars, sequenceCount }) => {
    const NativeWorker = window.Worker;
    window.Worker = class WorkerProxy {
      constructor(url, options) {
        if (!String(url).includes('jslast-runner-worker.js')) return new NativeWorker(url, options);
        this.onmessage = null;
        this.onerror = null;
        this.pendingDone = null;
      }

      postMessage(message) {
        if (message?.type === 'tab-preview-painted') {
          window.__largePreviewPaintAckLength = document.querySelector('#tab')?.textContent.length || 0;
          const done = this.pendingDone;
          this.pendingDone = null;
          setTimeout(() => this.emit(done), 0);
          return;
        }
        if (message?.type !== 'run') return;
        window.__largeRunOptions = {
          lastdbArgs: message.payload?.lastdbArgs,
          lastalArgs: message.payload?.lastalArgs,
          useDefaultDbArgs: message.payload?.useDefaultDbArgs,
        };
        setTimeout(() => {
          try {
            const fastaRefs = Array.from({ length: sequenceCount }, (_, index) => ({
              name: `ref_${sequenceCount - index - 1}`,
              length: 1000,
            }));
            const fastaQueries = Array.from({ length: sequenceCount }, (_, index) => ({
              name: `qry_${sequenceCount - index - 1}`,
              length: 1000,
            }));
            const lines = Array.from({ length: sequenceCount }, (_, index) => (
              `100\tref_${index}\t0\t10\t+\t1000\tqry_${index}\t0\t10\t+\t1000\t10\tEG2=0\tE=0\n`
            ));
            const block = lines.join('');
            const repeats = Math.floor(alignmentCount / sequenceCount);
            const remainder = alignmentCount % sequenceCount;
            const tabText = `# synthetic large result\n${block.repeat(repeats)}${lines.slice(0, remainder).join('')}# Query sequences=${sequenceCount} normal letters=${sequenceCount * 1000}\n`;
            const previewEnd = tabText.lastIndexOf('\n', previewChars) + 1;
            this.emit({
              type: 'progress',
              event: {
                stage: 'configured',
                message: 'Synthetic large result',
                refs: fastaRefs,
                queries: fastaQueries,
                runtime: 'wasm-pthreads',
                threads: 8,
                warnings: [],
              },
            });
            this.emit({ type: 'progress', event: { stage: 'done', message: `Done: ${alignmentCount} alignments` } });
            this.emit({
              type: 'tab-preview',
              text: tabText.slice(0, previewEnd),
              truncated: true,
              totalChars: tabText.length,
            });
            const tabBytes = new TextEncoder().encode(tabText);
            this.pendingDone = {
              type: 'done',
              result: {
                alignmentCount,
                refs: fastaRefs,
                queries: fastaQueries,
                threads: 8,
                runtime: 'wasm-pthreads',
                warnings: [],
                dbThreadStats: { spawned: 0, maxRunning: 0 },
                alThreadStats: { spawned: 7, maxRunning: 7 },
                searchWorkerStats: { started: 0, completed: 0 },
                searchProgressStats: null,
                timing: {
                  totalMs: 1000,
                  dbMs: 100,
                  searchMs: 900,
                  refLetters: 1000,
                  queryLetters: 1000,
                  searchWorkers: 8,
                },
                dbLog: [],
                alLog: [],
                tabByteLength: tabBytes.byteLength,
              },
              tabBuffer: tabBytes.buffer,
            };
          } catch (error) {
            this.onerror?.({ message: error?.message || String(error) });
          }
        }, 0);
      }

      emit(data) {
        this.onmessage?.({ data });
      }

      terminate() {}
    };
  }, { alignmentCount, previewChars: TAB_PREVIEW_CHARS, sequenceCount });

  page.on('pageerror', error => console.error('[large-result browser]', error.stack || error));
  await page.goto(`http://127.0.0.1:${port}/`, { waitUntil: 'domcontentloaded' });
  const formattedAxisLabels = await page.evaluate(async () => {
    const { formatAxisNumbers } = await import('/web/regl-dotplot.js?v=axis-label-test');
    return {
      twentyKilobase: formatAxisNumbers([1320000, 1340000, 1360000]),
      fiveKilobase: formatAxisNumbers([1400000, 1405000, 1410000]),
      oneBaseAtGigabase: formatAxisNumbers([1400000000, 1400000001]),
    };
  });
  if (formattedAxisLabels.twentyKilobase.join(',') !== '1.32M,1.34M,1.36M'
      || formattedAxisLabels.fiveKilobase.join(',') !== '1.400M,1.405M,1.410M'
      || formattedAxisLabels.oneBaseAtGigabase.join(',') !== '1.400000000G,1.400000001G') {
    throw new Error(`zoomed axis labels lost precision: ${JSON.stringify(formattedAxisLabels)}`);
  }
  if (await page.inputValue('#plotLengthCutoff') !== '1000') {
    throw new Error('dot-plot hidden alignment length does not default to 1000 bp');
  }
  if (await page.inputValue('#seedHitLimit') !== '100'
      || await page.inputValue('#queryBatchSize') !== '64M'
      || !await page.isChecked('#useFastaOrder')
      || await page.locator('#lastdbArgs').count() !== 0
      || await page.locator('#lastalArgs').count() !== 0
      || await page.locator('#safeMode').count() !== 0) {
    throw new Error('simplified LAST options are not configured correctly');
  }
  const fasta = path.join(rootDir, 'test', 'huma.fa');
  await page.setInputFiles('#refFasta', fasta);
  await page.setInputFiles('#qryFasta', fasta);
  await page.fill('#workerCount', '8');
  await page.fill('#queryBatchSize', '32m');
  await page.fill('#plotLengthCutoff', '9');
  await page.click('#runBtn');
  await page.waitForFunction(() => ['Done', 'Error'].includes(document.querySelector('#status')?.textContent), null, { timeout: 300000 });
  if (await page.textContent('#status') === 'Error') {
    throw new Error(`large-result browser run failed:\n${await page.textContent('#log')}`);
  }

  const state = await page.evaluate(() => ({
    tabText: document.querySelector('#tab')?.textContent || '',
    tabSummary: document.querySelector('#tabSummary')?.textContent || '',
    downloadEnabled: !document.querySelector('#downloadTabBtn')?.disabled,
    log: document.querySelector('#log')?.textContent || '',
    previewPaintAckLength: window.__largePreviewPaintAckLength || 0,
    runOptions: window.__largeRunOptions,
  }));
  const plotCanvasCount = await page.locator('#dotplot canvas').count();
  if (!state.tabText.includes('Preview limited to the first 20 KiB') || state.tabText.length > 22 * 1024) {
    throw new Error(`large TAB preview was not limited to 20 KiB: ${state.tabText.length}`);
  }
  if (state.previewPaintAckLength < 1) {
    throw new Error('full TAB processing started before the preview paint acknowledgement');
  }
  if (state.runOptions?.lastdbArgs !== ''
      || state.runOptions?.lastalArgs !== '-m100 -i32M'
      || state.runOptions?.useDefaultDbArgs !== true) {
    throw new Error(`simplified LAST options produced the wrong arguments: ${JSON.stringify(state.runOptions)}`);
  }
  if (!state.log.includes('dotplot sequence-order=fasta')) {
    throw new Error(`FASTA sequence order was not enabled:\n${state.log}`);
  }
  if (!/Complete TAB · \d+(?:\.\d+)? (?:KiB|MiB)/.test(state.tabSummary) || !state.downloadEnabled) {
    throw new Error(`large TAB download was not prepared: ${JSON.stringify(state)}`);
  }
  if (!state.log.includes(`dotplot-parser workers=8 alignments=${alignmentCount} rendered=${alignmentCount} hide-length<=9bp`)) {
    throw new Error(`large dot plot did not render every alignment:\n${state.log}`);
  }
  const profile = state.log.split('\n')
    .filter(line => line.startsWith('dotplot-timing') || line.startsWith('dotplot-longtasks') || line.startsWith('dotplot-memory'));
  console.log(`Synthetic ${alignmentCount}-alignment/${sequenceCount}-sequence profile:\n${profile.join('\n')}`);
  if (plotCanvasCount < 2) throw new Error('large full-resolution regl canvases were not created');
  await page.locator('.dotplot-stage').scrollIntoViewIfNeeded();
  const stage = await page.locator('.dotplot-stage').boundingBox();
  if (!stage) throw new Error('dot-plot interaction stage was not created');
  if (alignmentCount <= 250000) {
    const plotWidth = stage.width - 68 - 18;
    const plotHeight = stage.height - 22 - 54;
    const scale = Math.min(plotWidth / (sequenceCount * 1000), plotHeight / (sequenceCount * 1000)) * 0.96;
    const targetCoordinate = (sequenceCount - 1) * 1000 + 5;
    const targetX = stage.x + 68 + plotWidth / 2 + (targetCoordinate - sequenceCount * 500) * scale;
    const targetY = stage.y + 22 + plotHeight / 2 - (targetCoordinate - sequenceCount * 500) * scale;
    await page.mouse.move(targetX, targetY);
    await page.mouse.wheel(0, -1000);
    await page.waitForTimeout(100);
    await page.mouse.move(targetX + 1, targetY + 1);
    await page.waitForTimeout(1000);
    if (process.env.LAST_WEB_DOTPLOT_SCREENSHOT) {
      await page.locator('#dotplot').screenshot({ path: process.env.LAST_WEB_DOTPLOT_SCREENSHOT });
    }
    const hoverState = await page.evaluate(() => {
      const tooltip = document.querySelector('.dotplot-tooltip');
      return {
        hidden: tooltip?.hidden,
        text: tooltip?.textContent,
      };
    });
    if (hoverState.hidden
        || !hoverState.text?.includes('Score: 100')
        || !hoverState.text?.includes('Ref: ref_0:1-10')
        || !hoverState.text?.includes('Qry: qry_0:1-10')) {
      throw new Error(`GPU picking did not return the alignment: ${JSON.stringify({ targetX, targetY, hoverState })}`);
    }
  }
  if (alignmentCount <= 250000 && await page.locator('.dotplot-stage').getAttribute('data-view-mode') !== 'custom') {
    throw new Error('wheel zoom did not update the dot-plot view');
  }
  await page.click('.dotplot-fit');
  await page.waitForTimeout(100);
  if (await page.locator('.dotplot-stage').getAttribute('data-view-mode') !== 'fit') {
    throw new Error('Fit button did not reset the dot-plot view');
  }
  await page.setViewportSize({ width: 390, height: 844 });
  await page.locator('.dotplot-stage').scrollIntoViewIfNeeded();
  await page.waitForTimeout(100);
  const mobilePlot = await page.locator('#dotplot').boundingBox();
  const mobileCanvases = await page.locator('#dotplot canvas').evaluateAll(canvases => canvases.map(canvas => ({
    width: canvas.width,
    height: canvas.height,
    cssWidth: canvas.clientWidth,
    cssHeight: canvas.clientHeight,
  })));
  if (!mobilePlot || mobilePlot.width > 390 || mobileCanvases.length !== 2
      || mobileCanvases.some(canvas => canvas.width < 300 || canvas.height < 300 || canvas.cssWidth < 300 || canvas.cssHeight < 300)) {
    throw new Error(`mobile dot plot did not resize correctly: ${JSON.stringify({ mobilePlot, mobileCanvases })}`);
  }
  if (process.env.LAST_WEB_DOTPLOT_MOBILE_SCREENSHOT) {
    await page.locator('#dotplot').screenshot({ path: process.env.LAST_WEB_DOTPLOT_MOBILE_SCREENSHOT });
  }

  await page.fill('#plotLengthCutoff', '10');
  await page.click('#runBtn');
  await page.waitForFunction(() => ['Done', 'Error'].includes(document.querySelector('#status')?.textContent), null, { timeout: 300000 });
  const filteredStatus = await page.textContent('#status');
  const filteredLog = await page.textContent('#log');
  const filteredPlot = await page.textContent('#dotplot');
  if (filteredStatus !== 'Done'
      || !filteredLog.includes(`dotplot-parser workers=8 alignments=${alignmentCount} rendered=0 hide-length<=10bp`)
      || !filteredPlot.includes('No alignments longer than 10 bp.')
      || await page.locator('#dotplot canvas').count() !== 0) {
    throw new Error(`dot-plot minimum length filter failed:\n${filteredLog}`);
  }
  await page.close();
}

function makeQueryFasta(rootDir) {
  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'lastwasm-query-'));
  const queryPath = path.join(tmpDir, 'query.fa');
  const seed = fs.readFileSync(path.join(rootDir, 'test', 'galGal3-M-32.fa'), 'utf8');
  const firstRecord = seed.split(/^>/m).find(record => record.trim());
  const sequence = firstRecord.split(/\r?\n/).slice(1).join('').slice(0, 2000);
  const records = [];
  for (let i = 0; i < 512; i++) {
    records.push(`>query_${i + 1}\n${sequence}\n`);
  }
  fs.writeFileSync(queryPath, records.join(''));
  return queryPath;
}

function fastaRecordCount(filePath) {
  const matches = fs.readFileSync(filePath, 'utf8').match(/^>/gm);
  return Math.max(1, matches?.length || 0);
}

function tabPreview(tabText) {
  if (tabText.length <= TAB_PREVIEW_CHARS) return tabText;
  const newline = tabText.lastIndexOf('\n', TAB_PREVIEW_CHARS);
  const previewEnd = newline > 0 ? newline + 1 : TAB_PREVIEW_CHARS;
  return `${tabText.slice(0, previewEnd)}\n# Preview limited to the first 20 KiB. Download TAB contains the complete result.\n`;
}

function firstMismatch(left, right) {
  const length = Math.min(left.length, right.length);
  for (let index = 0; index < length; index++) {
    if (left[index] !== right[index]) return index;
  }
  return length;
}

function comparableTab(text) {
  return text.replace(/^# LAST version .*$/m, '# LAST version <build>');
}

function buildWasmBaseline(rootDir, refFastaPath, queryFastaPath) {
  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'lastwasm-baseline-'));
  const dbPath = path.join(tmpDir, 'refdb');
  try {
    runNode(rootDir, [
      'cli/lastdb.js',
      '-P', '1',
      '--bits=4',
      '-R00',
      '-uNEAR',
      '-w1',
      '-W1',
      '-S1',
      '-C1',
      '-v',
      dbPath,
      refFastaPath,
    ]);
    return runNode(rootDir, [
      'cli/lastal.js',
      '-f', 'TAB',
      '-P', '1',
      '-i', '64M',
      '-m100',
      dbPath,
      queryFastaPath,
    ]).stdout;
  } finally {
    fs.rmSync(tmpDir, { recursive: true, force: true });
  }
}

function runNode(rootDir, args, env = {}) {
  const result = spawnSync(process.execPath, args, {
    cwd: rootDir,
    encoding: 'utf8',
    env: { ...process.env, ...env },
    maxBuffer: 64 * 1024 * 1024,
  });
  if (result.status !== 0) {
    throw new Error(`${args.join(' ')} failed with ${result.status}\n${result.stderr || result.stdout}`);
  }
  return result;
}

main().catch(error => {
  console.error('E2E FAILED:', error && error.stack || error);
  process.exit(1);
});
