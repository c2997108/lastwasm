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
    let pathname = req.url.split('?')[0];
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
      res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
      res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
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

    await page.waitForFunction(() => document.querySelector('#status')?.textContent === 'Done', null, { timeout: 300000 });
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

    if (completeTab !== expectedTab) {
      throw new Error(`downloaded TAB differs from WASM baseline: expected ${expectedTab.length} chars, got ${completeTab.length}`);
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
    if (!/lastal-pthreads spawned=7 maxRunning=7/.test(log)) {
      throw new Error(`pthread completion log not found:\n${log}`);
    }
    if (!/dotplot-parser workers=8 alignments=\d+ plotted=\d+/.test(log)) {
      throw new Error(`parallel dot-plot parser log not found:\n${log}`);
    }
    if (!/dotplot-timing total=\d+\.\dms parser-total=\d+\.\dms plotly-newPlot=\d+\.\dms/.test(log)) {
      throw new Error(`detailed dot-plot timing log not found:\n${log}`);
    }
    if (!statusHistory.some(status => /^Searching\.\.\. (?:[1-9]|[1-9]\d)%$/.test(status))) {
      throw new Error(`intermediate pthread progress was not displayed: ${statusHistory.join(' | ')}`);
    }
    const resultReadyIndex = statusHistory.findIndex(status => status.startsWith('Result ready:'));
    const renderingIndex = statusHistory.indexOf('Rendering dot plot...');
    if (resultReadyIndex < 0 || renderingIndex <= resultReadyIndex || resultReadyTabLength < 1) {
      throw new Error(`TAB was not made ready before dot plot rendering: ${statusHistory.join(' | ')}`);
    }
    if (plotCanvasCount < 1) throw new Error('Plotly WebGL canvas was not created');
    if (wasmRequests.length < 1) throw new Error(`expected worker WASM requests, got: ${wasmRequests.join(', ')}`);

    const singleFasta = path.join(rootDir, 'test', 'huma.fa');
    const singleExpectedTab = buildWasmBaseline(rootDir, singleFasta, singleFasta);
    await page.setInputFiles('#refFasta', singleFasta);
    await page.setInputFiles('#qryFasta', singleFasta);
    await page.fill('#workerCount', '8');
    await page.click('#runBtn');
    await page.waitForFunction(() => document.querySelector('#status')?.textContent === 'Done', null, { timeout: 300000 });

    const singleTab = await page.textContent('#tab');
    const singleLog = await page.textContent('#log');
    const singleExpectedPreview = tabPreview(singleExpectedTab);
    if (singleTab !== singleExpectedPreview) {
      throw new Error(`single-record TAB preview differs from WASM baseline: expected ${singleExpectedPreview.length} chars, got ${singleTab.length}`);
    }
    if (!/runtime=compiled LAST wasm-pthreads threads=1/.test(singleLog)) {
      throw new Error(`single-record WASM runtime log not found:\n${singleLog}`);
    }
    if (!/Search threads were limited to 1 query FASTA record/.test(singleLog)) {
      throw new Error(`single-record thread limit warning not found:\n${singleLog}`);
    }

    await testLargeResultPipeline(context, rootDir, port);

    const largeFasta = process.env.LAST_WEB_LARGE_FASTA;
    if (largeFasta) {
      const largePath = path.resolve(largeFasta);
      await page.setInputFiles('#refFasta', largePath);
      await page.setInputFiles('#qryFasta', largePath);
      await page.fill('#workerCount', '8');
      await page.click('#runBtn');
      await page.waitForFunction(() => {
        const status = document.querySelector('#status')?.textContent;
        return status === 'Done' || status === 'Error';
      }, null, { timeout: 900000 });

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
      if (!String(largeLog || '').includes(`lastal-pthreads spawned=${childThreads} maxRunning=${childThreads}`)) {
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

async function testLargeResultPipeline(context, rootDir, port) {
  const alignmentCount = Math.max(1, Number(process.env.LAST_WEB_SYNTHETIC_ALIGNMENTS) || 100001);
  const sequenceCount = Math.max(1, Number(process.env.LAST_WEB_SYNTHETIC_SEQUENCES) || 1);
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
        setTimeout(() => {
          try {
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
                refs: [{ name: 'ref', length: 1000 }],
                queries: [{ name: 'qry', length: 1000 }],
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
                refs: [{ name: 'ref', length: 1000 }],
                queries: [{ name: 'qry', length: 1000 }],
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
  const fasta = path.join(rootDir, 'test', 'huma.fa');
  await page.setInputFiles('#refFasta', fasta);
  await page.setInputFiles('#qryFasta', fasta);
  await page.fill('#workerCount', '8');
  await page.click('#runBtn');
  await page.waitForFunction(() => document.querySelector('#status')?.textContent === 'Done', null, { timeout: 300000 });

  const state = await page.evaluate(() => ({
    tabText: document.querySelector('#tab')?.textContent || '',
    tabSummary: document.querySelector('#tabSummary')?.textContent || '',
    downloadEnabled: !document.querySelector('#downloadTabBtn')?.disabled,
    log: document.querySelector('#log')?.textContent || '',
    previewPaintAckLength: window.__largePreviewPaintAckLength || 0,
  }));
  const plotCanvasCount = await page.locator('#dotplot canvas').count();
  if (!state.tabText.includes('Preview limited to the first 20 KiB') || state.tabText.length > 22 * 1024) {
    throw new Error(`large TAB preview was not limited to 20 KiB: ${state.tabText.length}`);
  }
  if (state.previewPaintAckLength < 1) {
    throw new Error('full TAB processing started before the preview paint acknowledgement');
  }
  if (!/Complete TAB · \d+\.\d MiB/.test(state.tabSummary) || !state.downloadEnabled) {
    throw new Error(`large TAB download was not prepared: ${JSON.stringify(state)}`);
  }
  const sampleEvery = Math.max(1, Math.ceil(alignmentCount / 100000));
  const plottedCount = Math.ceil(alignmentCount / sampleEvery);
  if (!state.log.includes(`dotplot-parser workers=8 alignments=${alignmentCount} plotted=${plottedCount}`)) {
    throw new Error(`large dot plot was not sampled as expected:\n${state.log}`);
  }
  const profile = state.log.split('\n')
    .filter(line => line.startsWith('dotplot-timing') || line.startsWith('dotplot-longtasks') || line.startsWith('dotplot-memory'));
  console.log(`Synthetic ${alignmentCount}-alignment/${sequenceCount}-sequence profile:\n${profile.join('\n')}`);
  if (plotCanvasCount < 1) throw new Error('large sampled Plotly canvas was not created');
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
