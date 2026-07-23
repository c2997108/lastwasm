// Headless smoke test for the browser page, including the threaded runtime.

const { chromium } = require('playwright');
const path = require('path');
const os = require('os');
const fs = require('fs');
const { spawnSync } = require('child_process');

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

    if (tab !== expectedTab) {
      throw new Error(`TAB output differs from WASM baseline: expected ${expectedTab.length} chars, got ${tab.length}`);
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
    if (singleTab !== singleExpectedTab) {
      throw new Error(`single-record TAB output differs from WASM baseline: expected ${singleExpectedTab.length} chars, got ${singleTab.length}`);
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
      if (largeStatus !== 'Done') {
        throw new Error(`large FASTA browser run failed:\n${largeLog}`);
      }
      if (!/runtime=compiled LAST wasm-pthreads threads=8/.test(largeLog || '')) {
        throw new Error(`large FASTA pthread runtime log not found:\n${largeLog}`);
      }
      if (!/lastal-pthreads spawned=7 maxRunning=7/.test(largeLog || '')) {
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
  const page = await context.newPage();
  await page.addInitScript(({ alignmentCount, previewChars }) => {
    const NativeWorker = window.Worker;
    window.Worker = class WorkerProxy {
      constructor(url, options) {
        if (!String(url).includes('jslast-runner-worker.js')) return new NativeWorker(url, options);
        this.onmessage = null;
        this.onerror = null;
      }

      postMessage() {
        setTimeout(() => {
          try {
            const line = '100\tref\t0\t10\t+\t1000\tqry\t0\t10\t+\t1000\t10\tEG2=0\tE=0\n';
            const tabText = `# synthetic large result\n${line.repeat(alignmentCount)}# Query sequences=1 normal letters=1000\n`;
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
            this.emit({
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
            });
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
  }, { alignmentCount: 100001, previewChars: 2 * 1024 * 1024 });

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
  }));
  const plotCanvasCount = await page.locator('#dotplot canvas').count();
  if (!state.tabText.includes('Preview limited to the first 2 MiB')) {
    throw new Error('large TAB preview was not limited');
  }
  if (!/Complete TAB · \d+\.\d MiB/.test(state.tabSummary) || !state.downloadEnabled) {
    throw new Error(`large TAB download was not prepared: ${JSON.stringify(state)}`);
  }
  if (!/dotplot-parser workers=8 alignments=100001 plotted=50001/.test(state.log)) {
    throw new Error(`large dot plot was not sampled as expected:\n${state.log}`);
  }
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
