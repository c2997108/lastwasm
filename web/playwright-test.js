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
    await page.setInputFiles('#refFasta', refFasta);
    await page.setInputFiles('#qryFasta', queryFasta);
    await page.fill('#workerCount', '8');
    await page.click('#runBtn');

    await page.waitForFunction(() => document.querySelector('#status')?.textContent === 'Done', null, { timeout: 300000 });
    const tab = await page.textContent('#tab');
    const log = await page.textContent('#log');
    const timeEstimate = await page.textContent('#timeEstimate');

    if (tab !== expectedTab) {
      throw new Error(`TAB output differs from WASM baseline: expected ${expectedTab.length} chars, got ${tab.length}`);
    }
    if (!/Elapsed \d+:\d{2}/.test(timeEstimate || '')) {
      throw new Error(`elapsed time display not found: ${timeEstimate}`);
    }
    if (!log || !/runtime=compiled LAST wasm-worker-pool threads=8/.test(log)) {
      throw new Error(`worker-pool runtime log not found:\n${log}`);
    }
    if (!/lastal-workers started=8 completed=8/.test(log)) {
      throw new Error(`worker-pool completion log not found:\n${log}`);
    }
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
    if (!/runtime=compiled LAST asmjs threads=1/.test(singleLog)) {
      throw new Error(`single-record fallback runtime log not found:\n${singleLog}`);
    }
    if (!/Search workers were limited to 1 query FASTA record/.test(singleLog)) {
      throw new Error(`single-record worker limit warning not found:\n${singleLog}`);
    }
  } finally {
    await browser.close();
    server.close();
    if (queryFasta) fs.rmSync(path.dirname(queryFasta), { recursive: true, force: true });
  }

  console.log('E2E OK');
}

function makeQueryFasta(rootDir) {
  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'lastwasm-query-'));
  const queryPath = path.join(tmpDir, 'query.fa');
  const seed = fs.readFileSync(path.join(rootDir, 'test', 'galGal3-M-32.fa'), 'utf8');
  const firstRecord = seed.split(/^>/m).find(record => record.trim());
  const sequence = firstRecord.split(/\r?\n/).slice(1).join('').slice(0, 2000);
  const records = [];
  for (let i = 0; i < 8; i++) {
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
      '-m1000',
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
