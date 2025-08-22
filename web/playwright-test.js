// Headless smoke test for the WASM page using Playwright.
// - Serves ./web via Python http.server
// - Opens the page, selects sample FASTA files
// - Runs lastdb-only, verifies it completes
// - Runs full lastdb+lastal, verifies completion

const { chromium } = require('playwright');
const { spawn } = require('child_process');
const path = require('path');

async function waitForServer(url, timeoutMs = 15000) {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    try {
      const res = await fetch(url);
      if (res.ok) return;
    } catch {}
    await new Promise(r => setTimeout(r, 200));
  }
  throw new Error('server did not start in time');
}

async function main() {
  const webDir = path.resolve(__dirname);
  // Start a simple static server on an ephemeral port (no conflicts)
  const http = require('http');
  const fs = require('fs');
  const mime = {
    '.html': 'text/html', '.js': 'application/javascript', '.mjs': 'application/javascript',
    '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json',
    '.txt': 'text/plain', '.map': 'application/json'
  };
  function serve(req, res) {
    let p = req.url.split('?')[0];
    if (p === '/') p = '/index.html';
    const f = path.join(webDir, p);
    fs.readFile(f, (err, data) => {
      if (err) {
        res.writeHead(404); res.end('Not found'); return;
      }
      const ext = path.extname(f);
      res.setHeader('Content-Type', mime[ext] || 'application/octet-stream');
      res.end(data);
    });
  }
  const srv = http.createServer(serve);
  await new Promise((resolve) => srv.listen(0, '127.0.0.1', resolve));
  const port = srv.address().port;

  const browser = await chromium.launch({ headless: true });
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  page.on('console', msg => console.log('[browser]', msg.type(), msg.text()));

  const url = `http://127.0.0.1:${port}/`;
  await page.goto(url, { waitUntil: 'domcontentloaded' });

  // Fill inputs
  const ref = path.join(webDir, 'samples', 'sample_ref.fasta');
  const qry = path.join(webDir, 'samples', 'sample_qry.fasta');
  await page.setInputFiles('#refFasta', ref);
  await page.setInputFiles('#qryFasta', qry);

  // Enable safe mode
  await page.check('#safeMode');

  // Run lastdb only
  await page.click('#runDbBtn');
  try {
    await page.waitForFunction(() => document.querySelector('#status')?.textContent.includes('lastdb完了'), null, { timeout: 20000 });
  } catch (e) {
    const status = await page.textContent('#status');
    const log = await page.textContent('#log');
    console.error('lastdb-only timeout. status=', status, '\nlog=\n', log);
    throw e;
  }
  const logText = await page.textContent('#log');
  if (!/生成ファイル:\s*(?!\(なし\))/.test(logText || '')) {
    throw new Error('lastdb did not generate files');
  }

  // Now run full pipeline (TAB only)
  await page.click('#runBtn');
  try {
    await page.waitForFunction(() => document.querySelector('#status')?.textContent.includes('完了'), null, { timeout: 30000 });
  } catch (e) {
    const status = await page.textContent('#status');
    const log = await page.textContent('#log');
    const tab = await page.textContent('#tab');
    console.error('full run timeout. status=', status, '\nTAB=\n', tab, '\nlog=\n', log);
    throw e;
  }
  const tab = await page.textContent('#tab');
  if (!tab || !tab.includes('$ lastal')) {
    throw new Error('TAB output not found');
  }

  await browser.close();
  srv.close();
  console.log('E2E OK');
}

main().catch(err => {
  console.error('E2E FAILED:', err && err.stack || err);
  process.exit(1);
});
