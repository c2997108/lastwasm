// Smoke-test one browser WASM instance using LAST's pthread-based -P option.

const { chromium } = require('playwright');
const fs = require('fs');
const http = require('http');
const path = require('path');

const THREAD_COUNTS = [1, 8];
const QUERY_COPIES = 1000;

async function main() {
  const rootDir = path.resolve(__dirname, '..');
  const server = createServer(rootDir);
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const baseUrl = `http://127.0.0.1:${server.address().port}/`;

  const browser = await chromium.launch({ headless: true });
  const cdp = await browser.newBrowserCDPSession();

  try {
    const control = await runCpuControl(browser, cdp, baseUrl);
    console.log(`8-worker CPU control: ${control.effectiveCores.toFixed(2)} effective cores ` +
      `(${control.cpuSeconds.toFixed(3)} CPU seconds / ${control.elapsedSeconds.toFixed(3)} wall seconds)`);

    const results = [];
    for (const threads of THREAD_COUNTS) {
      results.push(await runCase(browser, cdp, baseUrl, threads));
    }

    console.table(results.map(result => ({
      threads: result.threads,
      elapsedSeconds: result.elapsedSeconds.toFixed(3),
      rendererCpuSeconds: result.rendererCpuSeconds.toFixed(3),
      effectiveCores: result.effectiveCores.toFixed(2),
      pthreadSpawned: result.pthreadSpawned,
      pthreadMaxRunning: result.pthreadMaxRunning,
      outputLines: result.outputLines,
    })));

    const single = results.find(result => result.threads === 1);
    const threaded = results.find(result => result.threads === 8);
    const speedup = single.elapsedSeconds / threaded.elapsedSeconds;
    console.log(`P8 speedup: ${speedup.toFixed(2)}x`);

    if (control.effectiveCores < 2) {
      throw new Error(`CPU control reached only ${control.effectiveCores.toFixed(2)} cores; host parallelism could not be validated`);
    }
    if (threaded.pthreadMaxRunning < 2) {
      throw new Error(`pthreadMaxRunning=${threaded.pthreadMaxRunning}; concurrent pthread execution was not observed`);
    }
    if (threaded.effectiveCores < 1.25) {
      throw new Error(`effectiveCores=${threaded.effectiveCores.toFixed(2)}; renderer CPU time does not show useful parallelism`);
    }
    if (threaded.outputLines !== single.outputLines) {
      throw new Error(`output line count differs: P1=${single.outputLines}, P8=${threaded.outputLines}`);
    }

    console.log('PTHREAD SMOKE TEST OK');
  } finally {
    await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
}

async function runCpuControl(browser, cdp, baseUrl) {
  const page = await browser.newPage();
  try {
    await page.goto(baseUrl, { waitUntil: 'domcontentloaded' });
    const before = await rendererCpuSnapshot(cdp);
    const startedAt = Date.now();
    await page.evaluate(async () => {
      const source = `self.onmessage = () => {
        const end = performance.now() + 2000;
        let value = 0.5;
        while (performance.now() < end) value = Math.sin(value + 1.23456789);
        self.postMessage(value);
      };`;
      const url = URL.createObjectURL(new Blob([source], { type: 'application/javascript' }));
      const workers = Array.from({ length: 8 }, () => new Worker(url));
      await Promise.all(workers.map(worker => new Promise((resolve, reject) => {
        worker.onmessage = resolve;
        worker.onerror = reject;
        worker.postMessage(null);
      })));
      for (const worker of workers) worker.terminate();
      URL.revokeObjectURL(url);
    });
    const elapsedSeconds = (Date.now() - startedAt) / 1000;
    const after = await rendererCpuSnapshot(cdp);
    const cpuSeconds = cpuDelta(before, after);
    return { elapsedSeconds, cpuSeconds, effectiveCores: cpuSeconds / elapsedSeconds };
  } finally {
    await page.close();
  }
}

async function runCase(browser, cdp, baseUrl, threads) {
  const page = await browser.newPage();
  page.on('console', message => console.log(`[browser P${threads}]`, message.type(), message.text()));
  page.on('pageerror', error => console.error(`[browser P${threads}]`, error.stack || error));

  try {
    await page.goto(baseUrl, { waitUntil: 'domcontentloaded' });
    const setup = await page.evaluate(setupCase, { threads, queryCopies: QUERY_COPIES });
    if (!setup.crossOriginIsolated || !setup.sharedMemory) {
      throw new Error(`P${threads} browser page does not have shared WASM memory`);
    }

    const before = await rendererCpuSnapshot(cdp);
    const startedAt = Date.now();
    const result = await page.evaluate(runPreparedCase, threads);
    const elapsedSeconds = (Date.now() - startedAt) / 1000;
    const after = await rendererCpuSnapshot(cdp);
    const rendererCpuSeconds = cpuDelta(before, after);

    return {
      threads,
      elapsedSeconds,
      rendererCpuSeconds,
      effectiveCores: rendererCpuSeconds / elapsedSeconds,
      ...result,
    };
  } finally {
    await page.close();
  }
}

async function setupCase({ threads, queryCopies }) {
  const version = '20260724-plot-profile';
  const lastdbUrl = new URL(`/web/lastdb.js?v=${version}`, location.href).href;
  const lastalUrl = new URL(`/web/lastal.js?v=${version}`, location.href).href;
  const createLastdbModule = (await import(lastdbUrl)).default;
  const createLastalModule = (await import(lastalUrl)).default;
  const refText = await fetch('/test/galGal3-M-32.fa').then(response => response.text());

  const dbStderr = [];
  const lastdb = await createLastdbModule({
    noInitialRun: true,
    useAsmJs: false,
    INITIAL_MEMORY: 134217728,
    locateFile: file => new URL(`/web/${file}`, location.href).href,
    mainScriptUrlOrBlob: lastdbUrl,
    pthreadPoolSize: 0,
    print: () => {},
    printErr: line => dbStderr.push(String(line)),
  });
  prepareFs(lastdb.FS);
  lastdb.FS.writeFile('ref.fa', refText);
  callMain(lastdb, [
    '-P', '1', '--bits=4', '-R00', '-uNEAR', '-w1', '-W1', '-S1', '-C1',
    'refdb', 'ref.fa',
  ], dbStderr);

  const dbFiles = lastdb.FS.readdir('.')
    .filter(name => name.startsWith('refdb.'))
    .map(name => ({ name, data: lastdb.FS.readFile(name) }));
  const sequence = refText.split(/\r?\n/)
    .filter(line => line && !line.startsWith('>'))
    .join('');
  const queryParts = [];
  for (let index = 0; index < queryCopies; index++) {
    queryParts.push(`>query_${index + 1}\n${sequence}\n`);
  }

  window.__pthreadSmokeOutputLines = 0;
  const stderr = [];
  const lastalMemory = new WebAssembly.Memory({
    initial: 268435456 / 65536,
    maximum: 32768,
    shared: true,
  });
  const lastal = await createLastalModule({
    noInitialRun: true,
    useAsmJs: false,
    INITIAL_MEMORY: 268435456,
    wasmMemory: lastalMemory,
    locateFile: file => new URL(`/web/${file}`, location.href).href,
    mainScriptUrlOrBlob: lastalUrl,
    pthreadPoolSize: Math.max(0, threads - 1),
    print: () => { window.__pthreadSmokeOutputLines += 1; },
    printErr: line => stderr.push(String(line)),
  });
  prepareFs(lastal.FS);
  for (const file of dbFiles) lastal.FS.writeFile(file.name, file.data);
  lastal.FS.writeFile('query.fa', queryParts.join(''));

  window.__pthreadSmoke = { lastal, stderr };
  return {
    crossOriginIsolated,
    sharedMemory: lastalMemory.buffer instanceof SharedArrayBuffer,
    queryLetters: sequence.length * queryCopies,
  };

  function prepareFs(fsApi) {
    if (!fsApi.analyzePath('/work').exists) fsApi.mkdir('/work');
    fsApi.chdir('/work');
  }

  function callMain(module, args, stderrLines) {
    try {
      module.callMain(args);
    } catch (error) {
      if (error?.name === 'ExitStatus' && error?.status === 0) return;
      throw new Error(`${error?.stack || error}\n${stderrLines.join('\n')}`);
    }
  }
}

function runPreparedCase(threads) {
  const { lastal, stderr } = window.__pthreadSmoke;
  try {
    lastal.callMain([
      '-f', 'TAB', '-P', String(threads), '-i', '67108864', '-m', '10',
      'refdb', 'query.fa',
    ]);
  } catch (error) {
    if (!(error?.name === 'ExitStatus' && error?.status === 0)) {
      throw new Error(`${error?.stack || error}\n${stderr.join('\n')}`);
    }
  }
  return {
    pthreadSpawned: Number(lastal.pthreadSpawnCount || 0),
    pthreadMaxRunning: Number(lastal.pthreadMaxRunning || 0),
    outputLines: window.__pthreadSmokeOutputLines,
  };
}

async function rendererCpuSnapshot(cdp) {
  const { processInfo } = await cdp.send('SystemInfo.getProcessInfo');
  return new Map(processInfo
    .filter(process => process.type === 'renderer')
    .map(process => [process.id, process.cpuTime]));
}

function cpuDelta(before, after) {
  let total = 0;
  for (const [processId, cpuTime] of after) {
    total += Math.max(0, cpuTime - (before.get(processId) || 0));
  }
  return total;
}

function createServer(rootDir) {
  const mime = {
    '.html': 'text/html',
    '.js': 'application/javascript',
    '.wasm': 'application/wasm',
    '.fa': 'text/plain',
  };

  return http.createServer((request, response) => {
    let pathname = request.url.split('?')[0];
    if (pathname === '/') pathname = '/index.html';
    const file = path.resolve(rootDir, `.${pathname}`);
    if (!file.startsWith(rootDir + path.sep)) {
      response.writeHead(403);
      response.end('Forbidden');
      return;
    }
    fs.readFile(file, (error, data) => {
      if (error) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
      response.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
      response.setHeader('Content-Type', mime[path.extname(file)] || 'application/octet-stream');
      response.end(data);
    });
  });
}

main().catch(error => {
  console.error('PTHREAD SMOKE TEST FAILED:', error?.stack || error);
  process.exitCode = 1;
});
