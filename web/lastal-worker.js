import createLastalModule from './lastal.js?v=20260724-wasm-pthreads';

self.onmessage = async event => {
  const {
    index,
    queryText,
    queryLetters,
    dbFiles,
    dbName,
    alArgs,
    useAsmJs,
    memory,
  } = event.data || {};

  const stdout = [];
  const stderr = [];
  const startedAt = performance.now();
  try {
    self.postMessage({ type: 'progress', index, queryLetters });
    const module = await createLastalModule({
      noInitialRun: true,
      useAsmJs,
      INITIAL_MEMORY: memory,
      locateFile: path => new URL(path, import.meta.url).href,
      print: line => stdout.push(line),
      printErr: line => stderr.push(line),
    });

    const fs = module.FS;
    if (!fs.analyzePath('/work').exists) fs.mkdir('/work');
    fs.chdir('/work');

    for (const file of dbFiles || []) {
      fs.writeFile(file.name, file.data);
    }
    fs.writeFile('query.fa', queryText || '');

    runMain(module, ['-f', 'TAB', '-P', '1', ...(alArgs || []), dbName, 'query.fa'], 'lastal', stderr);
    self.postMessage({
      ok: true,
      index,
      tabText: stdout.length > 0 ? `${stdout.join('\n')}\n` : '',
      alLog: stderr,
      elapsedMs: performance.now() - startedAt,
    });
  } catch (error) {
    const detail = error?.stack || error?.message || String(error);
    self.postMessage({
      ok: false,
      index,
      error: detail,
      alLog: stderr,
    });
  }
};

function runMain(module, argv, program, stderrLines) {
  try {
    module.callMain(argv);
  } catch (error) {
    const status = typeof error?.status === 'number' ? error.status : undefined;
    if (error?.name === 'ExitStatus' && status === 0) return;
    const rawError = error?.stack || error?.message || String(error);
    const detail = stderrLines.length > 0 ? `\n${stderrLines.join('\n')}` : `\n${rawError}`;
    const exitInfo = typeof status === 'number' ? ` exited with code ${status}` : ' failed';
    throw new Error(`${program}${exitInfo}.${detail}`);
  }
}
