Dotplot with LAST-JS on Browser (https://c2997108.github.io/lastwasm/)

LAST: find & align related regions of sequences
===============================================

LAST is designed for moderately large data (e.g. genomes, DNA reads,
proteomes).  It's especially good at:

* Finding rearrangements and recombinations: we believe last-split_
  does that more rigorously than anything else.

* Finding DNA-versus-protein related regions, especially protein_
  fossils_.

* Unusual data, e.g. AT-rich DNA, because we can fit_ parameters to
  the data and calculate significance_.

* Sensitive DNA-DNA search, due to fitting_, sensitive seeding_, and
  calculating significance_.

It can also: indicate the confidence/uncertainty of each column in an
alignment, and use sequence quality data in a rigorous fashion.

Usage
-----

Please see the cookbook_.  **Warning:** this documentation may not
apply to older versions of LAST!  You can see your version with::

  lastal --version

Install
-------

You can install it from bioconda_ or `Debian Med`_, or like this...

Download the highest version number from
https://gitlab.com/mcfrith/last/-/tags.  Using the command line, go
into the downloaded directory and type::

  make

This assumes you have a C++ compiler.  On Linux, you might need to
install a package called "g++".  On Mac, you might need to install
command-line developer tools.  On Windows, you might need to install
Cygwin.  You might also need to install something like "zlib-devel".

For ARM CPUs, the default "make" seems to work in some cases but not
others (sigh).  This seems to be good for ARM::

  make CXXFLAGS="-mcpu=native -O3 -pthread"

It's possible to specify a compiler like this: ``make CXX=MyOtherCompiler``.
If you re-run ``make`` in different ways, it may be good to do ``make clean``
first, to remove any previously-made files.

The programs are in the ``bin`` directory.  For convenient usage, set
up your computer to find them automatically.  Some possible ways:

* Copy the programs to a standard directory: ``sudo make install``
  (using "sudo" to request administrator permissions).

* Copy the programs to your personal bin directory: ``make install prefix=~``

* Adjust your `PATH variable`_.

You might have to log out and back in before your computer recognizes
the new programs.

Further info
------------

Details & citation: `LAST papers`_

LAST is distributed under the GNU General Public License, either
version 3 of the License, or (at your option) any later version.

LAST is brought to you by:

* `Computational Omics Research Team`_, AIRC_
* GSFS_, `University of Tokyo`_
* `AIST-Waseda University CBBD-OIL`_

.. _fit:
.. _fitting: doc/last-train.rst
.. _last-split: doc/last-split.rst
.. _seeding: doc/last-seeds.rst
.. _significance: doc/last-evalues.rst
.. _cookbook: doc/last-cookbook.rst
.. _LAST papers: doc/last-papers.rst
.. _protein: https://doi.org/10.1109/TCBB.2022.3177855
.. _fossils: https://doi.org/10.1093/molbev/msac068
.. _bioconda: https://bioconda.github.io/
.. _Debian Med: https://www.debian.org/devel/debian-med/
.. _PATH variable: https://en.wikipedia.org/wiki/PATH_(variable)
.. _Computational Omics Research Team: https://www.airc.aist.go.jp/en/cort/
.. _AIRC: https://www.airc.aist.go.jp/en/
.. _GSFS: https://www.k.u-tokyo.ac.jp/index.html.en
.. _University of Tokyo: https://www.u-tokyo.ac.jp/en/
.. _AIST-Waseda University CBBD-OIL: https://unit.aist.go.jp/cbbd-oil/en/

JavaScript/Web
--------------

The browser UI in ``index.html`` runs the same LAST algorithm as the WASM build
without special-casing input files.  ``web/lastdb.asm.js`` is the database
builder's Binaryen ``wasm2js`` module.  ``web/lastal.js`` is the pthread WASM
loader, while ``web/lastal-asm.js`` contains the single-thread asm.js fallback
generated from the same LAST source.  The browser builds the database once, then runs one
``lastal`` WASM instance with a shared linear memory and LAST's pthreads.  Query
records are loaded in 100 KiB batches so search progress updates throughout
large runs.  ``lastal -P`` divides each batch across the pthread pool, though
smaller batches can reduce parallel efficiency.  Without cross-origin isolation,
the browser falls back to one asm.js search worker generated from the same LAST
build, avoiding duplicated 256 MiB search runtimes.

Open the page through a local HTTP server.  Direct ``file://`` browsing is not
supported because browsers block ES module loading from local files.  The
included ``coi-serviceworker.js`` enables COOP/COEP on static hosts so browsers
can use the faster shared-memory WASM runtime.  The app waits for service-worker
activation before starting and automatically reloads when isolation must be
applied.

- Run locally:

  * ``python -m http.server 8000``
  * Open ``http://localhost:8000/``

- Browser implementation files:

  * ``web/main.js``: UI and dot-plot pipeline coordination
  * ``web/regl-dotplot.js``: full-count regl/WebGL2 renderer, axes, pan/zoom,
    and GPU picking
  * ``web/regl.min.js``: locally bundled regl 2.1.1 runtime
  * ``web/jslast.js``: high-level runner that calls compiled ``lastdb`` and
    dispatches compiled ``lastal`` search jobs
  * ``web/jslast-runner-worker.js``: background runner so indexing/search does
    not block browser UI updates
  * ``web/plot-parser-worker.js``: parallel TAB parser and full-count coordinate
    builder using shared input and transferable typed arrays
  * ``web/lastal-worker.js``: asm.js fallback Web Worker entry point
  * ``web/lastdb.js``: Emscripten browser glue for the database builder
  * ``web/lastal.js``: Emscripten pthread WASM browser glue
  * ``web/lastal-asm.js``: Emscripten asm.js fallback runtime
  * ``web/lastdb.asm.js``: Binaryen ``wasm2js`` database builder
  * ``web/build_asm.sh``: regenerates the single-thread asm.js fallback

- Regression test:

  * ``npm run test:web`` generates a fresh TAB baseline with the Node WASM CLI
    and compares it byte-for-byte with an 8-thread browser search.
  * ``npm run test:pthread`` checks that one shared WASM instance executes
    LAST's pthread search concurrently and reports effective CPU usage.

The requested thread count is capped by the number of query FASTA records.  A
single query record therefore uses one search thread.  The UI shows elapsed time
and pthread search progress based on completed query chunks.  It also shows an
approximate remaining time using a small browser-local timing model from previous
runs when available.  Other LAST options are passed through to the compiled LAST
programs rather than reinterpreted by a separate JavaScript search engine.

Large TAB results are transferred as a UTF-8 buffer.  The page immediately
shows the first 20 KiB before full-result processing and provides the complete
output through ``Download TAB``.
The dot plot uses up to eight parser workers and renders matching alignments as
WebGL2 line geometry; it does not sample the result.  The minimum plotted
hidden alignment-length cutoff is configurable in the input panel and defaults
to 1000 bp; alignments of 1000 bp or less remain in TAB output but are omitted
from the plot.  CPU metadata and
GPU geometry are kept in compact typed arrays instead of per-alignment JavaScript
objects.  Up to 100,000 deterministic alignment lines are emphasized to match
the former Plotly view, while all remaining alignments are still drawn at lower
opacity.  Axes and dotted sequence boundaries are drawn on a lightweight canvas
overlay, and pan, wheel zoom, fit, and hover details are handled without
rebuilding the TAB output.
Reference and query axes follow the original FASTA sequence order by default,
including sequences without plotted alignments.  This can be disabled in the
input panel to use TAB alignment-output order instead.
Detailed ``dotplot-timing`` and ``dotplot-longtasks`` log lines separate TAB
parsing, worker messaging, regl buffer upload, GPU completion, and final browser
paint time.

WASM/Node (legacy experimental)
-------------------------------

This repository still contains experimental WebAssembly builds for Node.js.

- Build for Node:

  * ``npm run build:node`` (builds ``node/lastdb.js`` + ``.wasm``)
  * ``npm run build:node:al`` (builds ``node/lastal.js`` + ``.wasm``)

- CLI wrappers (installed via ``npm i`` or using ``npx``):

  * ``lastdb-wasm -v -C1 /tmp/mydb test/huma.fa``
  * ``lastal-wasm -v /tmp/mydb test/ttttt.fa > out.maf``

Notes:

- ``lastdb`` for Node uses ``-sNODERAWFS=1`` and reads/writes directly to the host filesystem.
- 並列化: Emscripten pthreads を用いた並列化（``-P N``）に対応しています。
  デフォルトでは安定性のため ``-P`` を 1 に丸めますが、環境変数
  ``LASTDB_ALLOW_THREADS=1`` を付与すればそのまま使用できます。
  例: ``LASTDB_ALLOW_THREADS=1 lastdb-wasm -v -P 2 -C1 /tmp/mydb test/huma.fa``
- ``lastal`` wrapper copies the database and query files into an in-memory FS and writes results to stdout.
