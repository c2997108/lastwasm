WebAssembly builds
==================

This section gives a quick overview of building LAST tools for
execution inside a web browser using WebAssembly.

Requirements
------------

* Emscripten 3 or later (``apt-get install emscripten``)
* Python environment if you want to run ``maf-convert`` with Pyodide.

Basic setup
-----------

1. Copy ``/usr/share/emscripten/.emscripten`` to ``~/.emscripten`` and
   change ``FROZEN_CACHE = True`` to ``FROZEN_CACHE = False``.  Set the
   ``EM_CONFIG`` environment variable to point to this file.

2. Build ``lastal`` and ``lastdb`` with ``emmake``::

     emmake make -C src lastal lastdb \\
       CXXFLAGS="-s USE_ZLIB=1" LDFLAGS="-s USE_ZLIB=1"

   This produces ``src/lastal`` and ``src/lastdb`` in WebAssembly
   module form (``lastal.wasm`` and ``lastdb.wasm``) together with their
   JavaScript wrappers.

   The other C/C++ programs in ``src`` such as ``last-split``,
   ``last-merge-batches``, and ``last-pair-probs`` can be built in the
   same way if you need them inside the browser::

     emmake make -C src last-split last-merge-batches last-pair-probs \
       CXXFLAGS="-s USE_ZLIB=1" LDFLAGS="-s USE_ZLIB=1"

Testing
-------

After building, verify that the modules run correctly under Node.js by
executing ``test/webassembly-test.sh``::

     ./test/webassembly-test.sh

This script runs ``lastdb`` and ``lastal`` on a tiny example and checks
that they produce output.

Running maf-convert in the browser
----------------------------------

``maf-convert`` is a Python program.  You can run it in the browser by
loading it into Pyodide_::

     import micropip
     await micropip.install('maf-convert')
     import maf_convert

This allows conversion of the alignment files that ``lastal`` outputs.

.. _Pyodide: https://pyodide.org/
