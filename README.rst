LAST: Genome-Scale Sequence Comparison
======================================

LAST finds and aligns similar regions between sequences.  It's
designed for moderately large data (e.g. genomes, DNA reads,
proteomes).  It's especially good at:

* Finding sequence rearrangements and recombinations: we believe
  last-split_ does that more rigorously than anything else.

* Finding related parts of DNA and proteins, especially protein
  fossils (but introns are not considered).

* Unusual data, e.g. AT-rich DNA, because we can fit_ parameters to
  the data and calculate significance_.

* Sensitive DNA-DNA search, due to fitting_, sensitive seeding_, and
  calculating significance_.

It can also: indicate the confidence/uncertainty of each column in an
alignment, and use sequence quality data in a rigorous fashion.

Requirements
------------

To handle mammalian genomes, it's best if you have at least 10-20
gigabytes of real memory, but you can get by with 2 gigabytes.

To install the software, you need a C++ compiler.  On Linux, you might
need to install a package called "g++".  On Mac, you might need to
install command-line developer tools.  On Windows, you might need to
install Cygwin.

Setup
-----

Using the command line, go into the top-level LAST directory.  To
compile the programs, type::

  make

You might get some harmless warning messages.  It's possible to
specify a C++ compiler like this::

  make CXX=MyOtherCompiler

Install (optional)
------------------

You can copy the programs and scripts to a standard "bin" directory
(using "sudo" to request administrator permissions)::

  sudo make install

Or copy them to your personal bin directory::

  make install prefix=~

You might have to log out and back in before your computer recognizes
the new programs.

Usage
-----

Please see the other files in the doc directory, especially
`<doc/last-tutorial.rst>`_.

Detailed info & citation
------------------------

Please see: `<doc/last-papers.rst>`_

License
-------

LAST (including the scripts) is distributed under the GNU General
Public License, either version 3 of the License, or (at your option)
any later version.

.. _fit:
.. _fitting: doc/last-train.rst
.. _last-split: doc/last-split.rst
.. _seeding: doc/last-seeds.rst
.. _significance: doc/last-evalues.rst
