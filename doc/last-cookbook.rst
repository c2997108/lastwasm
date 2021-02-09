LAST Cookbook
=============

LAST_ is used by running commands in a terminal / command line.  It
has many options: unfortunately, the LAST developers don't know the
best options for every possible alignment task.  Here are some
reasonable starting points.  Feel free to optimize (and share!) search
protocols.

A minimal example: compare human & fugu mitochondrial genomes
-------------------------------------------------------------

Let's find and align similar regions between the human and fugu
mitochondrial genomes.  These FASTA-format files are in LAST's
examples directory: ``humanMito.fa`` and ``fuguMito.fa``.  The
simplest possible usage is::

  lastdb humdb humanMito.fa
  lastal humdb fuguMito.fa > myalns.maf

The lastdb_ command creates several files whose names begin with
"humdb".  The lastal_ command then compares fuguMito.fa to the humdb
files, and writes the alignments to a file called "myalns.maf".

Understanding the output
------------------------

The output has very long lines, so you need to view it without
line-wrapping.  For example, you can use::

  less -S myalns.maf

Each alignment looks like this (MAF_ format)::

  a score=27 EG2=4.7e+04 E=2.6e-05
  s humanMito 2170 145 + 16571 AGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTT...
  s fuguMito  1648 142 + 16447 AGTAGGCTTAGAAGCAGCCACCA--CAAGAAAGCGTT...

The score is a measure of how significant the similarity is.  EG2 and
E are explained at last-evalues_.  Lines starting with "s" contain:
the sequence name, the start coordinate of the alignment, the number
of bases spanned by the alignment, the strand, the sequence length,
and the aligned bases.

The start coordinates are zero-based.  This means that, if the
alignment begins right at the start of a sequence, the coordinate is
0.  If the strand is "-", the start coordinate is the coordinate in
the reverse-complemented sequence (the same as if you were to
reverse-complement the sequence before giving it to LAST).

You can convert MAF to other formats with maf-convert_, or use lastal_
option ``-f`` to get a few other formats.

More accurate: learn substitution & gap rates
---------------------------------------------

We can get more accurate alignments between the human and fugu
mitochondrial genomes like this::

  lastdb humdb humanMito.fa
  last-train humdb fuguMito.fa > hufu.train
  lastal -p hufu.train humdb fuguMito.fa > myalns.maf

The last-train_ command finds the rates of deletion, insertion, and
each kind of substitution between these sequences, and writes them to
a file "hufu.train".  lastal's ``-p`` option uses this file to get
more-accurate alignments.

Comparing protein sequences
---------------------------

We can compare some query proteins to some reference proteins like
this::

  lastdb -p -c -R01 refdb ref-prots.fa
  lastal refdb query-prots.fa > prot-alns.maf

``-p`` tells it the sequences are proteins.  (If you forget ``-p`` and
the sequences look proteinaceous, you'll get a warning message.)

The other options suppress alignments caused by simple sequence such
as ``APPSPAPPSPAPPSPAPPSPAP``:

* ``-R01`` converts the sequence letters to uppercase, then finds
  simple regions and converts them to lowercase.  This will be done
  for both ref-prots and query-prots.

* ``-c`` omits alignments that lack a significant amount of
  uppercase-to-uppercase alignment.

You can also use last-train_, but we've hardly tested it for
protein-protein alignment, so we're not sure if it helps.

Find high-similarity, and short, protein alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we just want high-similarity alignments, we can use the PAM30 (or
PAM10) `scoring scheme`_::

  lastdb -p -c -R01 refdb ref-prots.fa
  lastal -p PAM30 refdb query-prots.fa > prot-alns.maf

This has two advantages:

* It omits low-similarity alignments, or alignment parts.

* It can find short similarities, which would be deemed insignificant
  (likely to occur by chance between random sequences) unless we focus
  the search on high-similarity.

Comparing DNA to proteins
-------------------------

We can find related regions of DNA and proteins, allowing for nonsense
mutations and frameshifts.  For example, let's find DNA regions
related to transposon proteins::

  lastdb -q -c -R01 trandb transposon-prots.fa
  last-train --codon trandb dna.fa > codon.train
  lastal -p codon.train -m100 -D1e9 -K1 trandb dna.fa > out.maf

``-q`` appends ``*`` meaning STOP to each protein, and treats ``*`` as
a 21st protein letter.

``--codon`` makes it do DNA-versus-protein.  Here, last-train_ tries
to learn 21x64 substitution rates, so it needs a fairly large amount
of data (e.g. a chromosome).

``-m100`` makes it more slow-and-sensitive than the default (which is
``-m10``), see lastal_.

``-D1e9`` sets a strict significance_ threshold.  It means: only
report strong similarities that would be expected to occur by chance,
between random sequences, at most once per 10^9 base-pairs.  The
default is 1e6.

``-K1`` streamlines the output by omitting any alignment whose DNA
range lies in that of a higher-scoring alignment.

Another possibility is to add last-train_ option ``-X1``, which treats
matches to ``X`` (unknown) amino acids as neutral instead of
disfavored.

You can reuse ``last-train`` output for different alignment tasks, if
you expect the rates to be roughly the same.

Aligning high-indel-error long DNA reads to a genome
----------------------------------------------------

Suppose we have DNA reads in either FASTA or FASTQ_ format.  This is
sensitive but slow::

  lastdb -P8 -uNEAR -R01 mydb genome.fa
  last-train -P8 -Q0 mydb reads.fastq > reads.train
  lastal -P8 -p reads.train mydb reads.fastq | last-split > out.maf

``-P8`` uses 8 parallel threads, adjust as appropriate for your
computer.  This has no effect on the results.

``-uNEAR`` selects a `seeding scheme`_ that's better at finding
alignments with few substitutions and/or many gaps.

``-Q0`` makes it discard the fastq_ quality information (or you can
keep-but-ignore it with ``-Qkeep``).

last-split_ finds a unique best alignment for each part of each read.
It gives each alignment a `mismap probability`_, which is high if that
part of the read is almost equally similar to several parts of the
genome.

Here we used ``-R01`` to lowercase simple sequence like
``cacacacacacacacacacacaca``.  But we didn't suppress it with ``-c``,
so as not to hide anything from last-split_.  If desired, you can
filter lowercase with last-postmask_.

You can go faster by sacrificing a bit of sensitivity.  It depends on
your aim, e.g. slow-and-sensitive seems necessary to find intricate
rearrangements of repeats.  Suggested ways to go faster:

* `Mask repeats`_.  This has often worked well.

* Add lastal_ option ``-k8`` (or ``-k16`` etc).  This makes it faster,
  by only finding initial matches starting at every 8th (or 16th etc)
  position in the reads.

* Replace ``-uNEAR`` with ``-uRY32`` (or ``-uRY16``, ``-uRY8``,
  ``-uRY4``).  This makes it check for initial matches starting at
  only ~1/32 (or ~1/16 etc) of positions, in both reads and genome.
  Compared to ``-k``: this harms sensitivity slightly more, but
  reduces memory use and makes lastdb faster.

Which genome version to use?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some genome versions (e.g. for human) have artificial
exactly-duplicated regions, which makes it hard to align reads
uniquely.  To avoid that, look for a genome version called something
like "analysis set".

You can use multiple genomes, which will be treated like one big
genome::

  lastdb -P8 -uNEAR -R01 mydb human.fa virus.fa other-genomes.fa

Aligning low-error long DNA reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can do this the same way as for high-error reads, but perhaps
accelerate more aggressively (e.g. ``-uRY32``).

If repeats are not masked, lastal_ option ``-C2`` may reduce run time
with little effect on accuracy.

Aligning potentially-spliced RNA or cDNA long reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See here_.  (For low-error reads, you can probably omit ``-d90`` and
``-m20``.)

Aligning Illumina DNA reads to a genome
---------------------------------------

::

  lastdb -P8 -uNEAR -R01 -C2 mydb genome.fasta
  last-train -P8 -Q1 mydb reads.fastq.gz > reads.train
  lastal -P8 -p reads.train mydb reads.fastq.gz | last-split | gzip > out.maf.gz

Most LAST commands accept ``.gz`` compressed files, and you can
compress output with ``gzip`` as above.

lastdb_ option ``-C2`` makes the alignment a bit faster, but uses more
memory.  This has no effect on the results.  (You could use it in the
other examples too, but it might not be faster.)

``-Q1`` makes it use the fastq_ quality information to improve the
training and alignment.  LAST **assumes** that the qualities reflect
substitution errors, not insertion/deletion errors.  (For long
non-Illumina reads, we suspect this assumption doesn't hold, so we
didn't use this option.)

This recipe may be excessively slow-and-sensitive.  Adding lastal_
option ``-C2`` may make it faster with negligible accuracy loss.  You
can accelerate with e.g. ``-uRY16`` or ``-k16`` as above.

Finding very short DNA alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, LAST only reports significant_ alignments that will rarely
occur by chance.  In the preceding example, the minimum alignment
length is about 26 bases for a human-size genome (less for smaller
genomes).  To find shorter alignments, add lastal_ option ``-D100``
(say), to get alignments that could occur by chance once per hundred
query letters (the default is once per million.)  This makes the
minimum alignment length about 20 bases for a human-size genome.

Aligning paired-end Illumina DNA reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You could use last-pair-probs_.  It has a disadvantage: it doesn't
allow different parts of one read (i.e. one "end") to align to
different parts of the genome.  Alternatively, you could align the
reads individually, ignoring the pair relationships::

  fastq-interleave reads1.fq reads2.fq | lastal -P8 -p reads.train mydb | last-split > out.maf

``fastq-interleave`` ensures that each read has a unique name (by
appending "/1" and "/2" if necessary).

Aligning potentially-spliced Illumina reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See last-split_ (and last-pair-probs_).

Aligning human & chimp genomes
------------------------------

The aim of genome-genome alignment is discussed in `our paper`_.  Here
is a slow-and-sensitive recipe::

  lastdb -P8 -uNEAR -R01 humdb human_no_alt_analysis_set.fa
  last-train -P8 --revsym -E0.05 -C2 humdb chimp.fa > humchi.train
  lastal -E0.05 -C2 -p humchi.train humdb chimp.fa | last-split -fMAF+ > humchi1.maf

``--revsym`` makes the substitution rates the same on both strands.
For example, it makes A→G equal T→C (because A→G on one strand means
T→C on the other strand).  This is usually appropriate for
genome-genome comparison (but maybe not for mitochondria which have
asymmetric "heavy" and "light" strands).

``-E0.05`` means only get significant_ alignments that would be
expected to occur by chance at a rate ≤ 0.05 times per pair of random
sequences of length 1 billion each.

``-fMAF+`` makes it show `per-base mismap probabilities`_: the
probability that each query (chimp) base should be aligned to a
different part of the reference (human).

The result so far is asymmetric: each part of the chimp genome is
aligned to at most one part of the human genome, but not vice-versa.
We can get one-to-one alignments like this::

  maf-swap humchi1.maf | last-split > humchi2.maf

Here, last-split_ gets parts of the humchi1 alignments.  It uses the
humchi1 *per-base* mismap probabilities to get the humchi2
*per-alignment* mismap probabilities.

Then we can discard less-confident alignments, and convert_ to a
compact tabular format::

  last-postmask humchi2.maf | maf-convert -n tab | awk -F= '$2 <= 1e-5' > humchi.tab

last-postmask_ discards alignments caused by simple sequence.  The
``awk`` command gets alignments with `mismap probability`_ ≤ 10^-5.
Finally, we can make a dotplot_::

  last-dotplot humchi.tab humchi.png

**To go faster** with minor accuracy loss: replace ``-uNEAR`` with
``-uRY32`` and/or `mask repeats`_.

To squeeze out the last 0.000...1% of accuracy: add ``-m50`` to the
lastal_ options.

Aligning human & mouse genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can do this in the same way as human/chimp, except that ``-uNEAR``
should be omitted.  To increase sensitivity, but also time and memory
use, add lastdb seeding_ option ``-uMAM4`` or or ``-uMAM8``.  To
increase them even more, add lastal_ option ``-m100`` (or as high as
you can bear).

Large reference sequences
-------------------------

If the sequences that you give to lastdb exceed ~4 billion letters,
consider using 8-byte LAST (lastdb8_ and lastal8_).  Ordinary (4-byte)
LAST can't handle so much sequence at once, so lastdb_ splits it into
"volumes", which may be inefficient.  8-byte LAST avoids voluming, but
uses more memory.  So lastdb8_ works well with a memory-reducing
option: ``-uRY`` or ``-w`` or ``-W``.

Moar faster
-----------

* `Using multiple CPUs / cores <doc/last-parallel.rst>`_
* `Various speed & memory options <doc/last-tuning.rst>`_

Ambiguity of alignment columns
------------------------------

Consider this alignment::

  TGAAGTTAAAGGTATATGAATTCCAATTCTTAACCCCCCTATTAAACGAATATCTTG
  |||||||| ||||||  |  ||  | |  |    || ||||||   |||||||||||
  TGAAGTTAGAGGTAT--GGTTTTGAGTAGT----CCTCCTATTTTTCGAATATCTTG

The middle section has such weak similarity that its precise alignment
cannot be confidently inferred.  We can see the confidence of each
alignment column with lastal_ option ``-j4``::

  lastal -j4 -p hufu.train humdb fuguMito.fa > myalns.maf

The output looks like this::

  a score=17 EG2=9.3e+09 E=5e-06
  s seqX 0 57 + 57 TGAAGTTAAAGGTATATGAATTCCAATTCTTAACCCCCCTATTAAACGAATATCTTG
  s seqY 0 51 + 51 TGAAGTTAGAGGTAT--GGTTTTGAGTAGT----CCTCCTATTTTTCGAATATCTTG
  p                %*.14442011.(%##"%$$$$###""!!!""""&'(*,340.,,.~~~~~~~~~~~

The "p" line indicates the probability that each column is wrongly
aligned, using a compact code (based on ASCII_, the same as fastq_
format):

======  =================   ======  =================
Symbol  Error probability   Symbol  Error probability
======  =================   ======  =================
``!``   0.79 -- 1           ``0``   0.025 -- 0.032
``"``   0.63 -- 0.79        ``1``   0.02  -- 0.025
``#``   0.5  -- 0.63        ``2``   0.016 -- 0.02
``$``   0.4  -- 0.5         ``3``   0.013 -- 0.016
``%``   0.32 -- 0.4         ``4``   0.01  -- 0.013
``&``   0.25 -- 0.32        ``5``   0.0079 -- 0.01
``'``   0.2  -- 0.25        ``6``   0.0063 -- 0.0079
``(``   0.16 -- 0.2         ``7``   0.005  -- 0.0063
``)``   0.13 -- 0.16        ``8``   0.004  -- 0.005
``*``   0.1  -- 0.13        ``9``   0.0032 -- 0.004
``+``   0.079 -- 0.1        ``:``   0.0025 -- 0.0032
``,``   0.063 -- 0.079      ``;``   0.002  -- 0.0025
``-``   0.05  -- 0.063      ``<``   0.0016 -- 0.002
``.``   0.04  -- 0.05       ``=``   0.0013 -- 0.0016
``/``   0.032 -- 0.04       ``>``   0.001  -- 0.0013
======  =================   ======  =================

Note that each alignment is grown from a "core" region, and the
ambiguity estimates assume that the core is correctly aligned.  The
core is indicated by "~" symbols, and it contains exact matches only.

.. _last: README.rst
.. _lastdb8:
.. _lastdb: doc/lastdb.rst
.. _lastal8:
.. _lastal: doc/lastal.rst
.. _dotplot: doc/last-dotplot.rst
.. _last-pair-probs: doc/last-pair-probs.rst
.. _last-postmask: doc/last-postmask.rst
.. _per-base mismap probabilities:
.. _mismap probability:
.. _last-split: doc/last-split.rst
.. _last-train: doc/last-train.rst
.. _convert:
.. _maf-convert: doc/maf-convert.rst
.. _scoring scheme: doc/last-matrices.rst
.. _seeding scheme:
.. _seeding: doc/last-seeds.rst
.. _last-evalues:
.. _significant:
.. _significance: doc/last-evalues.rst
.. _fastq: https://doi.org/10.1093/nar/gkp1137
.. _here:
.. _mask repeats: https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md
.. _MAF: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
.. _ASCII: https://en.wikipedia.org/wiki/ASCII
.. _our paper: https://doi.org/10.1186/s13059-015-0670-9
