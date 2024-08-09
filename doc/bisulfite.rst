Aligning bisulfite-converted DNA reads to a genome
==================================================

Bisulfite is used to detect methylated cytosines.  It converts
unmethylated Cs to Ts, but it leaves methylated Cs intact.  If we then
sequence the DNA and align it to a reference genome, we can infer
cytosine methylation.

We developed a recipe_ for this in 2011, used in Bisulfighter_.  Since
then, computer memories have gotten bigger, allowing simpler recipes.

Suppose we have bisulfite-converted DNA reads in a file "reads.fastq",
and the genome in "mygenome.fa".  Assume all the reads are from
converted DNA strands, not reverse strands (so they have C→T
conversions and not G→A conversions).  We can align the reads like this::

   lastdb -P8 -uBISF -S2 mydb mygenome.fa
   last-train -P8 -S0 -Q1 mydb reads.fastq > my.train
   lastal -P8 --split -p my.train mydb reads.fastq > out.maf

* ``-uBISF`` selects a `seeding scheme`_ that tolerates c:t mismatches.

* ``-S2`` uses both strands of the genome.  (So we can use just the
  forward strands of the reads).

* ``-S0`` makes no difference here, but matters for reverse-strand reads.

If the DNA reads are from the reverse of converted strands (so they
have G→A conversions), add ``-s0`` to the ``last-train`` options.

Paired DNA reads
----------------

Suppose we have paired reads in 2 files: ``reads1.fastq`` from
converted strands, and ``reads2.fastq`` from the reverse of converted
strands.  We can align the reads like this::

   lastdb -P8 -uBISF -S2 mydb mygenome.fa
   last-train -P8 -S0 -Q1 mydb reads1.fastq > my.train
   lastal -P8 --split -p my.train -sFR mydb reads1.fastq reads2.fastq > out.maf

* ``-sFR`` makes it use forward strands of reads1 and reverse strands
  of reads2.


Avoiding bias
-------------

There is a potential bias: unconverted Cs are easier to align
unambiguously.  We can avoid this bias by converting all Cs in the
reads to Ts::

   perl -pe 'y/Cca-z/ttA-Z/ if $. % 4 == 2' reads.fastq |

   lastal -P8 --split -p my.train mydb |

   perl -F'(\s+)' -ane '$F[12] =~ y/ta/CG/ if /^s/ and $s++ % 2; print @F' > out.maf

The first line converts C to lowercase t, and all other letters to
uppercase.  The last line un-converts lowercase letters in the DNA reads.

.. _recipe: https://doi.org/10.1093/nar/gks275
.. _Bisulfighter: https://github.com/yutaka-saito/Bisulfighter
.. _cookbook: doc/last-cookbook.rst
.. _seeding scheme: doc/last-seeds.rst
