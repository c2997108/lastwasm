maf-linked
==========

maf-linked reads pair-wise sequence alignments in MAF_ format, and
omits isolated alignments.  It keeps groups of alignments that are
nearby in both sequences::

  maf-linked alignments.maf > out.maf

It can read alignments from a pipe like this::

  ... | maf-linked - > out.maf

It may be useful for genome-to-genome alignments: It removes
alignments between non-homologous insertions of homologous transposons
(Frith2022_).

It considers two alignments to be "linked" if, in both sequences, they
are separated by at most D base-pairs and by at most T other
alignments.  It keeps groups of at least C alignments that are linked
directly or indirectly.

Options
-------

-h, --help          show a help message with default option values
-c C, --count C     minimum number of linked alignments
-d D, --distance D  maximum distance between linked alignments
-t T, --tween T     maximum number of other alignments between linked alignments

.. _MAF: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
.. _Frith2022: https://doi.org/10.1093/molbev/msac068
