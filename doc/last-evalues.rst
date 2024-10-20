LAST E-values
=============

It's useful to know whether a similarity is likely to occur by chance
between random sequences (shuffled letters).  Let's see an example::

  a score=37 EG2=72 E=0.00001
  s chr3 9 23 + 939557 TCTGTGAGTTGAAGTTTTCGCCCTAG
  s seqA 2 21 +   2500 TCTGGGAGTTGAAGGTT--GCCCCAG

The score reflects the likelihood that these sequence regions are
related rather than random.  (It doesn't reflect percent identity: a
long, weak similarity can have the same score as a short, strong
similarity.)

The "E-value" is the expected number of similarities with score ≥ 37
between random sequences.  This depends on how long the random
sequences are.

``EG2=`` is the E-value for two random sequences of length 1 billion
each.  (Expected alignments per square giga.)

The meaning of ``E=`` depends on whether you use lastal_ option ``-H``.
``E=`` is the E-value for random sequences with the same lengths as

* all the reference sequences, and all the query sequences (with ``-H``),
* all the reference sequences, and that one query sequence (without ``-H``).

Experience suggests that ``EG2=`` isn't useful, ``E=`` with ``-H`` is
useful, and ``E=`` without ``-H`` is misleading.

Default threshold
-----------------

By default, LAST reports similarities that are expected by chance at
most once per million query letters (for a given database).

Setting a threshold
-------------------

You can make lastal report similarities that are expected by chance at
most once per (say) thousand query letters, with option ``-D``::

  lastal -D1000 humdb fuguMito.fa > myalns.maf

Alternatively, you can use ``-H``::

  lastal -H100 humdb fuguMito.fa > myalns.maf

This will get similarities with ``E=`` ≤ 100.

Advanced issues
---------------

Letter frequencies
~~~~~~~~~~~~~~~~~~

LAST's E-values are for random sequences *with specific letter
frequencies*.  Those frequencies are determined by the substitution
score matrix.  You can see them by running lastal with option -v
(verbose).

If your sequences have very different letter frequencies (e.g. very
AT-rich DNA), the E-values will be misleading.  The best solution is
to use a suitable score matrix, such as AT77 or ATMAP.

Calculation of EG2
~~~~~~~~~~~~~~~~~~

EG2 is calculated from the score like this::

  EG2 = K * exp(-lambda * score) * (1 billion) * (1 billion)

The parameters lambda and K are printed in the header of lastal's
output.

Calculation of E
~~~~~~~~~~~~~~~~

``E=`` is calculated like this.  With ``-H``::

  E = K * exp(-lambda * score) * strands * area(m, n, score) * M/m * N/n

Without ``-H``::

  E = K * exp(-lambda * score) * strands * area(q, n, score) * N/n

* ``strands`` is either 1 or 2 (if both DNA strands are searched).
* ``m`` is the length of the longest query sequence.
* ``n`` is the length of the longest sequence in the database.
* ``M`` is the total length of all query sequences.
* ``N`` is the total length of all database sequences.
* ``q`` is the length of that one query sequence.
* ``area(x, y, score)`` is `slightly less than
  <https://doi.org/10.1186/1756-0500-5-286>`_ ``x * y``.

Effect of option -D
~~~~~~~~~~~~~~~~~~~

Option -D sets the minimum alignment score, by this formula::

  K * exp(-lambda * score) * strands * D * length(n, score) * N / n  ≤  1,

where ``length(n, score)`` is slightly less than ``n``.

Bit score
~~~~~~~~~

Some people like to use "bit scores"::

  Bit score = (lambda * score - ln[K]) / ln[2]
            = log2( 1e18 / EG2 )

Limitations
~~~~~~~~~~~

* E-values cannot be calculated for scoring schemes with weak mismatch
  or gap costs (e.g. match score = 99, mismatch cost = 1).  In such
  cases, lastal requires you to set a score threshold with option -e.

* There may be a long startup time to calculate the E-value parameters
  (lambda, K, and finite-size correction parameters).  This is
  avoided, for particular scoring schemes, by storing the parameters
  in LAST's source code.  Parameters for other scoring schemes can be
  added on request.  (Or you can do it yourself: run lastal with
  option -v, so it prints the E-value parameters, then copy them into
  LastEvaluerData.hh.)

* The E-values do not account for lastal option -c.  If you use this
  option, the E-values will be under-estimates.

* The E-values are for local alignment.  So if you use lastal option
  -T1 (overlap alignment), the E-values will be over-estimates.

Other resources
~~~~~~~~~~~~~~~

For more flexible E-value calculation, try ALP and FALP:
http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/index/software.html

.. _lastal: doc/lastal.rst
