LAST: Genome-Scale Sequence Comparison
======================================

LAST finds similar regions between sequences, and aligns them.  LAST
is similar to BLAST, but it copes better with giga-scale biological
sequences.  It can also use sequence quality data, and it can indicate
the ambiguity of each column in an alignment.

Requirements
------------

To handle mammalian genomes, it's best if you have at least 10-20
gigabytes of real memory, but you can get by with 2 gigabytes.  To
install the software, you need a C++ compiler.

Setup
-----

Using the command line, go into the src directory and type 'make'.
(If you checked it out using subversion, then type 'make' in the
top-level directory, not the src directory.)  This should make three
programs: lastdb, lastal, and lastex.  If you like, you can copy these
programs to some location in your search path.

Miscellaneous
-------------

For usage instructions, please see the doc directory, especially
last-tutorial.txt.

LAST is distributed under the GNU General Public License, either
version 3 of the License, or (at your option) any later version.  For
details, see COPYING.txt.

LAST includes public domain code kindly provided by Yi-Kuo Yu and
Stephen Altschul at the NCBI.

If you have any questions, comments, or problems, please email:
last (ATmark) cbrc (dot) jp.
If reporting a problem, please describe exactly how to trigger it.

References
----------

If you would like to understand LAST in detail, please look at these
articles.  If you use LAST in your research, please cite one (or more)
of them.

The main algorithms used by LAST:
  Adaptive seeds tame genomic sequence comparison
  SM Kielbasa, R Wan, K Sato, P Horton, MC Frith: Genome Research (in press).

How LAST uses sequence quality data:
  Incorporating sequence quality data into alignment improves DNA read mapping
  MC Frith, R Wan, P Horton: Nucleic Acids Research 2010 38(7):e100.

Choice of score parameters, ambiguity of alignment columns, and
gamma-centroid alignment:
  Parameters for Accurate Genome Alignment
  MC Frith, M Hamada, P Horton: BMC Bioinformatics 2010 11:80.
