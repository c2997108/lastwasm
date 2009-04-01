#! /bin/sh

# This script demonstrates using LAST and maf-join.py to construct a
# multiple alignment of the human, mouse, chicken, and fugu
# mitochondrial genomes.

# You have to run this from inside the examples directory, else it
# won't find the files.

# Put the LAST programs and scripts into the command search path:
PATH=$PATH:../src:../scripts

# Make a LAST database of the human sequence:
lastdb -c -m110 humanMito humanMito.fa

# Align the mouse sequence to the human sequence:
# Let's use a score threshold of 25.  The accompanying E-value tables
# show that this is a reasonable threshold for sequences of this size.
lastal -u2 -e25 -j4 humanMito mouseMito.fa > hm.maf

# Remove paralogs (if any):
# this also sorts the alignments into the right order for maf-join.py
last-reduce-alignments.sh -d hm.maf > hm2.maf

# Align the chicken sequence to the human sequence:
lastal -u2 -e25 -j4 humanMito chickenMito.fa > hc.maf
last-reduce-alignments.sh -d hc.maf > hc2.maf

# Align the fugu sequence to the human sequence:
lastal -u2 -e25 -j4 humanMito fuguMito.fa > hf.maf
last-reduce-alignments.sh -d hf.maf > hf2.maf

# Join the pairwise alignments into a multiple alignment:
maf-join.py hm2.maf hc2.maf hf2.maf

# Clean up the intermediate files that we made:
rm humanMito*.??? h*.maf
