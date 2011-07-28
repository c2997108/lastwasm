#! /bin/sh

# Exercise LAST programs, and compare the output to a reference
# output.  More tests should be added!

cd $(dirname $0)

alias lastdb=../src/lastdb
alias lastal=../src/lastal

dnaSeq=galGal3-M-32.fa
protSeq=Q2LCP8.fa
fastq=SRR001981-1k.fastq
mat=../examples/HOXD70
gc=../examples/vertebrateMito.gc
seed=../examples/yass.seed
db=/tmp/last-test

{
    echo TEST 1  # spaced seeds, soft-masking, centroid alignment, matrix file
    lastdb -c -m110 $db $dnaSeq
    lastal -u1 -j5 -p $mat -x3400 -e2500 $db $dnaSeq
    echo

    echo TEST 2  # multiple volumes & query batches
    lastdb -m1 -s1 $db $dnaSeq
    lastal -f0 -i1 -w0 $db $dnaSeq
    echo

    echo TEST 3  # match-counting, with multiple query batches
    lastal -j0 -i1 -s0 $db $dnaSeq
    echo

    echo TEST 4  # FASTQ quality scores
    lastal -Q1 -e90 -a9 $db $fastq
    echo

    echo TEST 5  # translated alignment & genetic code file
    lastdb -p $db $protSeq
    lastal -F12 -e40 -G $gc $db $dnaSeq
    echo

    echo TEST 6  # subset seed file, soft-masking
    lastdb -c -u $seed $db $dnaSeq
    lastal -s0 -f0 -e18 $db $dnaSeq
    echo

    echo TEST 7  # asymmetric scoring matrix
    lastal -s0 -f0 -p asymmetric.mat -e2000 $db $dnaSeq
    echo

    echo TEST 8  # FASTQ-Illumina quality scores
    lastdb -m1111110 $db $dnaSeq
    lastal -Q3 -e110 $db illumina100.txt
    echo

    echo TEST 9  # PRB-format quality data
    lastal -Q4 -e90 $db mouse_tss_prb.txt
    echo

    echo TEST 10  # probabilistic alignment with quality scores
    lastal -Q1 -j6 -e90 -a9 $db $fastq
    echo
} |
grep -v version |  # omit header lines with the LAST version number
diff last-test.out -

# Test: lastdb, lastal, last-reduce-alignments.sh (maf-sort.sh,
# last-remove-dominated.py, maf-swap.py), maf-join.py
cd ../examples
./multiMito.sh | diff multiMito.maf -

rm $db*
