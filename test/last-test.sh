#! /bin/sh

# Exercise LAST programs, and compare the output to a reference
# output.  More tests should be added!

cd $(dirname $0)
PATH=$PATH:../src

seq=galGal3-M-32.fa
mat=../examples/HOXD70
db=/tmp/last-test

{
    echo TEST 1  # spaced seeds, soft-masking, centroid alignment, matrix file
    lastdb -c -m110 $db $seq
    lastal -u2 -j5 -p $mat -x3400 -e2500 $db $seq
    echo

    echo TEST 2  # multiple volumes & query batches
    lastdb -s1 $db $seq
    lastal -f0 -i1 -w0 $db $seq
    echo

    echo TEST 3  # match-counting, with multiple query batches
    lastal -j0 -i1 -s0 $db $seq
    echo
} |
grep -v version |  # omit header lines with the LAST version number
diff last-test.out -

rm $db*
