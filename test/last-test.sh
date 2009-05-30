#! /bin/sh

# Exercise LAST programs, and compare the output to a reference
# output.  More tests should be added!

cd $(dirname $0)
PATH=$PATH:../src

seq=galGal3-M-32.fa
mat=../examples/HOXD70
ref=last-test.out
db=/tmp/last-test

{
    echo TEST 1
    lastdb -c -m110 $db $seq
    lastal -u2 -j5 -p $mat -x3400 -e2500 $db $seq
    echo

    echo TEST 2
    lastdb -s1 $db $seq
    lastal -f0 -i1 -w0 $db $seq
    echo
} | diff $ref -

rm $db*
