#! /bin/sh

# Exercise LAST programs, and compare the output to a reference
# output.  More tests should be added!

cd $(dirname $0)
PATH=$PATH:../src

seq=galGal3-M-32.fa
mat=../examples/HOXD70
ref=last-test.maf
db=/tmp/last-test

lastdb -c -m110 $db $seq

lastal -u2 -j5 -p $mat -x3400 -e2500 $db $seq | diff $ref -

rm $db*
