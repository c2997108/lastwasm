#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../src:$PATH

maf=SRR359290-1k.maf

{
    last-split -h

    last-split -m0.01 $maf

    last-split -m0.01 -c0 $maf

    last-split -m0.01 -t0 $maf

    last-split -m0.01 -M7.5 -S2 $maf

    last-split -m0.001 -s180 $maf

    last-split -m0.01 -n $maf

    last-split -d0 -m0.001 -s180 $maf
    last-split -d1 -m0.001 -s180 $maf
    last-split -d2 -m0.001 -s180 $maf

    grep -v '^q' $maf | last-split -m0.001 -s180

    last-split 102.maf
    last-split -fMAF 102.maf

    last-split -d1 spliceWithGap.maf
} | diff -u last-split-test.out -
