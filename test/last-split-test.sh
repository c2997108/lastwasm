#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../src:$PATH

maf=SRR359290-1k.maf

{
    last-split -h

    last-split $maf

    last-split -c0 $maf

    last-split -t0 $maf

    last-split -M7.5 -S2 $maf

    last-split -m0.001 -s180 $maf

    last-split -n $maf
} |
diff last-split-test.out -
