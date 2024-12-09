#! /bin/sh

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=../bin:$PATH

db=/tmp/$(basename $0 .sh)

trap 'rm -f $db*' EXIT

{
    lastdb --bits=4 -R10 $db ../examples/humanMito.fa
    try "last-train -m1 $db < ../examples/mouseMito.fa"
    try last-train -m1 -C2 --revsym $db ../examples/mouseMito.fa
    try last-train -m1 -k16 --matsym --gapsym $db ../examples/mouseMito.fa

    lastdb -uBISF -S2 $db ../examples/humanMito.fa
    try last-train -S0 -Q1 $db bs100.fastq
} 2>&1 |
    grep -v '^# lastal' | diff -u last-train-test.out -
