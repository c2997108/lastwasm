#! /bin/sh

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=../src:../scripts:$PATH

db=/tmp/$(basename $0 .sh)

trap 'rm -f $db*' EXIT

{
    lastdb $db ../examples/humanMito.fa
    try last-train -m1 $db ../examples/mouseMito.fa
    try last-train -m1 --revsym $db ../examples/mouseMito.fa
    try last-train -m1 --matsym --gapsym $db ../examples/mouseMito.fa
} |
grep -v '^# lastal' |
diff -u $(basename $0 .sh).out -
