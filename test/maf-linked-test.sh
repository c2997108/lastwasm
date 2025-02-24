#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    maf-linked homSap-monDom.maf
    maf-linked -c2 -d100000 homSap-monDom.maf
    maf-linked 90089.maf
} | diff -u maf-linked-test.txt -
