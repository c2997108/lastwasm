#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    maf-cut chr1:152413-158286 split1.maf

    maf-cut chrM:150-200 bs100.maf bs100.fastq

    awk '(NR-1) % 4 < 2' bs100.fastq | tr '@' '>' |
	maf-cut chrM:150-200 bs100.maf -

    maf-cut chrUn_KI270748v1:2579-2609 frameshift-new.maf

    maf-cut UN-L1MA6_pol#LINE/L1:884-892 frameshift-new.maf
} | diff -u maf-cut-test.out -
