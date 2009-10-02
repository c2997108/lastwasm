#! /usr/bin/env python
# Convert MAF format alignments to tabular format
import fileinput, optparse

op = optparse.OptionParser(usage="%prog my-alignments.maf")
(opts, args) = op.parse_args()

def gap_string(gap1, gap2):
    return str(gap1) + ':' + str(gap2)

def aln2blocks(seq1, seq2):
    if len(seq1) != len(seq2): raise Exception('bad MAF data')
    blocks = []
    size = gap1 = gap2 = 0
    for i, j in zip(seq1, seq2):
        if i == '-' or j == '-':
            if size != 0:
                blocks.append(str(size))
                size = 0
            if i != '-': gap1 += 1
            if j != '-': gap2 += 1
        else:
            if gap1 != 0 or gap2 != 0:
                blocks.append(gap_string(gap1, gap2))
                gap1 = gap2 = 0
            size += 1
    if size != 0:
        blocks.append(str(size))
    elif gap1 != 0 or gap2 != 0:
        blocks.append(gap_string(gap1, gap2))
    return blocks

def doit(aln):
    if len(aln) == 0: return
    if len(aln) != 3: raise Exception('bad pairwise MAF data')
    score = aln[0][1].split('=')[1]
    seq1 = aln[1][6]
    seq2 = aln[2][6]
    blocks = aln2blocks(seq1, seq2)
    blkstr = ','.join(blocks)
    print '\t'.join([score] + aln[1][1:6] + aln[2][1:6] + [blkstr])

aln = []

for line in fileinput.input(args):
    words = line.split()
    if line.startswith('a'):
        doit(aln)
        aln = [words]
    elif line.startswith('s'):
        if len(aln) == 0: raise Exception('bad MAF data')
        aln.append(words)

doit(aln)
