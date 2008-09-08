#! /usr/bin/env python

# Change the order of sequences in MAF alignments, then make sure the
# top sequence is on the + strand.  The MAF spacing might get messed up.

import fileinput, string, re

complement = string.maketrans('ACGTNSWRYKMBDHVacgtnswrykmbdhv',
                              'TGCANSWYRMKVHDBtgcanswyrmkvhdb')

def revcomp(seq):
    return seq[::-1].translate(complement)

def flipstrand(strand):
    return '+' if strand == '-' else '-'

def flip(line):
    words = re.split(r'(\s+)', line)  # keep the spaces
    start = int(words[4])
    alnsize = int(words[6])
    seqsize = int(words[10])
    newstart = seqsize - start - alnsize
    words[4] = str(newstart)
    words[8] = flipstrand(words[8])
    words[12] = revcomp(words[12])
    return ''.join(words)

def write(lines):
    if lines:
        lines.reverse()
        topstrand = lines[0].split()[4]
        for line in lines:
            if topstrand == '-':
                print flip(line),
            else:
                print line,

lines = []

for line in fileinput.input():
    if line.startswith('s'):
        lines.append(line)
    else:
        write(lines)
        lines = []
        print line,

write(lines)
