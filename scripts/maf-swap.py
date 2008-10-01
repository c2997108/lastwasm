#! /usr/bin/env python

# Change the order of sequences in MAF alignments, then make sure the
# top sequence is on the + strand.  The MAF spacing might get messed up.

import fileinput, string, re

complement = string.maketrans('ACGTNSWRYKMBDHVacgtnswrykmbdhv',
                              'TGCANSWYRMKVHDBtgcanswyrmkvhdb')

def revcomp(seq):
    return seq[::-1].translate(complement)

def flipstrand(strand):
    if strand == '-': return '+'
    else:             return '-'

def flip_s_line(line):
    words = re.split(r'(\s+)', line)  # keep the spaces
    start = int(words[4])
    alnsize = int(words[6])
    seqsize = int(words[10])
    newstart = seqsize - start - alnsize
    words[4] = str(newstart)
    words[8] = flipstrand(words[8])
    words[12] = revcomp(words[12])
    return ''.join(words)

def flip_p_line(line):
    words = re.split(r'(\s+)', line)  # keep the spaces
    words[2:-2] = words[-3:1:-1]
    return ''.join(words)

def flip_line(line):
    if   line.startswith('s'): return flip_s_line(line)
    elif line.startswith('p'): return flip_p_line(line)
    else:                      return line

def write(lines):
    if not lines: return
    lines[1], lines[2] = lines[2], lines[1]
    if lines[1].split()[4] == '-':
        lines = map(flip_line, lines)
    print ''.join(lines)  # this prints a blank line at the end

lines = []

for line in fileinput.input():
    if line.startswith('#'):
        print line,
    elif line.isspace():
        write(lines)
        lines = []
    else:
        lines.append(line)

write(lines)
