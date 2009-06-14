#! /usr/bin/env python

# Change the order of sequences in MAF alignments, then make sure the
# top sequence is on the + strand.  The MAF spacing might get messed up.

import optparse, fileinput, string, re

parser = optparse.OptionParser()
parser.add_option("-n", type="int", default=2,
                  help="move the Nth sequence to the top (default: %default)")
(opts, args) = parser.parse_args()
if opts.n < 1: parser.error("option -n: should be >= 1")

complement = string.maketrans('ACGTNSWRYKMBDHVacgtnswrykmbdhv',
                              'TGCANSWYRMKVHDBtgcanswyrmkvhdb')

def revcomp(seq):
    return seq[::-1].translate(complement)

def flipstrand(strand):
    if strand == '-': return '+'
    else:             return '-'

def flip_s_line(words):
    start = int(words[4])
    alnsize = int(words[6])
    seqsize = int(words[10])
    newstart = seqsize - start - alnsize
    words[4] = str(newstart)
    words[8] = flipstrand(words[8])
    words[12] = revcomp(words[12])

def flip_p_line(words):
    words[2:-2] = words[-3:1:-1]

def flip_q_line(words):
    words[2] = words[2][::-1]

def flip_line(line):
    words = re.split(r'(\s+)', line)  # keep the spaces
    if   line.startswith('s'): flip_s_line(words)
    elif line.startswith('p'): flip_p_line(words)
    elif line.startswith('q'): flip_q_line(words)
    return ''.join(words)

def indexOfNthSequence(lines, n):
    for i, line in enumerate(lines):
        if line.startswith('s'):
            if n == 1: return i
            n -= 1
    parser.error("option -n: should be <= the number of sequences")

def write(lines):
    if not lines: return
    indexToMove = indexOfNthSequence(lines, opts.n)
    lines.insert(1, lines.pop(indexToMove))
    if lines[1].split()[4] == '-':
        lines = map(flip_line, lines)
    print ''.join(lines)  # this prints a blank line at the end

lines = []

for line in fileinput.input(args):
    if line.startswith('#'):
        print line,
    elif line.isspace():
        write(lines)
        lines = []
    else:
        lines.append(line)

write(lines)
