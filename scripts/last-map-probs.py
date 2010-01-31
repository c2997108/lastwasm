#! /usr/bin/env python

# Read tag-genome alignments in LAST tabular format: write them with
# mapping probabilities in an extra column.  These probabilities make
# the risky assumption that one of the alignments reported for each
# tag is correct.

# The input must be readable twice, so not a pipe.  (This is because
# the algorithm needs 2 passes over the input, and we want to cope
# with input that is too big to all fit in memory.  An alternative
# would be to write it to a temporary file.)

import sys, os, fileinput, math, optparse, signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # stop spurious error message

progName = os.path.basename(sys.argv[0])

op = optparse.OptionParser(usage="%prog last-tabular-output")
(opts, args) = op.parse_args()

if not args: op.error("please give me a filename")

if '-' in args: op.error("sorry, can't use '-' (standard input)")

def likelihoodRatio(score, t): return math.exp(float(score) / t)

temperature = -1
denominators = {}

for line in fileinput.input(args):
    words = line.split()
    if line.startswith('#') or not words:
        for i in words:
            if i.startswith('t='): temperature = float(i[2:])
    else:
        if temperature < 0:
            sys.exit(progName + 
                     ": I need a header line with: t=(a positive value)")
        tagName = words[6]
        lr = likelihoodRatio(words[0], temperature)
        denominators.setdefault(tagName, 0.0)
        denominators[tagName] += lr

for line in fileinput.input(args):
    words = line.split()
    if line.startswith('#') or not words:
        print line,
    else:
        tagName = words[6]
        lr = likelihoodRatio(words[0], temperature)
        prob = lr / denominators[tagName]
        print '\t'.join(words + ['%g' % prob])
