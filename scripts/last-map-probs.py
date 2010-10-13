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

def logsum(x, y):
    a = max(x, y)
    b = min(x, y)
    return a + math.log(1 + math.exp(b-a))

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
        x = float(words[0]) / temperature
        y = denominators.get(tagName, -1e9)
        denominators[tagName] = logsum(x, y)

for line in fileinput.input(args):
    words = line.split()
    if line.startswith('#') or not words:
        print line,
    else:
        tagName = words[6]
        x = float(words[0]) / temperature
        y = denominators[tagName]
        prob = math.exp(x - y)
        print '\t'.join(words + ['%g' % prob])
