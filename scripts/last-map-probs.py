#! /usr/bin/env python

# Copyright 2010 Martin C. Frith

# Read query-genome alignments: write them along with the probability
# that each alignment is not the true mapping of its query.  These
# probabilities make the risky assumption that one of the alignments
# reported for each query is correct.

# The input must be readable twice, so not a pipe.  (This is because
# the algorithm needs 2 passes over the input, and we want to cope
# with input that is too big to all fit in memory.  An alternative
# would be to write it to a temporary file.)

import sys, os, fileinput, math, optparse, signal

def logsum(x, y):
    """log(exp(x) + exp(y))."""
    a = max(x, y)
    b = min(x, y)
    return a + math.log(1 + math.exp(b-a))

def addScore(score, temperature, queryName, denominators):
    if temperature < 0:
        raise Exception("I need a header line with: t=(a positive value)")
    x = score / temperature
    y = denominators.get(queryName, -1e9)
    denominators[queryName] = logsum(x, y)

def mismapProb(score, temperature, queryName, denominators):
    x = score / temperature
    y = denominators[queryName]
    assert x <= y
    prob = 1 - math.exp(x - y)
    assert prob >= 0
    return prob

def mafScore(words):
    for word in words:
        if word.startswith("score="):
            return float(word[6:])
    raise Exception("found an alignment without a score")

def lastMapProbs(opts, args):
    temperature = -1
    denominators = {}

    for line in fileinput.input(args):
        words = line.split()
        if not words: pass
        elif line.startswith('#'):
            for i in words:
                if i.startswith('t='): temperature = float(i[2:])
        elif words[0].isdigit():  # we have an alignment in tabular format
            score = float(words[0])
            queryName = words[6]
            addScore(score, temperature, queryName, denominators)
        elif words[0] == "a":  # we have maf format
            score = mafScore(words)
            sLineCount = 0
        elif words[0] == "s":  # we have maf format
            sLineCount += 1
            if sLineCount == 2:  # the query sequence is on the 2nd s line
                queryName = words[1]
                addScore(score, temperature, queryName, denominators)

    storedLines = None
    for line in fileinput.input(args):
        words = line.split()
        if not words: pass
        elif words[0].isdigit():  # we have an alignment in tabular format
            score = float(words[0])
            queryName = words[6]
            p = mismapProb(score, temperature, queryName, denominators)
            if score < opts.score or p > opts.mismap: continue
            newLineEnd = "\t%g\n" % p
            line = line.rstrip() + newLineEnd
        elif words[0] == "a":  # we have maf format
            score = mafScore(words)
            sLineCount = 0
            storedLines = []  # start storing lines instead of printing them
        elif words[0] == "s":  # we have maf format
            sLineCount += 1
            if sLineCount == 2:  # the query sequence is on the 2nd s line
                queryName = words[1] 
                p = mismapProb(score, temperature, queryName, denominators)
                if score < opts.score or p > opts.mismap: continue
                newLineEnd = " mismap=%g\n" % p
                storedLines[0] = storedLines[0].rstrip() + newLineEnd
                for i in storedLines: print i,
                storedLines = None
        if storedLines is None: print line,
        else: storedLines.append(line)
        if not words: storedLines = None  # reset at end of maf paragraph

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = """
  %prog --help
  %prog [options] lastal-alignments"""

    description = "Calculate a mismap probability for each alignment.  This is the probability that the alignment does not reflect the origin of the query sequence, assuming that one reported alignment does reflect the origin of each query."

    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-m", "--mismap", type="float", default=0.01, metavar="M",
                  help="don't write alignments with mismap probability > M (default: %default)")
    op.add_option("-s", "--score", type="float", default=0, metavar="S",
                  help="don't write alignments with score < S (default: %default)")
    (opts, args) = op.parse_args()
    if not args: op.error("please give me a filename")
    if '-' in args: op.error("sorry, can't use '-' (standard input)")

    try: lastMapProbs(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
