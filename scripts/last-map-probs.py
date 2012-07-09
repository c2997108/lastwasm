#! /usr/bin/env python

# Copyright 2010, 2011, 2012 Martin C. Frith

# Read query-genome alignments: write them along with the probability
# that each alignment is not the true mapping of its query.  These
# probabilities make the risky assumption that one of the alignments
# reported for each query is correct.

import sys, os, fileinput, math, optparse, signal

def logsum(x, y):
    """log(exp(x) + exp(y))."""
    a = max(x, y)
    b = min(x, y)
    return a + math.log(1 + math.exp(b-a))

def mafScore(aLine):
    for word in aLine.split():
        if word.startswith("score="):
            return float(word[6:])
    raise Exception("found an alignment without a score")

def namesAndScores(lines):
    queryNames = []
    scores = []
    for line in lines:
        if line[0] == "a":  # faster than line.startswith("a")
            s = mafScore(line)
            scores.append(s)
            sLineCount = 0
        elif line[0] == "s":
            sLineCount += 1
            if sLineCount == 2: queryNames.append(line.split(None, 2)[1])
        elif line[0].isdigit():  # we have an alignment in tabular format
            w = line.split(None, 7)
            scores.append(float(w[0]))
            queryNames.append(w[6])
    return queryNames, scores

def scoreTotals(queryNames, scores, temperature):
    denominators = {}
    for n, s in zip(queryNames, scores):
        d = denominators.get(n, -1e9)
        denominators[n] = logsum(d, s / temperature)
    return denominators

def writeOneBatch(lines, queryNames, scores, denominators, opts, temperature):
    isWanted = True
    i = 0
    for line in lines:
        if line[0] == "a":
            s = scores[i]
            p = 1.0 - math.exp(s / temperature - denominators[queryNames[i]])
            i += 1
            if s < opts.score or p > opts.mismap:
                isWanted = False
            else:
                newLineEnd = " mismap=%g\n" % p
                line = line.rstrip() + newLineEnd
        elif line[0].isdigit():  # we have an alignment in tabular format
            s = scores[i]
            p = 1.0 - math.exp(s / temperature - denominators[queryNames[i]])
            i += 1
            if s < opts.score or p > opts.mismap: continue
            newLineEnd = "\t%g\n" % p
            line = line.rstrip() + newLineEnd
        if isWanted: print line,
        if line.isspace(): isWanted = True  # reset at end of maf paragraph

def doOneBatch(lines, opts, temperature):
    queryNames, scores = namesAndScores(lines)
    denominators = scoreTotals(queryNames, scores, temperature)
    writeOneBatch(lines, queryNames, scores, denominators, opts, temperature)

def readHeaderOrDie(lines):
    for line in lines:
        if line.startswith("#") or line.isspace():
            print line,
            for i in line.split():
                if i.startswith("t="):
                    return float(i[2:])
        else:
            raise Exception("I need a header line with: t=(a positive value)")

def lastMapProbs(opts, args):
    f = fileinput.input(args)
    temperature = readHeaderOrDie(f)
    lines = []

    for line in f:
        if line.startswith("# batch"):
            doOneBatch(lines, opts, temperature)
            lines = []
        lines.append(line)
    doOneBatch(lines, opts, temperature)

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
    if not args and sys.stdin.isatty():
        op.print_help()
        op.exit()

    try: lastMapProbs(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
