#! /usr/bin/env python

# Copyright 2012 Martin C. Frith

# This script reads alignments of query sequences to a genome, and
# estimates which alignment parts represent the source of each query.
# It allows different parts of one query to come from different parts
# of the genome.

# The input should be in MAF format
# (http://genome.ucsc.edu/FAQ/FAQformat.html#format5).  It must
# include a header of the sort written by lastal, with score
# parameters and genome size.

# This script should work with Python 2.x, x>=4.

import fileinput, math, operator, optparse, os, signal, string, sys
from itertools import *

import numpy

minQualCode = 33
numQualCodes = 64
maxQualCode = minQualCode + numQualCodes - 1
maxQualSymbol = chr(maxQualCode)

# The "rescale"s are there to avoid numeric overflow.  Mathematically
# speaking, they could be removed (i.e. set to 1) with no effect.

def forward(baseProbs, midProbs, alignBegs, jumpProb):
    dpProbs = [0] * len(alignBegs)
    begProb = 1
    rescale = 1
    yield dpProbs, rescale

    for i, x, y in izip(count(), baseProbs, midProbs):
        probFromJump = sum(dpProbs) * jumpProb
        dpProbs = [(a * c + probFromJump + (0, begProb)[d == i]) * b / rescale
                   for a, b, c, d in izip(dpProbs, x, y, alignBegs)]
        begProb /= rescale
        rescale = sum(dpProbs) + 1
        yield dpProbs, rescale

def backward(baseProbs, midProbs, alignEnds, jumpProb, rescales):
    dpProbs = [0] * len(alignEnds)
    endProb = 1
    yield dpProbs

    r = xrange(len(midProbs))
    for i, x, y, z in izip(reversed(r), reversed(baseProbs), reversed(midProbs), reversed(rescales)):
        probFromJump = sum(dpProbs) * jumpProb
        dpProbs = [(a * c + probFromJump + (0, endProb)[d == i]) * b / z
                   for a, b, c, d in izip(dpProbs, x, y, alignEnds)]
        endProb /= z
        yield dpProbs

def prod(n):
    return reduce(operator.mul, n)

def partitionFunction(dpMatrix, alignEnds, rescales):
    return sum(dpMatrix[x][i] / prod(rescales[x:])
               for i, x in enumerate(alignEnds))

def marginalProbs(iBeg, j, qAln, f, b, z, baseProbs, midProbs, rescales):
    i = iBeg
    for x in qAln:
        if x == "-":
            yield f[i][j] * b[i][j] * midProbs[i][j] / rescales[i] / z
        else:
            yield f[i+1][j] * b[i][j] / baseProbs[i][j] / z
            i += 1

def viterbi(baseScores, midScores, alignBegs, jumpScore):
    dpScores = [-1e9] * len(alignBegs)
    yield dpScores

    for i, x, y in izip(count(), baseScores, midScores):
        scoreFromJump = max(dpScores) + jumpScore
        dpScores = [max(a + c, scoreFromJump, (-1e9, 0)[d == i]) + b
                    for a, b, c, d in izip(dpScores, x, y, alignBegs)]
        yield dpScores

def endpoint(dpMatrix, alignEnds):
    endScores = [dpMatrix[x][i] for i, x in enumerate(alignEnds)]
    sMax = max(endScores)
    jMax = endScores.index(sMax)  # breaks ties by minimizing j
    iMax = alignEnds[jMax]
    return sMax, iMax, jMax

def traceback(dpMatrix, midScores, alignBegs, jumpScore, iMax, jMax):
    j = jMax
    i = iMax
    iEnd = i

    while j > -1:
        i -= 1

        score = dpMatrix[i][j] + midScores[i][j]
        jNew = -2

        if alignBegs[j] == i and score <= 0:
            score = 0
            jNew = -1

        s = max(dpMatrix[i])
        if s + jumpScore > score:
            jNew = dpMatrix[i].index(s)

        if jNew > -2:
            yield j, i, iEnd
            iEnd = i
            j = jNew

def baseScoresFromAln(rAln, qAln, qQual, scoreTable, firstGapScore):
    """Yield the alignment score for each base of the query sequence."""
    for i, j, q in izip(rAln, qAln, qQual):
        if   j == "-": pass  # there is no query base (deletion in query)
        elif i == "-": yield firstGapScore  # treat all gaps as 1st gaps
        else:          yield scoreTable[i][j][q]  # substitution score

def midScoresFromAln(rAln, qAln, gapExistenceScore, gapExtensionScore):
    """Yield the alignment score between each pair of query bases."""
    delScore = 0
    iPrevious = None
    for i, j in izip(rAln, qAln):
        if j == "-":  # deletion in query
            assert i != "-"  # xxx (might happen, in bad MAF input)
            if delScore == 0: delScore = gapExistenceScore
            delScore += gapExtensionScore
        else:
            if i == "-" and iPrevious == "-":
                yield -gapExistenceScore  # correction for treating gaps as 1st
            else:
                yield delScore
            delScore = 0
        iPrevious = i
    yield delScore

def mafParts(maf):
    """Yield each "s" line, with its corresponding "q" line, if any."""
    sLine = None
    for line in maf:
        if line[0] == "s":
            if sLine: yield sLine, qLine
            sLine = line
            qLine = None
        elif line[0] == "q":
            qLine = line
    if sLine: yield sLine, qLine

def convertMafLine(w):  # modifies "w"
    if w[0] == "s":
        w[2] = int(w[2])
        w[3] = int(w[3])
        w[5] = int(w[5])
        w[3] += w[2]  # replace size with end

def unconvertMafLine(w):  # modifies "w"
    if w[0] == "s":
        w[3] -= w[2]  # replace end with size
        w[2] = str(w[2])
        w[3] = str(w[3])
        w[5] = str(w[5])

def flipMafStrands(maf):  # modifies "maf"
    for i in maf:
        if i[0] == "s":
            seqLen = i[5]
            i[2], i[3] = seqLen - i[3], seqLen - i[2]
            i[6] = i[6][::-1]
        elif i[0] == "q":
            i[2] = i[2][::-1]
        elif i[0] == "p":
            i[1] = i[1][::-1]

def addEdgeGaps(maf, qSeqLine, begGaps, endGaps):  # modifies "maf"
    for i in maf:
        if i is qSeqLine:
            i[2] -= begGaps
            i[3] += endGaps
            i[6] = "." * begGaps + i[6] + "." * endGaps
        elif i[0] == "s":
            i[6] = "-" * begGaps + i[6] + "-" * endGaps
        elif i[0] == "q":
            i[2] = "!" * begGaps + i[2] + "!" * endGaps
        elif i[0] == "p":
            i[1] = "!" * begGaps + i[1] + "!" * endGaps

def getMafData(maf):  # modifies "maf"
    r, q = mafParts(maf)
    rSeqLine, rQualLine = r
    qSeqLine, qQualLine = q

    if rQualLine:
        raise Exception("I can't handle quality data for the genomic sequence")

    qBeg, qEnd, qStrand, qSeqLen = qSeqLine[2:6]

    qFin = qSeqLen - qEnd
    addEdgeGaps(maf, qSeqLine, qBeg, qFin)

    rAln = rSeqLine[6]
    qAln = qSeqLine[6]

    if qQualLine: qQual = qQualLine[2]
    else:         qQual = maxQualSymbol * len(qAln)

    return qStrand, qBeg, qEnd, rAln, qAln, qQual

# *** Routines for getting slices of alignments:

def seqPosFromAlnPos(alnPos, aln):
    return alnPos - aln[:alnPos].count("-")

def mafSlice(maf, alnBeg, alnEnd):
    yield ["a"]
    for i in maf:
        if i[0] == "s":
            beg = i[2] + seqPosFromAlnPos(alnBeg, i[6])
            end = i[2] + seqPosFromAlnPos(alnEnd, i[6])
            aln = i[6][alnBeg:alnEnd]
            yield [i[0], i[1], beg, end, i[4], i[5], aln]
        if i[0] == "q":
            yield [i[0], i[1], i[2][alnBeg:alnEnd]]
        if i[0] == "p":
            yield [i[0], i[1][alnBeg:alnEnd]]

def nthBaseIndex(sequenceWithGaps, n):
    for i, x in enumerate(sequenceWithGaps):
        # simple loop is faster than messing with generators
        if x != "-":
            if n: n -= 1
            else: return i
    assert 0

def mafSliceBeg(rAln, qAln, qBeg):
    alnBeg = nthBaseIndex(qAln, qBeg)
    numInserts = nthBaseIndex(rAln[alnBeg:], 0)  # trim initial insertions
    return alnBeg + numInserts, qBeg + numInserts

def mafSliceEnd(rAln, qAln, qEnd, qLen):
    qFin = qLen - qEnd
    alnFin, qFin = mafSliceBeg(rAln[::-1], qAln[::-1], qFin)
    qEnd = qLen - qFin
    alnEnd = len(qAln) - alnFin
    return alnEnd, qEnd

# *** Routines for pretty-printing MAF format:

def maxlen(s):
    return max(imap(len, s))

def sLineFieldWidths(mafLines):
    sLines = (i for i in mafLines if i[0] == "s")
    sColumns = izip(*sLines)
    return map(maxlen, sColumns)

def joinedMafS(words, fieldWidths):
    formatParams = chain(*izip(fieldWidths, words))
    return "%*s %-*s %*s %*s %*s %*s %*s" % tuple(formatParams)

def joinedMafLine(words, fieldWidths):
    if words[0] == "s":
        return joinedMafS(words, fieldWidths)
    elif words[0] == "q":
        words = words[:2] + [""] * 4 + words[2:]
        return joinedMafS(words, fieldWidths)
    elif words[0] == "p":
        words = words[:1] + [""] * 5 + words[1:]
        return joinedMafS(words, fieldWidths)
    else:
        return " ".join(words)

def printMaf(maf):  # modifies "maf"
    for i in maf:
        unconvertMafLine(i)
    fieldWidths = sLineFieldWidths(maf)
    for i in maf:
        print joinedMafLine(i, fieldWidths)
    print

# *** End of routines for pretty-printing MAF format

def asciiFromProb(probRight):
    """Probability -> phred score in fastq-sanger ASCII representation."""
    probWrong = 1 - probRight
    e = max(probWrong, 1e-10)  # avoid overflow errors
    s = int(-10 * math.log10(e))  # phred score, with fractions rounded down
    return chr(min(s + 33, 126))

def doOneQueryWithOneAlignment(maf, opts):  # modifies "maf"
    aLine = maf[0]
    for i in aLine:
        if i.startswith("score="): score = int(i[6:])
    if "score" not in vars():
        raise Exception("found an alignment without a score")
    if score < opts.score: return
    aLine.append("mismap=1e-10")
    qSeqLine = queryLine(maf)
    q = "~" * len(qSeqLine[6])  # not quite correct: neglects self-jumps
    maf.append(["p", q])
    if qSeqLine[4] == "-": flipMafStrands(maf)
    printMaf(maf)

def doOneQuery(mafs, scoreTable, gop, gep, jumpCost, scale, opts):
    if len(mafs) == 1:  # this special-case is unnecessary, but faster
        return doOneQueryWithOneAlignment(mafs[0], opts)

    d = imap(getMafData, mafs)
    qStrands, qBegs, qEnds, rAlns, qAlns, qQuals = zip(*d)

    tBaseScores = map(list, imap(baseScoresFromAln, rAlns, qAlns, qQuals,
                                 repeat(scoreTable), repeat(-gop-gep)))
    baseScores = zip(*tBaseScores)

    tMidScores = map(list, imap(midScoresFromAln, rAlns, qAlns,
                                repeat(-gop), repeat(-gep)))
    midScores = zip(*tMidScores)

    qLen = len(baseScores)

    baseProbs = [[math.exp(j / scale) for j in i] for i in baseScores]
    midProbs = [[math.exp(j / scale) for j in i] for i in midScores]
    jumpProb = math.exp(-jumpCost / scale)

    f, rescales = zip(*forward(baseProbs, midProbs, qBegs, jumpProb))
    b = list(backward(baseProbs, midProbs, qEnds, jumpProb, rescales))
    b.reverse()
    z = partitionFunction(f, qEnds, rescales)

    if opts.no_split:
        tb = [(j, 0, qLen) for j in xrange(len(mafs))]
    else:
        dp = list(viterbi(baseScores, midScores, qBegs, -jumpCost))
        sMax, iMax, jMax = endpoint(dp, qEnds)
        tb = list(traceback(dp, midScores, qBegs, -jumpCost, iMax, jMax))
        tb.reverse()

    for j, iBeg, iEnd in tb:
        alnBeg, qBeg = mafSliceBeg(rAlns[j], qAlns[j], iBeg)
        alnEnd, qEnd = mafSliceEnd(rAlns[j], qAlns[j], iEnd, qLen)
        score = sum(tBaseScores[j][qBeg:qEnd])
        score += sum(tMidScores[j][qBeg+1:qEnd])
        if score < opts.score: continue
        qAln = qAlns[j][alnBeg:alnEnd]
        p = list(marginalProbs(qBeg, j, qAln, f, b, z, baseProbs, midProbs,
                               rescales))
        mismap = 1 - max(p)
        mismap = max(mismap, 1e-10)  # try to avoid floating-point inaccuracy
        if mismap > opts.mismap: continue
        maf = list(mafSlice(mafs[j], alnBeg, alnEnd))
        maf[0].append("score=%d" % score)
        maf[0].append("mismap=%.3g" % mismap)
        q = "".join(imap(asciiFromProb, p))
        #q = " ".join("%.2g" % i for i in p)
        maf.append(["p", q])
        if qStrands[j] == "-": flipMafStrands(maf)
        printMaf(maf)

def queryLine(maf):
    return [i for i in maf if i[0] == "s"][1]  # 2nd line starting with "s"

def queryName(maf):
    return queryLine(maf)[1]

def doOneBatch(mafs, scoreTable, gop, gep, jumpCost, scale, opts):
    mafs.sort(key=queryLine)
    for k, v in groupby(mafs, queryName):
        doOneQuery(list(v), scoreTable, gop, gep, jumpCost, scale, opts)

# *** Routines for parsing and precalculating substitution scores:

def isMatrixHead(w):
    return len(w) > 1 and w[0] == "#" and not [i for i in w if len(i) > 1]

def isMatrixBody(w):
    return len(w) > 2 and w[0] == "#" and len(w[1]) == 1

def probFromPhred(s):
    return 10.0 ** (-0.1 * s)

def generalizedScore(score, scale, phredScore, letterProb):
    r = math.exp(score / scale)
    p = probFromPhred(phredScore)
    u = p / (1 - letterProb)
    c = 1 - u
    return int(round(scale * math.log(c * r + u)))

def matrixLookup(matrix, rowNames, colNames, x, y):
    try:
        row = rowNames.index(x)
        col = colNames.index(y)
    except ValueError:
        return min(chain(*matrix))
    return matrix[row][col]

def tableFromMatrix(rowNames, colNames, matrix, scale):
    """FASTA-versus-FASTA score matrix -> FASTA-versus-FASTQ score table."""
    rowNames = map(str.upper, rowNames)
    colNames = map(str.upper, colNames)
    table = {}

    # Reverse-engineer the abundances of ACGT from the score matrix:
    a = [[math.exp(matrixLookup(matrix, rowNames, colNames, i, j) / scale)
          for j in "ACGT"] for i in "ACGT"]
    b = [1.0] * len(a)
    queryLetterProbs = numpy.linalg.solve(a, b).tolist()
    #print queryLetterProbs

    for i in string.ascii_letters:
        table[i] = {}
        x = i.upper()
        for j in string.ascii_letters:
            table[i][j] = {}
            y = j.upper()
            score = matrixLookup(matrix, rowNames, colNames, x, y)
            for q in range(numQualCodes):
                symbol = chr(q + minQualCode)
                if x in "ACGT" and y in "ACGT":
                    p = queryLetterProbs["ACGT".index(y)]
                    table[i][j][symbol] = generalizedScore(score, scale, q, p)
                else:
                    table[i][j][symbol] = score

    return table

# End of routines for parsing and precalculating substitution scores

def addMaf(mafs, maf):  # modifies "mafs" and "maf"
    if maf:
        queryStrand = queryLine(maf)[4]
        if queryStrand == "-": flipMafStrands(maf)
        mafs.append(maf)

def lastSplitProbs(opts, args):
    state = 0
    maf = []
    mafs = []
    for line in fileinput.input(args):
        w = line.split()
        if state == -1:  # we are reading the score matrix within the header
            if isMatrixBody(w):
                rowNames.append(w[1])
                matrix.append(map(int, w[2:]))
            else:
                state = 0
        if state == 0:  # we are reading the header
            if isMatrixHead(w):
                colNames = w[1:]
                rowNames = []
                matrix = []
                state = -1
            elif line.startswith("#"):
                for i in w:
                    if i.startswith("a="): gop = int(i[2:])  # gap open pen
                    if i.startswith("b="): gep = int(i[2:])  # gap extn pen
                    if i.startswith("t="): scale = float(i[2:])
                    if i.startswith("letters="): genomeSize = int(i[8:])
            elif not line.isspace():
                if "matrix" not in vars():
                    raise Exception("I need a header with score parameters")
                if "gop" not in vars():
                    raise Exception("I need a header line with a=")
                if "gep" not in vars():
                    raise Exception("I need a header line with b=")
                if "scale" not in vars():
                    raise Exception("I need a header line with t=")
                if "genomeSize" not in vars():
                    raise Exception("I need a header line with letters=")
                scoreTable = tableFromMatrix(rowNames, colNames, matrix, scale)
                jumpProb = opts.break_prob / (2 * genomeSize)  # 2 strands
                jumpCost = -int(round(scale * math.log(jumpProb)))
                print "# break prob=%s jump cost=%s" % (opts.break_prob,
                                                        jumpCost)
                print "#"
                state = 1
        if state == 1:  # we are reading alignments
            if line.startswith("# batch "):
                addMaf(mafs, maf)
                maf = []
                doOneBatch(mafs, scoreTable, gop, gep, jumpCost, scale, opts)
                mafs = []
            elif line.isspace():
                addMaf(mafs, maf)
                maf = []
            elif not line.startswith("#"):
                convertMafLine(w)
                maf.append(w)
        if line.startswith("#"): print line,
    addMaf(mafs, maf)
    if mafs: doOneBatch(mafs, scoreTable, gop, gep, jumpCost, scale, opts)

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = "%prog [options] LAST-alignments.maf"

    description = "Read alignments of query sequences to a genome, and estimate the genomic source of each part of each query, allowing different parts of one query to come from different parts of the genome."

    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-m", "--mismap", type="float", default=0.01, metavar="M",
                  help="don't write alignments with mismap probability > M "
                  "(default: %default)")
    op.add_option("-s", "--score", type="int", default=0, metavar="S", help=
                  "don't write alignments with score < S (default: %default)")
    op.add_option("-b", "--break-prob", type="float", default=1e-5,
                  metavar="P",
                  help="breakpoint probability per base (default: %default)")
    op.add_option("-n", "--no-split", action="store_true",
                  help="write the original alignments, not split alignments")
    (opts, args) = op.parse_args()

    try: lastSplitProbs(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
