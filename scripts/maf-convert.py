#! /usr/bin/env python

# Copyright 2010 Martin C. Frith

# Read MAF-format alignments: write them in other formats.

import sys, os, fileinput, optparse, itertools, signal

##### General-purpose routines: #####

def flatten(listOfLists): return sum(listOfLists, [])

def scoreFromMaf(mafLines):
    for word in mafLines[0]:
        if word.startswith("score="):
            return word.split('=')[1]
    raise Exception("encountered an alignment without a score")

##### Routines for converting to AXT format: #####

def axtDataFromMafLine(sLine):
    chr = sLine[1]
    mafBeg = int(sLine[2])
    mafLen = int(sLine[3])
    axtBeg = str(mafBeg + 1)  # convert from 0-based to 1-based coordinate
    axtEnd = str(mafBeg + mafLen)
    strand = sLine[4]
    return [chr, axtBeg, axtEnd, strand]

axtCounter = itertools.count()

def writeAxt(maf):
    count = str(axtCounter.next())
    score = scoreFromMaf(maf)
    sLines = [i for i in maf if i[0] == "s"]
    if not sLines: raise Exception("empty alignment")
    sequences = [i[6] for i in sLines]
    axtData = map(axtDataFromMafLine, sLines)
    dataHead, dataTail  = axtData[0], axtData[1:]
    if dataHead[3] != "+":
        raise Exception("for AXT, the 1st strand in each alignment must be +")
    print " ".join([count] + dataHead[:3] + flatten(dataTail) + [score])
    for i in sequences: print i
    print  # print a blank line at the end

##### Routines for converting to tabular format: #####

def gapString(gap1, gap2):
    return str(gap1) + ":" + str(gap2)

def symbolSize(symbol, letterSize):
    if symbol == "-": return 0
    elif symbol == "\\": return 1
    elif symbol == "/": return -1
    else: return letterSize

def matchAndGapSizes(sequence1, sequence2, letterSize1, letterSize2):
    if len(sequence1) != len(sequence2):
        raise Exception('aligned sequences have different lengths')
    matchSize = 0
    gap1 = gap2 = 0
    for i, j in zip(sequence1, sequence2):
        if '-' in (i, j):
            if matchSize != 0:
                yield str(matchSize)
                matchSize = 0
            gap1 += symbolSize(i, letterSize1)
            gap2 += symbolSize(j, letterSize2)
        else:
            if gap1 != 0 or gap2 != 0:
                yield gapString(gap1, gap2)
                gap1 = gap2 = 0
            matchSize += 1  # this is correct for translated alignments too
    if matchSize != 0:
        yield str(matchSize)
    elif gap1 != 0 or gap2 != 0:
        yield gapString(gap1, gap2)

def gapSizePerLetter(sLine):
    sequence = sLine[6]
    if "/" in sequence or "\\" in sequence: return 3
    gaps = sequence.count("-")
    nonGaps = len(sequence) - gaps
    mafLen = int(sLine[3])
    if mafLen == nonGaps * 3: return 3
    return 1

def writeTab(maf):
    score = scoreFromMaf(maf)
    sLines = [i for i in maf if i[0] == "s"]
    if len(sLines) != 2: raise Exception("pairwise alignments only, please")
    outWords = [i[1:6] for i in sLines]
    sequences = [i[6] for i in sLines]
    sizes = map(gapSizePerLetter, sLines)
    gapInfo = matchAndGapSizes(sequences[0], sequences[1], sizes[0], sizes[1])
    gapString = ','.join(gapInfo)
    print "\t".join([score] + flatten(outWords) + [gapString])

##### Routines for reading MAF format: #####

def mafInput(lines):
    maf = []
    for line in lines:
        if line.startswith("#"):
            print line,
        elif line.isspace():
            if maf: yield maf
            maf = []
        else:
            maf.append(line.split())
    if maf: yield maf

def formatName(formatString):
    s = formatString.lower()
    if   "axt".startswith(s): return "axt"
    elif "tabular".startswith(s): return "tab"
    else: raise Exception("unknown format: " + formatString)

def mafConvert(args):
    format = formatName(args[0])
    for maf in mafInput(fileinput.input(args[1])):
        if format == "axt": writeAxt(maf)
        else: writeTab(maf)

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = """
  %prog axt my-alignments.maf
  %prog tab my-alignments.maf"""

    op = optparse.OptionParser(usage=usage)
    (opts, args) = op.parse_args()
    if len(args) != 2: op.error("I need a format-name and a file-name")

    try: mafConvert(args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
