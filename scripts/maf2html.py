#! /usr/bin/env python

# Copyright 2009, 2010 Martin C. Frith

# Read MAF-format sequence alignments: write them in human-friendly
# HTML, with colours indicating alignment probabilities (from MAF
# lines starting with 'p'.)

import fileinput, operator, itertools, optparse, signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # stop spurious error message

op = optparse.OptionParser(usage="%prog [options] my-alignments.maf")
op.add_option('-l', '--linesize', type="int", default=100, metavar="CHARS",
              help="write CHARS characters per line (default: %default)")
(opts, args) = op.parse_args()
if opts.linesize <= 0: op.error("option -l: should be >= 1")

print '''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
 "http://www.w3.org/TR/html4/strict.dtd">
<html lang="en"><head>
<meta http-equiv="Content-type" content="text/html; charset=UTF-8">
<title>Reliable Alignments</title>
<style type="text/css">
.a {background-color: #3333FF}
.b {background-color: #9933FF}
.c {background-color: #FF66CC}
.d {background-color: #FF3333}
.e {background-color: #FF9933}
.f {background-color: #FFFF00}
</style>
</head><body>
'''

print '''
<div style="line-height:1">
<pre style="display:inline; margin-right:2em"><span class="a">  </span> prob &gt; 0.999</pre>
<pre style="display:inline; margin-right:2em"><span class="b">  </span> prob &gt; 0.99 </pre>
<pre style="display:inline; margin-right:2em"><span class="c">  </span> prob &gt; 0.95 </pre>
<pre style="display:inline; margin-right:2em"><span class="d">  </span> prob &gt; 0.9  </pre>
<pre style="display:inline; margin-right:2em"><span class="e">  </span> prob &gt; 0.5  </pre>
<pre style="display:inline; margin-right:2em"><span class="f">  </span> prob &le; 0.5  </pre>
</div>
'''

def mafInput(lines):
    '''Read lines and yield MAF blocks.'''
    maf = []
    for line in lines:
        if line.startswith('#'): continue
        if line.isspace():
            if maf: yield maf
            maf = []
        else:
            maf.append(line.split())
    if maf: yield maf

def columnSymbol(column):
    '''Return a star if all the letters are identical, or else a space.'''
    firstSymbol = column[0].lower()
    for i in column[1:]:
        if i.lower() != firstSymbol: return ' '
    return '*'

def probabilityClass(probabilityColumn):
    # here, imap seems to be a little bit faster than map:
    try: p = reduce(operator.mul, itertools.imap(float, probabilityColumn))
    except ValueError: return None
    if   p > 0.999: return 'a'
    elif p > 0.99: return 'b'
    elif p > 0.95: return 'c'
    elif p > 0.9: return 'd'
    elif p > 0.5: return 'e'
    else: return 'f'

def identicalRuns(s):  # I want to use groupby, but it's new in python 2.4
    '''Yield (start, end, item) for each run of identical items in s.'''
    for i, x in enumerate(s):
        if not i or x != old:
            if i: yield start, i, old
            start = i
            old = x
    if s: yield start, i+1, x

def htmlSpans(text, spans):
    for beg, end, key in spans:
        chunk = text[beg:end]
        if key: yield '<span class="%s">%s</span>' % (key, chunk)
        else: yield chunk

def cunningCoordinate(sLine):
    '''Parse a MAF start coordinate, dealing cunningly with reverse strands.'''
    coordinate = int(sLine[1])
    if sLine[3] == '+': return coordinate
    else: return coordinate - int(sLine[4]) - 1

def coordinateString(coordinate): return str(abs(coordinate))

def maxlen(s): return max(map(len, s))

def countOfNonGaps(s): return len(s) - s.count('-')

for maf in mafInput(fileinput.input(args)):
    aLines = [i[1:] for i in maf if i[0] == 'a']  # there should be 1 'a' line
    scores = [i for i in aLines[0] if i.startswith('score=')]
    if scores: print '<h3>Alignment %s:</h3>' % scores[0]
    else:      print '<h3>Alignment:</h3>'

    sLines = [i[1:] for i in maf if i[0] == 's']
    names = [i[0] for i in sLines]
    nameWidth = maxlen(names)
    coordinates = map(cunningCoordinate, sLines)
    sequences = [i[5] for i in sLines]
    alignmentLength = maxlen(sequences)  # ?

    pLines = [i[1:] for i in maf if i[0] == 'p']
    if pLines: classes = map(probabilityClass, zip(*pLines))
    else:      classes = itertools.repeat(None)

    print '<pre>'
    for x in range(0, alignmentLength, opts.linesize):
        y = x + opts.linesize
        seqChunks = [s[x:y] for s in sequences]
        classChunk = itertools.islice(classes, x, y)
        spans = list(identicalRuns(classChunk))
        columnSymbols = map(columnSymbol, zip(*seqChunks))

        begs = map(coordinateString, [c+1 for c in coordinates])
        begWidth = maxlen(begs)
        for i, s in enumerate(seqChunks): coordinates[i] += countOfNonGaps(s)
        ends = map(coordinateString, coordinates)
        endWidth = maxlen(ends)

        for n, b, s, e in zip(names, begs, seqChunks, ends):
            print '%-*s' % (nameWidth, n),
            print '%*s' % (begWidth, b),
            print ''.join(htmlSpans(s, spans)),
            print '%*s' % (endWidth, e)
        print ' ' * nameWidth, ' ' * begWidth, ''.join(columnSymbols)
        print
    print '</pre>'

print '</body></html>'
