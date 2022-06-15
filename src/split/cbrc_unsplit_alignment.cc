// Copyright 2013, 2014 Martin C. Frith

#include "cbrc_unsplit_alignment.hh"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <iostream>
#include <numeric>  // accumulate
#include <stdexcept>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static bool isGraph(char c) {
  return c > ' ';  // faster than std::isgraph
}

static bool isSpace(char c) {
  return c > 0 && c <= ' ';  // faster than std::isspace
}

static bool isDigit(char c) {
  return c >= '0' && c <= '9';
}

namespace cbrc {
    
static const char *readUint(const char *c, unsigned &x) {
  if (!c) return 0;

  // faster than std::strtoul
  while (isSpace(*c)) ++c;
  if (!isDigit(*c)) return 0;
  unsigned z = *c++ - '0';
  while (isDigit(*c)) {
    if (z > UINT_MAX / 10) return 0;
    unsigned digit = *c++ - '0';
    z = z * 10 + digit;
    if (z < digit) return 0;
  }

  x = z;
  return c;
}

static const char *readChar(const char *c, char &d) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  if (*c == 0) return 0;
  d = *c++;
  return c;
}

static const char *skipWord(const char *c) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  const char *e = c;
  while (isGraph(*e)) ++e;
  if (e == c) return 0;
  return e;
}

static const char *skipSpace(const char *c) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  return c;
}

static const char *rskipSpace(const char *c) {
  while (isSpace(*(c-1))) --c;
  return c;
}

struct Complement {
  unsigned char c[UCHAR_MAX+1];
  Complement() {
    static const char x[] = "ACGTRYSWKMBDHVN";
    static const char y[] = "TGCAYRSWMKVHDBN";
    for (unsigned i = 0; i < UCHAR_MAX+1; ++i) c[i] = i;
    for (unsigned i = 0; x[i] && y[i]; ++i) {
      c[toupper(x[i])] = toupper(y[i]);
      c[tolower(x[i])] = tolower(y[i]);
    }
  }
  unsigned char operator()(unsigned char x) {
    return c[x];
  }
};
static Complement complement;

static bool isRevQryStrand(StringIt linesBeg, StringIt linesEnd,
			   unsigned rankOfQrySeq) {
  char strand = 0;
  unsigned s = 0;
  for (StringIt i = linesBeg; i < linesEnd; ++i) {
    const char *c = i->c_str();
    if (*c == 's' && ++s == rankOfQrySeq) {
      c = skipWord(c);
      c = skipWord(c);
      c = skipWord(c);
      c = skipWord(c);
      c = readChar(c, strand);
    }
  }
  return strand == '-';
}

void UnsplitAlignment::init(bool isTopSeqQuery) {
  const unsigned rankOfQrySeq = 2 - isTopSeqQuery;
  const unsigned rankOfRefSeq = isTopSeqQuery + 1;

  char refStrand, qryStrand;
  unsigned refSpan, qrySpan;
  unsigned refSeqLen, qrySeqLen;
  qQual = 0;  // in case the input lacks sequence quality data
  unsigned sLineCount = 0;

  bool isRev = isRevQryStrand(linesBeg, linesEnd, rankOfQrySeq);

  for (StringIt i = linesBeg; i < linesEnd; ++i) {
    const char *c = i->c_str();
    const char *lineEnd = c + i->size();
    const char *d, *e, *f;
    const char *g = rskipSpace(lineEnd);
    if (*c == 's') {
      ++sLineCount;
      unsigned start = 0;
      unsigned len = 0;
      unsigned seqLen = 0;
      char strand = 0;
      d = skipWord(c);
      d = skipSpace(d);
      e = skipWord(d);
      f = readUint(e, start);
      f = readUint(f, len);
      f = readChar(f, strand);
      f = readUint(f, seqLen);
      f = skipSpace(f);
      if (!f || f >= g) err("bad MAF line: " + *i);
      (*i)[e - c] = 0;  // write a terminator for the sequence name
      if (g < lineEnd) (*i)[g - c] = 0;  // trim trailing whitespace
      if (sLineCount == rankOfRefSeq) {
	rstart = start;
	refSpan = len;
	refStrand = strand;
	refSeqLen = seqLen;
	rname = d;
	ralign = f;
      } else if (sLineCount == rankOfQrySeq) {
	qstart = start;
	qrySpan = len;
	qryStrand = strand;
	qrySeqLen = seqLen;
	qname = d;
	qalign = f;
      }
      if (isRev) {
	std::string::iterator beg = i->begin() + (f - c);
	std::string::iterator end = i->begin() + (g - c);
	reverse(beg, end);
	transform(beg, end, beg, complement);
      }
    } else if (*c == 'q') {
      if (sLineCount != rankOfQrySeq) {
	err("I can only handle quality data for the query sequence");
      }
      f = skipSpace(skipWord(skipWord(c)));
      if (!f || f >= g) err("bad MAF line: " + *i);
      if (g < lineEnd) (*i)[g - c] = 0;  // trim trailing whitespace
      if (qryStrand == '-') {
	reverse(i->begin() + (f - c), i->begin() + (g - c));
      }
      qQual = f;
    }
  }

  if (sLineCount != 2) err("I need 2 sequences per alignment");

  if (qryStrand == '-') {
    rstart = refSeqLen - rstart - refSpan;
    qstart = qrySeqLen - qstart - qrySpan;
  }

  qstrand = (qryStrand != refStrand) * 2 + (qryStrand == '-');

  rend = rstart + refSpan;
  qend = qstart + qrySpan;
}

static unsigned seqPosFromAlnPos(unsigned alnPos, const char *aln) {
  return alnPos - std::count(aln, aln + alnPos, '-');
}

static unsigned nthBasePrefix(const char* sequenceWithGapsBeg, unsigned n) {
  for (unsigned i = 0; /* noop */; ++i)
    if (sequenceWithGapsBeg[i] != '-') {
      if (n > 0) --n;
      else return i;
    }
}

static unsigned nthBaseSuffix(const char *sequenceWithGapsEnd, unsigned n) {
  for (unsigned i = 0; /* noop */; ++i)
    if (*(sequenceWithGapsEnd - 1 - i) != '-') {
      if (n > 0) --n;
      else return i;
    }
}

void mafSliceBeg(const char* rAln, const char* qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg) {
  if (qSliceBeg < qBeg) {
    qSliceBeg = qBeg;
    alnBeg = 0;
  } else {
    alnBeg = nthBasePrefix(qAln, qSliceBeg - qBeg);
  }
  unsigned numInserts = nthBasePrefix(rAln + alnBeg, 0);
  alnBeg += numInserts;
  qSliceBeg += numInserts;
}

void mafSliceEnd(const char* rAln, const char* qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd) {
  unsigned alnLength = strlen(qAln);
  if (qSliceEnd > qEnd) {
    qSliceEnd = qEnd;
    alnEnd = alnLength;
  } else {
    alnEnd = alnLength - nthBaseSuffix(qAln + alnLength, qEnd - qSliceEnd);
  }
  unsigned numInserts = nthBaseSuffix(rAln + alnEnd, 0);
  alnEnd -= numInserts;
  qSliceEnd -= numInserts;
}

static const char *updateFieldWidth(const char *p, int &width) {
  const char *b = p;
  while (isGraph(*p)) ++p;
  if (p - b > width) width = p - b;
  while (isSpace(*p)) ++p;
  return p;
}

static std::vector<int>
sLineFieldWidths(const std::vector<std::string>& maf) {
  std::vector<int> widths;
  for (unsigned i = 0; i < maf.size(); ++i) {
    const char* p = maf[i].c_str();
    if (*p != 's') continue;
    for (unsigned j = 0; *p; ++j) {
      if (widths.size() <= j) widths.push_back(0);
      p = updateFieldWidth(p, widths[j]);
    }
  }
  return widths;
}

// Copy the next field of src to dest, left-justified
static void sprintLeft(char*& dest, const char*& src, int width) {
  const char *end = dest + width;
  while (isSpace(*src)) ++src;
  while (isGraph(*src)) *dest++ = *src++;
  while (dest < end) *dest++ = ' ';
  ++dest;
}

// Copy the next field of src to dest, right-justified
static void sprintRight(char*& dest, const char*& src, int width) {
  while (isSpace(*src)) ++src;
  const char* s = src;
  while (isGraph(*s)) ++s;
  int w = s - src;
  while (w++ < width) *dest++ = ' ';
  while (isGraph(*src)) *dest++ = *src++;
  ++dest;
}

// Probability -> phred score in fastq-sanger ASCII representation
static char asciiFromProb(double probRight) {
  double probWrong = 1 - probRight;
  double e = std::max(probWrong, 1e-10);  // avoid overflow errors
  int s = floor(-10 * log10(e));  // phred score, rounded down
  return std::min(s + 33, 126);
}

std::vector<std::string> mafSlice(const UnsplitAlignment &aln,
				  unsigned alnBeg, unsigned alnEnd,
				  const double *probs) {
  std::vector<std::string> out;

  unsigned alnLen = alnEnd - alnBeg;

  for (StringCi i = aln.linesBeg; i < aln.linesEnd; ++i) {
    const char *c = i->c_str();
    const char *d, *e, *f;
    if (*c == 's') {
      unsigned beg = 0;
      unsigned len = 0;
      d = skipWord(skipWord(c));
      e = d + 1;  // skip over the string terminator
      e = readUint(e, beg);
      e = readUint(e, len);
      f = skipWord(e);
      f = skipWord(f);
      f = skipSpace(f);
      unsigned begPos = seqPosFromAlnPos(alnBeg, f);
      unsigned endPos = seqPosFromAlnPos(alnEnd, f);
      unsigned newBeg = aln.isFlipped() ? beg + len - endPos : beg + begPos;
      unsigned newLen = endPos - begPos;
      char buffer[64];
      sprintf(buffer, " %u %u", newBeg, newLen);
      out.push_back(std::string(c, d) + buffer +
		    std::string(e, f) + std::string(f + alnBeg, f + alnEnd));
      if (aln.isFlipped()) {
	reverse(out.back().end() - alnLen, out.back().end());
	transform(out.back().end() - alnLen, out.back().end(),
		  out.back().end() - alnLen, complement);
      }
    } else if (*c == 'q') {
      d = skipSpace(skipWord(skipWord(c)));
      out.push_back(std::string(c, d) + std::string(d + alnBeg, d + alnEnd));
      if (aln.isFlipped()) {
	reverse(out.back().end() - alnLen, out.back().end());
      }
    } else if (*c == 'p') {
      d = skipSpace(skipWord(c));
      const char *g = rskipSpace(c + i->size());
      out.push_back(std::string(c, d) +
		    (aln.isFlipped()
		     ? std::string(g - alnEnd, g - alnBeg)
		     : std::string(d + alnBeg, d + alnEnd)));
    }
  }

  std::string s(2 + alnLen, ' ');
  s[0] = 'p';
  transform(probs, probs + alnLen, s.begin() + 2, asciiFromProb);
  if (aln.isFlipped()) reverse(s.begin() + 2, s.end());
  out.push_back(s);

  return out;
}

void printMaf(const std::vector<std::string>& maf) {
  std::vector<int> w = sLineFieldWidths(maf);
  unsigned lineLength = std::accumulate(w.begin(), w.end(), w.size());
  std::vector<char> line(lineLength, ' ');
  line[lineLength - 1] = '\n';

  for (unsigned i = 0; i < maf.size(); ++i) {
    const char* src = maf[i].c_str();
    char* dest = &line[0];
    if (*src == 's') {
      sprintLeft(dest, src, w[0]);
      sprintLeft(dest, src, w[1]);
      sprintRight(dest, src, w[2]);
      sprintRight(dest, src, w[3]);
      sprintLeft(dest, src, w[4]);
      sprintRight(dest, src, w[5]);
      sprintLeft(dest, src, w[6]);
      std::cout.write(&line[0], lineLength);
    } else if (*src == 'q') {
      sprintLeft(dest, src, w[0]);
      sprintLeft(dest, src, w[1] + w[2] + w[3] + w[4] + w[5] + 4);
      sprintLeft(dest, src, w[6]);
      std::cout.write(&line[0], lineLength);
    } else if (*src == 'p') {
      sprintLeft(dest, src, w[0] + w[1] + w[2] + w[3] + w[4] + w[5] + 5);
      sprintLeft(dest, src, w[6]);
      std::cout.write(&line[0], lineLength);
    } else {
      std::cout << src << '\n';
    }
  }

  std::cout << '\n';
}

double pLinesToErrorProb(const char *line1, const char *line2) {
  double maxGoodProb = 0;
  const char *i = skipSpace(skipWord(line1));
  const char *j = skipSpace(skipWord(line2));
  while (isGraph(*i) && isGraph(*j)) {
    double x = pow(0.1, (*i - 33) * 0.1);  // error probability
    double y = pow(0.1, (*j - 33) * 0.1);  // error probability
    double z = (1 - x) * (1 - y);  // probability that neither is an error
    if (z > maxGoodProb) maxGoodProb = z;
    ++i;
    ++j;
  }
  return 1 - maxGoodProb;  // minimum combined error probability
}

}
