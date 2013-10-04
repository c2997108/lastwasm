// Copyright 2013 Martin C. Frith

#include "cbrc_unsplit_alignment.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <numeric>  // accumulate
#include <sstream>
#include <stdexcept>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

namespace cbrc {
    
static const char *strBeg(const std::string& s) {
  return s.c_str();
}

static const char *strEnd(const std::string& s) {
  return s.c_str() + s.size();
}

void flipMafStrands(std::vector<std::string>& maf) {
  std::string a, b, c, d;
  for (unsigned i = 0; i < maf.size(); ++i) {
    std::string& line = maf[i];
    std::istringstream iss(line);
    if (line[0] == 's') {
      unsigned x, y, z;
      iss >> a >> b >> x >> y >> c >> z >> d;
      if (!iss) err("bad MAF line: " + line);
      x = z - x - y;
      reverse(d.begin(), d.end());
      std::ostringstream oss;
      oss << a << ' ' << b << ' ' << x << ' ' << y << ' '
	  << c << ' ' << z << ' ' << d;
      line = oss.str();
    } else if (line[0] == 'q') {
      iss >> a >> b >> c;
      if (!iss) err("bad MAF line: " + line);
      reverse(c.begin(), c.end());
      line = a + ' ' + b + ' ' + c;
    } else if (line[0] == 'p') {
      iss >> a >> b;
      if (!iss) err("bad MAF line: " + line);
      reverse(b.begin(), b.end());
      line = a + ' ' + b;
    }
  }
}

static void canonicalizeMafStrands(std::vector<std::string>& maf) {
  unsigned s = 0;
  for (unsigned i = 0; i < maf.size(); ++i) {
    std::string& line = maf[i];
    if (line[0] == 's') ++s;
    if (s == 2) {
      std::istringstream iss(line);
      std::string a, b, c, d, e;
      iss >> a >> b >> c >> d >> e;
      if (!iss) err("bad MAF line: " + line);
      if (e == "-") flipMafStrands(maf);
      return;
    }
  }
  err("bad MAF data");
}

void UnsplitAlignment::init() {
  canonicalizeMafStrands(lines);

  std::string a, b, c;
  unsigned s = 0;
  for (unsigned i = 0; i < lines.size(); ++i) {
    std::string& line = lines[i];
    std::istringstream iss(line);
    if (line[0] == 'a') {
      while (iss >> a >> std::ws && getline(iss, b, '='))
	if (b == "score" && iss >> score)
	  break;
    } else if (line[0] == 's') {
      ++s;
      unsigned len;
      if (s == 1) {
	iss >> a >> rname >> rstart >> len >> b >> c >> ralign;
	rend = rstart + len;
      } else if (s == 2) {
	iss >> a >> qname >> qstart >> len >> qstrand >> qfullend >> qalign;
	qend = qstart + len;
      }
    } else if (line[0] == 'q') {
      if (s == 1)
        err("I can't handle quality data for the genomic sequence");
      if (s == 2)
        iss >> a >> b >> qQual;
    }
    if (!iss) err("bad MAF line: " + line);
  }
}

static unsigned seqPosFromAlnPos(unsigned alnPos, const std::string& aln) {
  return alnPos - count(aln.begin(), aln.begin() + alnPos, '-');
}

std::vector<std::string> mafSlice(const std::vector<std::string>& maf,
				  unsigned alnBeg, unsigned alnEnd) {
  std::vector<std::string> out;
  std::string a, b, c, d, e, f;
  for (unsigned i = 0; i < maf.size(); ++i) {
    const std::string& line = maf[i];
    std::istringstream iss(line);
    if (line[0] == 's') {
      unsigned x;
      iss >> a >> b >> x >> c >> d >> e >> f;
      unsigned beg = x + seqPosFromAlnPos(alnBeg, f);
      unsigned end = x + seqPosFromAlnPos(alnEnd, f);
      unsigned len = end - beg;
      std::ostringstream oss;
      oss << a << ' ' << b << ' ' << beg << ' ' << len << ' '
          << d << ' ' << e << ' ' << f.substr(alnBeg, alnEnd - alnBeg);
      out.push_back(oss.str());
    } else if (line[0] == 'q') {
      iss >> a >> b >> c;
      out.push_back(a + ' ' + b + ' ' + c.substr(alnBeg, alnEnd - alnBeg));
    } else if (line[0] == 'p') {
      iss >> a >> b;
      out.push_back(a + ' ' + b.substr(alnBeg, alnEnd - alnBeg));
    }
  }
  return out;
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

void mafSliceBeg(const std::string& rAln, const std::string& qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg) {
  if (qSliceBeg < qBeg) {
    qSliceBeg = qBeg;
    alnBeg = 0;
  } else {
    alnBeg = nthBasePrefix(strBeg(qAln), qSliceBeg - qBeg);
  }
  unsigned numInserts = nthBasePrefix(strBeg(rAln) + alnBeg, 0);
  alnBeg += numInserts;
  qSliceBeg += numInserts;
}

void mafSliceEnd(const std::string& rAln, const std::string& qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd) {
  if (qSliceEnd > qEnd) {
    qSliceEnd = qEnd;
    alnEnd = qAln.size();
  } else {
    alnEnd = qAln.size() - nthBaseSuffix(strEnd(qAln), qEnd - qSliceEnd);
  }
  unsigned numInserts = nthBaseSuffix(strBeg(rAln) + alnEnd, 0);
  alnEnd -= numInserts;
  qSliceEnd -= numInserts;
}

static std::vector<unsigned>
sLineFieldWidths(const std::vector<std::string>& maf) {
  std::vector<unsigned> widths;
  for (unsigned i = 0; i < maf.size(); ++i) {
    const char* p = maf[i].c_str();
    if (*p != 's') continue;
    for (unsigned j = 0; *p; ++j) {
      const char* pOld = p;
      while (std::isgraph(*p)) ++p;
      unsigned width = p - pOld;
      if (widths.size() <= j) widths.push_back(width);
      else widths[j] = std::max(widths[j], width);
      while (std::isspace(*p)) ++p;
    }
  }
  return widths;
}

// Copy the next field of src to dest, left-justified
static void sprintLeft(char*& dest, const char*& src, unsigned width) {
  while (std::isspace(*src)) ++src;
  const char* s = src;
  while (std::isgraph(*src)) *dest++ = *src++;
  unsigned w = src - s;
  while (w++ < width) *dest++ = ' ';
  ++dest;
}

// Copy the next field of src to dest, right-justified
static void sprintRight(char*& dest, const char*& src, unsigned width) {
  while (std::isspace(*src)) ++src;
  const char* s = src;
  while (std::isgraph(*s)) ++s;
  unsigned w = s - src;
  while (w++ < width) *dest++ = ' ';
  while (std::isgraph(*src)) *dest++ = *src++;
  ++dest;
}

void printMaf(const std::vector<std::string>& maf) {
  std::vector<unsigned> w = sLineFieldWidths(maf);
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

// Probability -> phred score in fastq-sanger ASCII representation
static char asciiFromProb(double probRight) {
  double probWrong = 1 - probRight;
  double e = std::max(probWrong, 1e-10);  // avoid overflow errors
  int s = std::floor(-10 * std::log10(e));  // phred score, rounded down
  return std::min(s + 33, 126);
}

std::string pLineFromProbs(const std::vector<double>& p) {
  std::string s(2 + p.size(), ' ');
  s[0] = 'p';
  transform(p.begin(), p.end(), s.begin() + 2, asciiFromProb);
  return s;
}

}
