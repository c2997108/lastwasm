// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

#include "last-pair-probs.hh"

#include <algorithm>
#include <cctype>  // isdigit, etc
#include <cerrno>
#include <cmath>
#include <cstdlib>  // atof, atol
#include <cstring>  // strncmp
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <limits.h>
#include <cfloat>
#include <map>
#include <stddef.h>  // size_t

struct Alignment {
  std::string genomeStrand;
  long c;
  long rSize;
  double scaledScore;
  std::vector<std::string> text;
  std::string qName;
  bool operator<( const Alignment& aln ) const {
    return genomeStrand < aln.genomeStrand;
  }
};

typedef std::multimap<std::string, Alignment> MMAP;

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

static const char *readLong(const char *c, long &x) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  // this doesn't read negative numbers:
  if (!isDigit(*c)) return 0;
  long z = *c++ - '0';
  while (isDigit(*c)) {
    if (z > LONG_MAX / 10) return 0;
    z *= 10;
    long digit = *c++ - '0';
    if (z > LONG_MAX - digit) return 0;
    z += digit;
  }
  x = z;
  return c;
}

static const char *readDouble(const char *c, double &x) {
  if (!c) return 0;
  errno = 0;
  char *e;
  double z = std::strtod(c, &e);
  if (e == c || errno == ERANGE) return 0;
  x = z;
  return e;
}

static const char *readChar(const char *c, char &d) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  if (*c == 0) return 0;
  d = *c++;
  return c;
}

static const char *readWord(const char *c, std::string &s) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  const char *e = c;
  while (isGraph(*e)) ++e;
  if (e == c) return 0;
  s.assign(c, e);
  return e;
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

static std::istream& openIn(const std::string& fileName, std::ifstream& ifs) {
  ifs.open(fileName.c_str());
  if (!ifs) err("can't open file: " + fileName);
  return ifs;
}

static double logSumExp(std::vector<double>& numbers) {
  // Adds numbers, in log space, to avoid overflow.
  if (numbers.empty()) return -1.0e99;  // should be -inf
  const double m = *std::max_element(numbers.begin(), numbers.end());
  double s = 0.0;
  for (std::vector<double>::iterator itr = numbers.begin(); itr != numbers.end(); itr++) {
    s += std::exp(*itr - m);
  }
  return std::log(s) + m;
}

class AlignmentParameters {
  // Parses the score scale factor, minimum score, and genome size.
  double t;
  double e;
  long g;

 public:
  AlignmentParameters() {	// dummy values:
    this->t = -1.0;	// score scale factor
    this->e = -1.0;	// minimum score
    this->g = -1;		// genome size
  }

  void update(const std::string& line) {
    std::stringstream ss(line);
    std::string i;
    while (ss >> i) {
      if (this->t == -1.0 && i.substr(0,2) == "t=") {
        this->t = std::atof(i.substr(2).c_str());
        if (this->t <= 0) {
          err("t must be positive");
        }
      }
      if (this->e == -1.0 && i.substr(0,2) == "e=") {
        this->e = std::atof(i.substr(2).c_str());
        if (this->e <= 0) {
          err("e must be positive");
        }
      }
      if (this->g == -1.0 && i.substr(0,8) == "letters=") {
        this->g = std::atol(i.substr(8).c_str());
        if (this->g <= 0) {
          err("letters must be positive");
        }
      }
    }
  }

  bool isValid() const {
    return this->t != -1.0 && this->e != -1.0 && this->g != -1;
  }

  void validate() const {
    if (this->t == -1.0) err("I need a header line with t=");
    if (this->e == -1.0) err("I need a header line with e=");
    if (this->g == -1) err("I need a header line with letters=");
  }

  double tGet() const {
    return this->t;
  }

  double eGet() const {
    return this->e;
  }

  long gGet() const {
    return this->g;
  }
};

static void printAlignmentWithMismapProb(const Alignment& alignment,
                                         double prob, std::string suf) {
  const std::vector<std::string>& lines = alignment.text;
  const std::string& qName = alignment.qName;
  size_t qNameLen = qName.length();
  if (qNameLen >= 2 && qName[qNameLen-2] == '/')
    if (qName[qNameLen-1] == '1' || qName[qNameLen-1] == '2')
      suf = "";
  char p[32];
  sprintf(p, "%.3g", prob);
  if (lines.size() == 1) {	// we have tabular format
    const char *c = lines[0].c_str();
    const char *d = c;
    for (int i = 0; i < 7; ++i) d = skipWord(d);
    std::cout.write(c, d - c);
    std::cout << suf << d << '\t' << p << '\n';
  }
  else {	// we have MAF format
    std::cout << lines[0] << " mismap=" << p << '\n';
    std::string pad(suf.length(), ' ');	// spacer to keep the alignment of MAF lines
    size_t rNameEnd = alignment.genomeStrand.length() + 1;  // where to insert the spacer
    size_t qNameEnd = qNameLen + 2;	// where to insert the suffix
    unsigned s = 0;
    for (std::vector<std::string>::const_iterator itr = lines.begin()+1; itr != lines.end(); itr++) {
      const char *c = itr->c_str();
      if (*c == 's' || *c == 'q') {
        if (*c == 's') s++;
        if (s == 1) {
	  std::cout.write(c, rNameEnd);
	  std::cout << pad << (c + rNameEnd) << '\n';
        } else {
	  std::cout.write(c, qNameEnd);
	  std::cout << suf << (c + qNameEnd) << '\n';
        }
      } else if (*c == 'p') {
        std::cout.write(c, 1);
        std::cout << pad << (c + 1) << '\n';
      } else {
        std::cout << c << '\n';
      }
    }
    std::cout << '\n';	// each MAF block should end with a blank line
  }
}

static long headToHeadDistance(const Alignment& alignment1, const Alignment& alignment2) {
  // The 5'-to-5' distance between 2 alignments on opposite strands.
  long length = alignment1.c + alignment2.c;
  if (length > alignment1.rSize) {
    length -= alignment1.rSize;	// for circular chroms
  }
  return length;
}

static std::vector<double> conjointScores(const Alignment& aln1,
                                          const std::pair<MMAP::iterator, MMAP::iterator>& alns2,
                                          double fraglen, double inner, bool isRna) {
  std::vector<double> cScores;
  for (MMAP::iterator itr = alns2.first; itr != alns2.second; itr++) {
    const long length = headToHeadDistance(aln1, itr->second);
    if (isRna) {	// use a log-normal distribution
      if (length <= 0) continue;
      const double loglen = std::log((double)length);
      cScores.push_back(itr->second.scaledScore + inner * std::pow(loglen-fraglen, 2.0) - loglen);
    }
    else {    	// use a normal distribution
      if ((length > 0) != (fraglen > 0.0)) continue;	// ?
      cScores.push_back(itr->second.scaledScore + inner * std::pow(length-fraglen, 2.0));
    }
  }
  return cScores;
}

static std::vector<double> probForEachAlignment(const std::vector<Alignment>& alignments1,
                                                const std::vector<Alignment>& alignments2,
                                                const LastPairProbsOptions& opts) {
  std::vector<double> zs;
  std::vector<double> numbers;
  std::vector<Alignment>::const_iterator itr;
  for (itr = alignments2.begin(); itr != alignments2.end(); itr++) {
    numbers.push_back(itr->scaledScore);
  }
  const double x = opts.disjointScore + logSumExp(numbers);
  MMAP mapAlns2;
  for (itr = alignments2.begin(); itr != alignments2.end(); itr++) {
    mapAlns2.insert(std::pair<std::string, Alignment>(itr->genomeStrand, *itr));
  }

  std::vector<Alignment>::const_iterator aln1;
  for (aln1 = alignments1.begin(); aln1 != alignments1.end(); aln1++) {
    // get the items in alignments2 that have the same genomeStrand:
    std::pair<MMAP::iterator, MMAP::iterator> alns2 = mapAlns2.equal_range(aln1->genomeStrand);
    std::vector<double> cScores = conjointScores(*aln1, alns2, opts.fraglen, opts.inner, opts.rna);
    if (cScores.size()) {
      std::vector<double> xy;
      double y = opts.outer + logSumExp(cScores);
      xy.push_back(x);
      xy.push_back(y);
      zs.push_back(aln1->scaledScore + logSumExp(xy));
    } else {	// no items in alignments2 have the same genomeStrand
      zs.push_back(aln1->scaledScore + x);
    }
  }
  return zs;
}

static void printAlnsForOneRead(const std::vector<Alignment>& alignments1,
                                const std::vector<Alignment>& alignments2,
                                const LastPairProbsOptions& opts,
                                double maxMissingScore, const std::string suf) {
  std::vector<Alignment>::const_iterator itr;
  std::vector<double> zs;
  double w;

  if (alignments2.size()) {
    zs = probForEachAlignment(alignments1, alignments2, opts);
    double w0 = -DBL_MAX;
    for (itr = alignments2.begin(); itr != alignments2.end(); itr++) {
      w0 = (itr->scaledScore > w0) ? itr->scaledScore : w0;
    }
    w = maxMissingScore + w0;
  } else {
    for (itr = alignments1.begin(); itr != alignments1.end(); itr++) {
      zs.push_back(itr->scaledScore + opts.disjointScore);
    }
    w = maxMissingScore;
  }

  std::vector<double> zw0;
  zw0.push_back(logSumExp(zs));
  zw0.push_back(w);
  const double zw = logSumExp(zw0);
  
  std::vector<double>::iterator itrzs;
  for (itr = alignments1.begin(), itrzs = zs.begin(); itr != alignments1.end() && itrzs != zs.end(); itr++, itrzs++) {
    const double prob = 1.0 - std::exp(*itrzs - zw);
    if (prob <= opts.mismap) printAlignmentWithMismapProb(*itr, prob, suf);
  }
}

static void unambiguousFragmentLengths(const std::vector<Alignment>& alignments1,
                                       const std::vector<Alignment>& alignments2,
                                       std::vector<long>& lengths) {
  // Returns the fragment length implied by alignments of a pair of reads.
  long oldLength = LONG_MAX;
  std::vector<Alignment>::const_iterator itr1, itr2;
  for (itr1 = alignments1.begin(); itr1 != alignments1.end(); itr1++) {
    for (itr2 = alignments2.begin(); itr2 != alignments2.end(); itr2++) {
      if (itr1->genomeStrand == itr2->genomeStrand) {
        long newLength = headToHeadDistance(*itr1, *itr2);
        if (oldLength == LONG_MAX) {
          oldLength = newLength;
        } else if (newLength != oldLength) {
          return;  // the fragment length is ambiguous
        }
      }
    }
  }
  if (oldLength != LONG_MAX) {
    lengths.push_back(oldLength);
  }
}

static AlignmentParameters readHeaderOrDie(std::istream& lines) {
  std::string line;
  AlignmentParameters params;
  while (getline(lines, line)) {
    if (line.substr(0,1) == "#") {
      params.update(line);
      if (params.isValid()) {
        return params;
      }
    }
    else if (line.find_first_not_of(" ") != std::string::npos) {
      break;
    }
  }
  params.validate();	// die
  return params;			// dummy
}

static Alignment parseAlignment(double score, const std::string& rName,
                                long rStart, long rSpan, long rSize,
                                const std::string& qName, char qStrand,
                                const std::vector<std::string>& text,
                                char strand, double scale,
                                const std::set<std::string>& circularChroms) {
  const std::string genomeStrand = (qStrand == strand) ? rName + "+" : rName + "-";

  long c = -rStart;
  if (qStrand == '-') {
    c = rStart + rSpan;
    if (circularChroms.find(rName) != circularChroms.end() ||
        circularChroms.find(".") != circularChroms.end()) c += rSize;
  }

  const double scaledScore = score / scale;	// needed in 2nd pass
  Alignment parse;
  parse.genomeStrand = genomeStrand;
  parse.c = c;
  parse.rSize = rSize;
  parse.scaledScore = scaledScore;
  parse.text = text;
  parse.qName = qName;
  return parse;
}

static double parseMafScore(const char *aLine) {
  const char *c = aLine;
  while ((c = skipWord(c))) {
    c = skipSpace(c);
    if (std::strncmp(c, "score=", 6) == 0) {
      double score;
      c = readDouble(c + 6, score);
      if (!c) err("bad score");
      return score;
    }
  }
  err("missing score");
  return 0.0;	// dummy;
}

static Alignment parseMaf(const std::vector<std::string>& lines,
			  char strand, double scale,
			  const std::set<std::string>& circularChroms) {
  const double score = parseMafScore(lines[0].c_str());
  std::string rName, qName;
  char qStrand;
  long rStart, rSpan, rSize;
  unsigned n = 0;
  for (std::vector<std::string>::const_iterator itr = lines.begin(); itr != lines.end(); itr++) {
    const char *c = itr->c_str();
    if (*c == 's') {
      if (n == 0) {
	c = skipWord(c);
        c = readWord(c, rName);
        c = readLong(c, rStart);
        c = readLong(c, rSpan);
        c = skipWord(c);
	c = readLong(c, rSize);
      } else if (n == 1) {
	c = skipWord(c);
	c = readWord(c, qName);
	c = skipWord(c);
	c = skipWord(c);
	c = readChar(c, qStrand);
      }
      n++;
    }
    if (!c) err("bad MAF line: " + *itr);
  }
  if (n < 2) err("bad MAF");
  return parseAlignment(score, rName, rStart, rSpan, rSize, qName, qStrand,
                        lines, strand, scale, circularChroms);
}

static Alignment parseTab(const std::string& line,
			  char strand, double scale,
			  const std::set<std::string>& circularChroms) {
  std::vector<std::string> lines;
  lines.push_back(line);
  double score;
  std::string rName, qName;
  char qStrand;
  long rStart, rSpan, rSize;
  const char *c = line.c_str();
  c = readDouble(c, score);
  c = readWord(c, rName);
  c = readLong(c, rStart);
  c = readLong(c, rSpan);
  c = skipWord(c);
  c = readLong(c, rSize);
  c = readWord(c, qName);
  c = skipWord(c);
  c = skipWord(c);
  c = readChar(c, qStrand);
  if (!c) err("bad line: " + line);
  return parseAlignment(score, rName, rStart, rSpan, rSize, qName, qStrand,
                        lines, strand, scale, circularChroms);
}

static bool readBatch(std::istream& input,
		      char strand, const double scale,
		      const std::set<std::string>& circularChroms,
		      std::vector<Alignment>& alns) {
  // Yields alignment data from MAF or tabular format.
  alns.clear();
  std::vector<std::string> maf;
  std::string line;
  while (std::getline(input, line)) {
    if (std::isdigit(line[0])) {
      alns.push_back(parseTab(line, strand, scale, circularChroms));
    } else if (std::isalpha(line[0])) {
      maf.push_back(line);
    } else if (!std::isgraph(line[0])) {
      if (maf.size()) {
        alns.push_back(parseMaf(maf, strand, scale, circularChroms));
        maf.clear();
      }
    }
    else if (line.substr(0,8) == "# batch ") {
      break;
    }
  }
  if (maf.size()) {
    alns.push_back(parseMaf(maf, strand, scale, circularChroms));
  }

  return input;
}

static std::vector<long> readQueryPairs1pass(std::istream& in1, std::istream& in2,
                                             const double scale1, const double scale2,
                                             const std::set<std::string>& circularChroms) {
  std::vector<long> lengths;
  std::vector<Alignment> a1, a2;
  while (1) {
    bool ok1 = readBatch(in1, '+', scale1, circularChroms, a1);
    bool ok2 = readBatch(in2, '-', scale2, circularChroms, a2);
    unambiguousFragmentLengths(a1, a2, lengths);
    if (!ok1 || !ok2) break;
  }
  return lengths;
}

static void readQueryPairs2pass(std::istream& in1, std::istream& in2,
                                double scale1, double scale2,
                                const LastPairProbsOptions& opts) {
  std::vector<Alignment> a1, a2;
  while (1) {
    bool ok1 = readBatch(in1, '+', scale1, opts.circular, a1);
    bool ok2 = readBatch(in2, '-', scale2, opts.circular, a2);
    printAlnsForOneRead(a1, a2, opts, opts.maxMissingScore1, "/1");
    printAlnsForOneRead(a2, a1, opts, opts.maxMissingScore2, "/2");
    if (!ok1 || !ok2) break;
  }
}

static double myRound(double myFloat) {
  char buf[32];
  sprintf(buf, "%g", myFloat);
  return std::atof(buf);
}

static void estimateFragmentLengthDistribution(std::vector<long> lengths,
                                               LastPairProbsOptions& opts) {
  if (lengths.empty()) {
    err("can't estimate the distribution of distances");
  }

  // Define quartiles in the most naive way possible:
  std::sort(lengths.begin(), lengths.end());
  const size_t sampleSize = lengths.size();
  const long quartile1 = lengths[sampleSize / 4];
  const long quartile2 = lengths[sampleSize / 2];
  const long quartile3 = lengths[sampleSize * 3 / 4];

  std::cerr << opts.progName << ": distance sample size: " << sampleSize << "\n";
  std::cerr << opts.progName << ": distance quartiles: " << quartile1 << " " << quartile2 << " " << quartile3 << "\n";

  if (opts.rna && quartile1 <= 0) {
    err("too many distances <= 0");
  }

  const char* thing = (opts.rna) ? "ln[distance]" : "distance";

  if (!opts.isFraglen) {
    if (opts.rna) {
      opts.fraglen = myRound(std::log((double)quartile2));
    } else {
      opts.fraglen = double(quartile2);
    }
    std::cerr << opts.progName << ": estimated mean " << thing << ": " << opts.fraglen << "\n";
  }

  if (!opts.isSdev) {
    const double iqr = (opts.rna) ? 
      std::log((double)quartile3) - std::log((double)quartile1) : quartile3 - quartile1;
    // Normal Distribution: sdev = iqr / (2 * qnorm(0.75))
    opts.sdev = myRound(iqr / 1.34898);
    std::cerr << opts.progName << ": estimated standard deviation of " << thing << ": " << opts.sdev << "\n";
  }
}

static double safeLog(const double x) {
  if (x == 0.0) {
    return -1.0e99;
  } else {
    return std::log(x);
  }
}

static void calculateScorePieces(LastPairProbsOptions& opts,
                                 const AlignmentParameters& params1,
                                 const AlignmentParameters& params2) {
  if (opts.sdev == 0.0) {
    opts.outer = opts.rna ? opts.fraglen : 0.0;
    opts.inner = -1.0e99;
  }
  else {	// parameters for a Normal Distribution (of fragment lengths):
    const double pi = atan(1.0) * 4.0;
    opts.outer = -std::log(opts.sdev * std::sqrt(2.0 * pi));
    opts.inner = -1.0 / (2.0 * std::pow(opts.sdev, 2.0));
  }
  opts.outer += safeLog(1.0 - opts.disjoint);

  if (params1.gGet() != params2.gGet()) err("unequal genome sizes");
  // Multiply genome size by 2, because it has 2 strands:
  opts.disjointScore = safeLog(opts.disjoint) - std::log(2.0 * params1.gGet());

  // Max possible influence of an alignment just below the score threshold:
  double maxLogPrior = opts.outer;
  if (opts.rna) maxLogPrior += std::pow(opts.sdev, 2.0) / 2.0 - opts.fraglen;
  opts.maxMissingScore1 = (params1.eGet() - 1.0) / params1.tGet() + maxLogPrior;
  opts.maxMissingScore2 = (params2.eGet() - 1.0) / params2.tGet() + maxLogPrior;
}

void lastPairProbs(LastPairProbsOptions& opts) {
  std::string fileName1 = opts.inputFileNames[0];
  std::string fileName2 = opts.inputFileNames[1];

  if (!opts.isFraglen || !opts.isSdev) {
    std::ifstream inFile1, inFile2;
    std::istream& in1 = openIn(fileName1, inFile1);
    std::istream& in2 = openIn(fileName2, inFile2);
    std::vector<long> lengths = readQueryPairs1pass(in1, in2, 1.0, 1.0, opts.circular);
    estimateFragmentLengthDistribution(lengths, opts);
  }

  if (!opts.estdist) {
    std::ifstream inFile1, inFile2;
    std::istream& in1 = openIn(fileName1, inFile1);
    std::istream& in2 = openIn(fileName2, inFile2);
    AlignmentParameters params1 = readHeaderOrDie(in1);
    AlignmentParameters params2 = readHeaderOrDie(in2);
    calculateScorePieces(opts, params1, params2);
    std::cout << "# fraglen=" << opts.fraglen
              << " sdev=" << opts.sdev
              << " disjoint=" << opts.disjoint
              << " genome=" << params1.gGet() << "\n";
    readQueryPairs2pass(in1, in2, params1.tGet(), params2.tGet(), opts);
  }
}
