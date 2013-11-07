// Copyright 2013 Martin C. Frith

#include "last-split.hh"

#include "cbrc_split_aligner.hh"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static std::istream& openIn(const std::string& fileName, std::ifstream& ifs) {
  if (fileName == "-") return std::cin;
  ifs.open(fileName.c_str());
  if (!ifs) err("can't open file: " + fileName);
  return ifs;
}

// Does the string start with the prefix?
static bool startsWith(const std::string& s, const char* prefix) {
  const char* t = s.c_str();
  for (;;) {
    if (*prefix == 0) return true;
    if (*prefix != *t) return false;
    ++t;
    ++prefix;
  }
}

// Does the string have no non-space characterts?
static bool isSpace(const std::string& s) {
  const char* t = s.c_str();
  for (;;) {
    if (*t == 0) return true;
    if (!std::isspace(*t)) return false;
    ++t;
  }
}

static int scoreFromProb(double prob, double scale) {
  return std::floor(scale * std::log(prob) + 0.5);
}

// Defines an ordering, for sorting.
static bool less(const cbrc::UnsplitAlignment& a,
		 const cbrc::UnsplitAlignment& b) {
  return
    a.qname != b.qname ? a.qname < b.qname :
    a.qstart != b.qstart ? a.qstart < b.qstart :
    a.qend != b.qend ? a.qend < b.qend :
    a.qstrand != b.qstrand ? a.qstrand < b.qstrand :
    a.qalign != b.qalign ? a.qalign < b.qalign :
    a.rname < b.rname;
}

static void doOneAlignmentPart(cbrc::SplitAligner& sa,
			       const cbrc::UnsplitAlignment& a,
			       unsigned alnNum,
			       unsigned qSliceBeg, unsigned qSliceEnd,
			       double forwardDirectionProb,
			       const LastSplitOptions& opts) {
  unsigned alnBeg, alnEnd;
  cbrc::mafSliceBeg(a.ralign, a.qalign, a.qstart, qSliceBeg, alnBeg);
  cbrc::mafSliceEnd(a.ralign, a.qalign, a.qend, qSliceEnd, alnEnd);

  int score = sa.segmentScore(alnNum, qSliceBeg, qSliceEnd);
  if (score < opts.score) return;

  std::vector<double> p;
  if (opts.direction != 0) {
    p = sa.marginalProbs(qSliceBeg, alnNum, alnBeg, alnEnd);
  }
  std::vector<double> pRev;
  if (opts.direction != 1) {
    sa.flipSpliceSignals();
    pRev = sa.marginalProbs(qSliceBeg, alnNum, alnBeg, alnEnd);
    sa.flipSpliceSignals();
  }
  if (opts.direction == 0) p.swap(pRev);
  if (opts.direction == 2) {
    double reverseDirectionProb = 1.0 - forwardDirectionProb;
    for (unsigned i = 0; i < p.size(); ++i) {
      p[i] = forwardDirectionProb * p[i] + reverseDirectionProb * pRev[i];
    }
  }

  assert(!p.empty());
  double mismap = 1.0 - *std::max_element(p.begin(), p.end());
  mismap = std::max(mismap, 1e-10);
  if (mismap > opts.mismap) return;

  std::cout << std::setprecision(3)
	    << "a score=" << score << " mismap=" << mismap << "\n"
	    << std::setprecision(6);
  std::vector<std::string> s = cbrc::mafSlice(a.lines, alnBeg, alnEnd);
  s.push_back(cbrc::pLineFromProbs(p));
  if (a.qstrand == "-") cbrc::flipMafStrands(s);
  cbrc::printMaf(s);
}

static void doOneQuery(std::vector<cbrc::UnsplitAlignment>::const_iterator beg,
		       std::vector<cbrc::UnsplitAlignment>::const_iterator end,
		       cbrc::SplitAligner& sa, const LastSplitOptions& opts) {
  if (opts.verbose) std::cerr << beg->qname;
  sa.initForOneQuery(beg, end);

  if (opts.direction != 0) {
    sa.forward();
    sa.backward();
  }
  if (opts.direction != 1) {
    sa.flipSpliceSignals();
    sa.forward();
    sa.backward();
    sa.flipSpliceSignals();
  }

  double forwardDirectionProb = -1;
  if (opts.direction == 2) {
    forwardDirectionProb = sa.spliceSignalStrandProb();
    if (opts.verbose) std::cerr << "\tforwardProb=" << forwardDirectionProb;
  }

  if (opts.no_split) {
    if (opts.verbose) std::cerr << "\n";
    for (unsigned i = 0; i < end - beg; ++i) {
      doOneAlignmentPart(sa, beg[i], i, 0, beg->qfullend,
			 forwardDirectionProb, opts);
    }
  } else {
    long viterbiScore = LONG_MIN;
    if (opts.direction != 0) {
      viterbiScore = sa.viterbi();
      if (opts.verbose) std::cerr << "\t" << viterbiScore;
    }
    long viterbiScoreRev = LONG_MIN;
    if (opts.direction != 1) {
      sa.flipSpliceSignals();
      viterbiScoreRev = sa.viterbi();
      sa.flipSpliceSignals();
      if (opts.verbose) std::cerr << "\t" << viterbiScoreRev;
    }
    std::vector<unsigned> alnNums;
    std::vector<unsigned> queryBegs;
    std::vector<unsigned> queryEnds;
    if (viterbiScore >= viterbiScoreRev) {
      sa.traceBack(viterbiScore, alnNums, queryBegs, queryEnds);
    } else {
      sa.flipSpliceSignals();
      sa.traceBack(viterbiScoreRev, alnNums, queryBegs, queryEnds);
      sa.flipSpliceSignals();
    }
    std::reverse(alnNums.begin(), alnNums.end());
    std::reverse(queryBegs.begin(), queryBegs.end());
    std::reverse(queryEnds.begin(), queryEnds.end());

    if (opts.verbose) std::cerr << "\n";
    for (unsigned k = 0; k < alnNums.size(); ++k) {
      unsigned i = alnNums[k];
      doOneAlignmentPart(sa, beg[i], i, queryBegs[k], queryEnds[k],
			 forwardDirectionProb, opts);
    }
  }
}

static void doOneBatch(std::vector<cbrc::UnsplitAlignment>& mafs,
                       cbrc::SplitAligner& sa, const LastSplitOptions& opts) {
  stable_sort(mafs.begin(), mafs.end(), less);
  std::vector<cbrc::UnsplitAlignment>::const_iterator b = mafs.begin();
  std::vector<cbrc::UnsplitAlignment>::const_iterator e = mafs.begin();
  while (e < mafs.end()) {
    ++e;
    if (e == mafs.end() || e->qname != b->qname) {
      doOneQuery(b, e, sa, opts);
      b = e;
    }
  }
}

static void printParameters(const LastSplitOptions& opts) {
  std::cout << std::setprecision(12) << "#"
	    << " m=" << opts.mismap
	    << " s=" << opts.score;
  if (opts.isSplicedAlignment) {
    std::cout << " d=" << opts.direction
	      << " c=" << opts.cis
	      << " t=" << opts.trans
	      << " M=" << opts.mean
	      << " S=" << opts.sdev;
  }
  std::cout << "\n" << std::setprecision(6);
}

static void addMaf(std::vector<cbrc::UnsplitAlignment>& mafs,
		   const std::vector<std::string>& maf) {
  if (maf.empty()) return;
  mafs.push_back(cbrc::UnsplitAlignment(maf));
}

void lastSplit(LastSplitOptions& opts) {
  cbrc::SplitAligner sa;
  std::vector< std::vector<int> > scoreMatrix;
  std::string rowNames, colNames;
  std::string line, word, name, key;
  int state = 0;
  int gapExistenceCost = -1;
  int gapExtensionCost = -1;
  int lastalScoreThreshold = -1;
  double scale = 0;
  double genomeSize = 0;
  std::vector<std::string> maf;
  std::vector<cbrc::UnsplitAlignment> mafs;

  for (unsigned i = 0; i < opts.inputFileNames.size(); ++i) {
    std::ifstream inFileStream;
    std::istream& input = openIn(opts.inputFileNames[i], inFileStream);
    while (getline(input, line)) {
      if (state == -1) {  // we are reading the score matrix within the header
	std::istringstream ls(line);
	std::vector<int> row;
	int score;
	ls >> word >> name;
	while (ls >> score) row.push_back(score);
	if (word == "#" && name.size() == 1 && !row.empty() && ls.eof()) {
	  rowNames.push_back(std::toupper(name[0]));
	  scoreMatrix.push_back(row);
	} else {
	  state = 0;
	}
      }
      if (state == 0) {  // we are reading the header
	std::istringstream ls(line);
	std::string names;
	ls >> word;
	while (ls >> name) {
	  if (name.size() == 1) names.push_back(std::toupper(name[0]));
	  else break;
	}
	if (word == "#" && !names.empty() && !ls && scoreMatrix.empty()) {
	  colNames = names;
	  state = -1;
	} else if (startsWith(line, "#")) {
	  std::istringstream ls(line);
	  while (ls >> word) {
	    std::istringstream ws(word);
	    getline(ws, key, '=');
	    if (key == "a") ws >> gapExistenceCost;
	    if (key == "b") ws >> gapExtensionCost;
	    if (key == "e") ws >> lastalScoreThreshold;
	    if (key == "t") ws >> scale;
	    if (key == "letters") ws >> genomeSize;
	  }
	} else if (!isSpace(line)) {
	  if (scoreMatrix.empty())
	    err("I need a header with score parameters");
	  if (gapExistenceCost < 0 || gapExtensionCost < 0 ||
	      lastalScoreThreshold < 0 || scale <= 0 || genomeSize <= 0)
	    err("can't read the header");
	  if (opts.score < 0)
	    opts.score = lastalScoreThreshold + scoreFromProb(1000, scale);
	  int restartCost =
	    opts.isSplicedAlignment ? -(INT_MIN/2) : opts.score - 1;
	  double jumpProb = opts.isSplicedAlignment
	    ? opts.trans / (2 * genomeSize)  // 2 strands
	    : 0.0;
	  int jumpCost =
	    (jumpProb > 0.0) ? -scoreFromProb(jumpProb, scale) : -(INT_MIN/2);
	  printParameters(opts);
	  sa.setParams(-gapExistenceCost, -gapExtensionCost,
		       -jumpCost, -restartCost, scale);
	  double splicePrior = opts.isSplicedAlignment ? opts.cis : 0.0;
	  sa.setSpliceParams(splicePrior, opts.mean, opts.sdev);
	  sa.setScoreMat(scoreMatrix, rowNames, colNames);
	  sa.setSpliceSignals();
	  if (!opts.genome.empty()) sa.readGenome(opts.genome);
	  sa.printParameters();
	  std::cout << "#\n";
	  state = 1;
	}
      }
      if (state == 1) {  // we are reading alignments
	if (startsWith(line, "# batch")) {
	  addMaf(mafs, maf);
	  maf.clear();
	  doOneBatch(mafs, sa, opts);
	  mafs.clear();
	} else if (isSpace(line)) {
	  addMaf(mafs, maf);
	  maf.clear();
	} else if (line[0] != '#') {
	  maf.push_back(line);
	}
      }
      if (startsWith(line, "#")) std::cout << line << "\n";
    }
  }
  addMaf(mafs, maf);
  doOneBatch(mafs, sa, opts);
}
