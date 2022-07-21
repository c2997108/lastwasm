// Author: Martin C. Frith 2013
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_last_splitter.hh"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <stdexcept>

// Defines an ordering, for sorting.
static bool less(const cbrc::UnsplitAlignment& a,
		 const cbrc::UnsplitAlignment& b) {
  int qnameCmp = strcmp(a.qname, b.qname);
  if (qnameCmp  != 0        ) return qnameCmp  < 0;
  if (a.qstart  != b.qstart ) return a.qstart  < b.qstart;
  if (a.qend    != b.qend   ) return a.qend    < b.qend;
  if (a.qstrand != b.qstrand) return a.qstrand < b.qstrand;
  int qalignCmp = strcmp(a.qalign, b.qalign);
  if (qalignCmp != 0        ) return qalignCmp < 0;
  int ralignCmp = strcmp(a.ralign, b.ralign);
  if (ralignCmp != 0        ) return ralignCmp < 0;
  return a.linesBeg < b.linesBeg;  // stabilizes the sort
}

static int printSense(char *out, double senseStrandLogOdds) {
  double b = senseStrandLogOdds / log(2.0);
  if (b < 0.1 && b > -0.1) b = 0;
  else if (b > 10) b = floor(b + 0.5);
  else if (b < -10) b = ceil(b - 0.5);
  int precision = (b < 10 && b > -10) ? 2 : 3;
  return sprintf(out, " sense=%.*g", precision, b);
}

namespace mcf {

void LastSplitter::doOneAlignmentPart(const LastSplitOptions &opts,
				      const cbrc::SplitAlignerParams &params,
				      bool isAlreadySplit,
				      const cbrc::UnsplitAlignment &a,
				      unsigned numOfParts, unsigned partNum,
				      unsigned alnNum,
				      unsigned qSliceBeg, unsigned qSliceEnd,
				      bool isSenseStrand,
				      double senseStrandLogOdds) {
  if (qSliceBeg >= a.qend || qSliceEnd <= a.qstart) {
    return;  // this can happen for spliced alignment!
  }

  unsigned qSliceBegTrimmed = qSliceBeg;
  unsigned qSliceEndTrimmed = qSliceEnd;
  unsigned alnBeg, alnEnd;
  cbrc::mafSliceBeg(a.ralign, a.qalign, a.qstart, qSliceBegTrimmed, alnBeg);
  cbrc::mafSliceEnd(a.ralign, a.qalign, a.qend,   qSliceEndTrimmed, alnEnd);

  if (qSliceBegTrimmed >= qSliceEndTrimmed) {
    return;  // I think this can happen for spliced alignment
  }

  int score =
    sa.segmentScore(alnNum, qSliceBeg, qSliceEnd) -
    sa.segmentScore(alnNum, qSliceBeg, qSliceBegTrimmed) -
    sa.segmentScore(alnNum, qSliceEndTrimmed, qSliceEnd);
  if (score < opts.score) return;

  std::vector<double> p;
  if (opts.direction != 0) {
    p = sa.marginalProbs(qSliceBegTrimmed, alnNum, alnBeg, alnEnd);
  }
  std::vector<double> pRev;
  if (opts.direction != 1) {
    sa.flipSpliceSignals(params);
    pRev = sa.marginalProbs(qSliceBegTrimmed, alnNum, alnBeg, alnEnd);
    sa.flipSpliceSignals(params);
  }
  if (opts.direction == 0) p.swap(pRev);
  if (opts.direction == 2) {
    double reverseProb = 1 / (1 + exp(senseStrandLogOdds));
    // the exp might overflow to inf, but that should be OK
    double forwardProb = 1 - reverseProb;
    for (unsigned i = 0; i < p.size(); ++i) {
      p[i] = forwardProb * p[i] + reverseProb * pRev[i];
    }
  }

  assert(!p.empty());
  double mismap = 1.0 - *std::max_element(p.begin(), p.end());
  mismap = std::max(mismap, 1e-10);
  if (mismap > opts.mismap) return;
  int mismapPrecision = 3;

  std::vector<char> slice;
  size_t lineLen = cbrc::mafSlice(slice, a, alnBeg, alnEnd, &p[0]);
  const char *sliceBeg = &slice[0];
  const char *sliceEnd = sliceBeg + slice.size();
  const char *pLine = sliceEnd - lineLen;
  const char *secondLastLine = pLine - lineLen;

  if (isAlreadySplit && secondLastLine[0] == 'p') {
    size_t backToSeq = alnEnd - alnBeg + 1;
    mismap = cbrc::pLinesToErrorProb(pLine - backToSeq, sliceEnd - backToSeq);
    if (mismap > opts.mismap) return;
    mismapPrecision = 2;
  }

  bool isLastalProbs = (*(secondLastLine - isAlreadySplit * lineLen) == 'p');
  if (opts.format == 'm' || (opts.format == 0 && !isLastalProbs)) {
    while (*(sliceEnd - lineLen) == 'p') sliceEnd -= lineLen;
  }

  const char *aLineOld = a.linesBeg[0];
  size_t aLineOldSize = a.linesBeg[1] - a.linesBeg[0] - 1;
  std::vector<char> aLine(aLineOldSize + 128);
  char *out = &aLine[0];
  if (opts.no_split && aLineOld[0] == 'a') {
    memcpy(out, aLineOld, aLineOldSize);
    out += aLineOldSize;
  } else {
    out += sprintf(out, "a score=%d", score);
  }
  out += sprintf(out, " mismap=%.*g", mismapPrecision, mismap);
  if (opts.direction == 2) out += printSense(out, senseStrandLogOdds);
  if (!opts.genome.empty() && !opts.no_split) {
    if (partNum > 0) {
      out = strcpy(out, isSenseStrand ? " acc=" : " don=") + 5;
      sa.spliceEndSignal(out, params, alnNum, qSliceBeg, isSenseStrand);
      out += 2;
    }
    if (partNum + 1 < numOfParts) {
      out = strcpy(out, isSenseStrand ? " don=" : " acc=") + 5;
      sa.spliceBegSignal(out, params, alnNum, qSliceEnd, isSenseStrand);
      out += 2;
    }
  }
  *out++ = '\n';

  outputText.insert(outputText.end(), &aLine[0], out);
  outputText.insert(outputText.end(), sliceBeg, sliceEnd);

  if (opts.no_split && a.linesEnd[-1][0] == 'c') {
    outputText.insert(outputText.end(), a.linesEnd[-1], a.linesEnd[0]);
    outputText.back() = '\n';
  }

  outputText.push_back('\n');
}

void LastSplitter::doOneQuery(const LastSplitOptions &opts,
			      const cbrc::SplitAlignerParams &params,
			      bool isAlreadySplit,
			      const cbrc::UnsplitAlignment *beg,
			      const cbrc::UnsplitAlignment *end) {
  if (opts.verbose) std::cerr << beg->qname << "\t" << (end - beg);
  sa.layout(params, beg, end);
  if (opts.verbose) std::cerr << "\tcells=" << sa.cellsPerDpMatrix();
  size_t bytes = sa.memory(params, !opts.no_split, opts.direction == 2);
  if (bytes > opts.bytes) {
    if (opts.verbose) std::cerr << "\n";
    std::cerr << "last-split: skipping sequence " << beg->qname
	      << " (" << bytes << " bytes)\n";
    return;
  }
  sa.initMatricesForOneQuery(params);

  long viterbiScore = LONG_MIN;
  long viterbiScoreRev = LONG_MIN;
  std::vector<unsigned> alnNums;
  std::vector<unsigned> queryBegs;
  std::vector<unsigned> queryEnds;

  if (opts.direction != 0) {
    sa.forwardBackward(params);
  }
  if (opts.direction != 1) {
    sa.flipSpliceSignals(params);
    sa.forwardBackward(params);
    sa.flipSpliceSignals(params);
  }

  double senseStrandLogOdds =
    (opts.direction == 2) ? sa.spliceSignalStrandLogOdds() : 0;

  if (opts.no_split) {
    if (opts.verbose) std::cerr << "\n";
    unsigned numOfParts = end - beg;
    for (unsigned i = 0; i < numOfParts; ++i) {
      doOneAlignmentPart(opts, params, isAlreadySplit, beg[i],
			 numOfParts, i, i, beg[i].qstart, beg[i].qend,
			 true, senseStrandLogOdds);
    }
  } else {
    if (opts.direction != 0) {
      viterbiScore = sa.viterbi(params);
      if (opts.verbose) std::cerr << "\t" << viterbiScore;
    }
    if (opts.direction != 1) {
      sa.flipSpliceSignals(params);
      viterbiScoreRev = sa.viterbi(params);
      sa.flipSpliceSignals(params);
      if (opts.verbose) std::cerr << "\t" << viterbiScoreRev;
    }
    bool isSenseStrand = (viterbiScore >= viterbiScoreRev);
    if (isSenseStrand) {
      sa.traceBack(params, viterbiScore, alnNums, queryBegs, queryEnds);
    } else {
      sa.flipSpliceSignals(params);
      sa.traceBack(params, viterbiScoreRev, alnNums, queryBegs, queryEnds);
      sa.flipSpliceSignals(params);
    }
    std::reverse(alnNums.begin(), alnNums.end());
    std::reverse(queryBegs.begin(), queryBegs.end());
    std::reverse(queryEnds.begin(), queryEnds.end());

    if (opts.verbose) std::cerr << "\n";
    unsigned numOfParts = alnNums.size();
    for (unsigned k = 0; k < numOfParts; ++k) {
      unsigned i = alnNums[k];
      doOneAlignmentPart(opts, params, isAlreadySplit, beg[i],
			 numOfParts, k, i, queryBegs[k], queryEnds[k],
			 isSenseStrand, senseStrandLogOdds);
    }
  }
}

void LastSplitter::split(const LastSplitOptions &opts,
			 const cbrc::SplitAlignerParams &params,
			 bool isAlreadySplit) {
  sort(mafs.begin(), mafs.end(), less);

  const cbrc::UnsplitAlignment *beg = mafs.empty() ? 0 : &mafs[0];
  const cbrc::UnsplitAlignment *end = beg + mafs.size();
  const cbrc::UnsplitAlignment *mid = beg;
  size_t qendMax = 0;

  while (mid < end) {
    if (mid->qend > qendMax) qendMax = mid->qend;
    ++mid;
    if (mid == end || strcmp(mid->qname, beg->qname) != 0 ||
	(mid->qstart >= qendMax && !opts.isSplicedAlignment)) {
      doOneQuery(opts, params, isAlreadySplit, beg, mid);
      beg = mid;
      qendMax = 0;
    }
  }

  mafs.clear();
}

void setLastSplitParams(cbrc::SplitAlignerParams &params,
			const LastSplitOptions &opts,
			const std::vector< std::vector<int> > &scoreMatrix,
			const char *rowNames, const char *colNames,
			int delOpenCost, int delGrowCost,
			int insOpenCost, int insGrowCost,
			double scale, double genomeSize, int sequenceFormat) {
  int restartCost =
    opts.isSplicedAlignment ? -(INT_MIN/2) : opts.score - 1;

  double jumpProb = opts.isSplicedAlignment
    ? opts.trans / (2 * genomeSize)  // 2 strands
    : 0.0;

  int jumpCost = (jumpProb > 0.0) ?
    -params.scoreFromProb(jumpProb, scale) : -(INT_MIN/2);

  if (sequenceFormat == 2 || sequenceFormat == 4 || sequenceFormat == 5) {
    throw std::runtime_error("unsupported last-split Q format");
  }

  int qualityOffset =
    (sequenceFormat == 1) ? 33 : (sequenceFormat == 3) ? 64 : 0;

  params.setParams(-delOpenCost, -delGrowCost, -insOpenCost, -insGrowCost,
		   -jumpCost, -restartCost, scale, qualityOffset);
  double splicePrior = opts.isSplicedAlignment ? opts.cis : 0.0;
  params.setSpliceParams(splicePrior, opts.mean, opts.sdev);
  params.setScoreMat(scoreMatrix, rowNames, colNames);
  params.setSpliceSignals();
  if (!opts.genome.empty()) params.readGenome(opts.genome);
}

}
