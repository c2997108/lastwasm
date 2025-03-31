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
  int rnameCmp = strcmp(a.rname, b.rname);
  if (rnameCmp  != 0        ) return rnameCmp  < 0;
  return a.linesBeg < b.linesBeg;  // stabilizes the sort
}

static int whichPairedSequence(const char *name) {
  size_t s = strlen(name);
  if (s > 1 && name[s-2] == '/') {
    if (name[s-1] == '1') return 1;
    if (name[s-1] == '2') return 2;
  }
  throw std::runtime_error("query sequence name should end in /1 or /2");
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

static void doOneSlice(SliceData &sd, cbrc::AlignmentPart &ap,
		       const cbrc::SplitAligner &sa,
		       const cbrc::UnsplitAlignment &a) {
  sd.score = -1;

  if (ap.queryBeg >= a.qend || ap.queryEnd <= a.qstart) {
    return;  // this can happen for spliced alignment!
  }

  unsigned qSliceBegOld = ap.queryBeg;
  unsigned qSliceEndOld = ap.queryEnd;
  cbrc::mafSliceBeg(a.ralign, a.qalign, a.qstart, ap.queryBeg, sd.alnBeg);
  cbrc::mafSliceEnd(a.ralign, a.qalign, a.qend,   ap.queryEnd, sd.alnEnd);

  if (ap.queryBeg >= ap.queryEnd) {
    return;  // I think this can happen for spliced alignment
  }

  sd.score =
    sa.segmentScore(ap.alnNum, qSliceBegOld, qSliceEndOld) -
    sa.segmentScore(ap.alnNum, qSliceBegOld, ap.queryBeg) -
    sa.segmentScore(ap.alnNum, ap.queryEnd, qSliceEndOld);
}

void LastSplitter::doOneAlignmentPart(const LastSplitOptions &opts,
				      const cbrc::SplitAlignerParams &params,
				      bool isAlreadySplit,
				      const cbrc::UnsplitAlignment &a,
				      unsigned numOfParts, unsigned partNum,
				      const SliceData &sd,
				      cbrc::AlignmentPart ap,
				      int rnaStrand, bool isSenseStrand,
				      double senseStrandLogOdds) {
  if (sd.score < opts.score) return;

  std::vector<double> columnProbabilities;
  size_t alnLen = sd.alnEnd - sd.alnBeg;
  columnProbabilities.resize(alnLen * (rnaStrand/2 + 1));
  double *probs = columnProbabilities.data();
  double *probsRev = probs + alnLen * rnaStrand/2;

  if (rnaStrand != 0) {
    sa.marginalProbs(probs, ap.queryBeg, ap.alnNum, sd.alnBeg, sd.alnEnd);
  }
  if (rnaStrand != 1) {
    sa.flipSpliceSignals(params);
    sa.marginalProbs(probsRev, ap.queryBeg, ap.alnNum, sd.alnBeg, sd.alnEnd);
    sa.flipSpliceSignals(params);
  }
  if (rnaStrand == 2) {
    double reverseProb = 1 / (1 + exp(senseStrandLogOdds));
    // the exp might overflow to inf, but that should be OK
    double forwardProb = 1 - reverseProb;
    for (size_t i = 0; i < alnLen; ++i) {
      probs[i] = forwardProb * probs[i] + reverseProb * probsRev[i];
    }
  }

  assert(alnLen > 0);
  double mismap = 1.0 - *std::max_element(probs, probs + alnLen);
  mismap = std::max(mismap, 1e-10);
  if (mismap > opts.mismap) return;

  bool isSplitProbs = (isAlreadySplit && a.linesEnd[-1][0] == 'p');
  bool isAlignProbs = (a.linesEnd[-1-isAlreadySplit][0] == 'p');
  int format = opts.format ? opts.format : "mM"[isAlignProbs];
  int mismapPrecision = 3 - isSplitProbs;

  bool isCopyFirstLine = (opts.no_split && a.linesBeg[0][0] == 'a');
  size_t firstLineSize = a.linesBeg[1] - a.linesBeg[0] - 1;
  size_t aLineSpace = firstLineSize * isCopyFirstLine + 128;

  std::vector<char> slice;
  size_t lineLen = cbrc::mafSlice(slice, a, sd.alnBeg, sd.alnEnd, probs);
  const char *sliceBeg = &slice[0];
  const char *sliceEnd = sliceBeg + slice.size();

  if (isSplitProbs) {
    size_t off = sd.alnEnd - sd.alnBeg + 1;
    mismap = cbrc::pLinesToErrorProb(sliceEnd - lineLen - off, sliceEnd - off);
    if (mismap > opts.mismap) return;
  }

  if (format == 'm') {
    while (*(sliceEnd - lineLen) == 'p') sliceEnd -= lineLen;
  }

  std::vector<char> aLine(aLineSpace);
  char *out = &aLine[0];
  if (isCopyFirstLine) {
    memcpy(out, a.linesBeg[0], firstLineSize);
    out += firstLineSize;
  } else {
    out += sprintf(out, "a score=%d", sd.score);
  }
  out += sprintf(out, " mismap=%.*g", mismapPrecision, mismap);
  if (rnaStrand == 2) out += printSense(out, senseStrandLogOdds);
  if (!opts.genome.empty() && !opts.no_split) {
    if (partNum > 0) {
      out = strcpy(out, isSenseStrand ? " acc=" : " don=") + 5;
      sa.spliceEndSignal(out, params, ap.alnNum, ap.queryBeg, isSenseStrand);
      out += 2;
    }
    if (partNum + 1 < numOfParts) {
      out = strcpy(out, isSenseStrand ? " don=" : " acc=") + 5;
      sa.spliceBegSignal(out, params, ap.alnNum, ap.queryEnd, isSenseStrand);
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
  int rnaStrand = opts.direction;
  if (rnaStrand == 3) rnaStrand = 2 - whichPairedSequence(beg->qname);
  if (rnaStrand == 4) rnaStrand = whichPairedSequence(beg->qname) - 1;

  sa.layout(params, beg, end);
  if (opts.verbose) std::cerr << "\tcells=" << sa.cellsPerDpMatrix();
  size_t bytes = sa.memory(params, rnaStrand == 2);
  if (bytes > opts.bytes) {
    if (opts.verbose) std::cerr << "\n";
    std::cerr << "last-split: skipping sequence " << beg->qname
	      << " (" << bytes << " bytes)\n";
    return;
  }
  sa.initMatricesForOneQuery(params, rnaStrand == 2);

  long viterbiScore = LONG_MIN;
  long viterbiScoreRev = LONG_MIN;
  std::vector<cbrc::AlignmentPart> alignmentParts;

  if (opts.no_split) {
    unsigned numOfParts = end - beg;
    alignmentParts.resize(numOfParts);
    for (unsigned i = 0; i < numOfParts; ++i) {
      alignmentParts[i].alnNum = i;
      alignmentParts[i].queryBeg = beg[i].qstart;
      alignmentParts[i].queryEnd = beg[i].qend;
    }
  } else {
    if (rnaStrand != 0) {
      viterbiScore = sa.viterbi(params);
      if (opts.verbose) std::cerr << "\t" << viterbiScore;
    }
    if (rnaStrand != 1) {
      sa.flipSpliceSignals(params);
      viterbiScoreRev = sa.viterbi(params);
      sa.flipSpliceSignals(params);
      if (opts.verbose) std::cerr << "\t" << viterbiScoreRev;
    }
    if (viterbiScore >= viterbiScoreRev) {
      sa.traceBack(params, viterbiScore, alignmentParts);
    } else {
      sa.flipSpliceSignals(params);
      sa.traceBack(params, viterbiScoreRev, alignmentParts);
      sa.flipSpliceSignals(params);
    }
    std::reverse(alignmentParts.begin(), alignmentParts.end());
  }

  unsigned numOfParts = alignmentParts.size();
  slices.resize(numOfParts);
  for (unsigned k = 0; k < numOfParts; ++k) {
    unsigned i = alignmentParts[k].alnNum;
    doOneSlice(slices[k], alignmentParts[k], sa, beg[i]);
  }

  sa.exponentiateScores(params);

  if (rnaStrand != 0) {
    sa.forwardBackward(params);
  }
  if (rnaStrand != 1) {
    sa.flipSpliceSignals(params);
    sa.forwardBackward(params);
    sa.flipSpliceSignals(params);
  }

  bool isSenseStrand = (viterbiScore >= viterbiScoreRev);

  double senseStrandLogOdds =
    (rnaStrand == 2) ? sa.spliceSignalStrandLogOdds() : 0;

  if (opts.verbose) std::cerr << "\n";

  for (unsigned k = 0; k < numOfParts; ++k) {
    unsigned i = alignmentParts[k].alnNum;
    doOneAlignmentPart(opts, params, isAlreadySplit, beg[i], numOfParts, k,
		       slices[k], alignmentParts[k],
		       rnaStrand, isSenseStrand, senseStrandLogOdds);
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
      if (opts.verbose) std::cerr << beg->qname << '\t' << beg->qstart << '\t'
				  << qendMax << '\t' << (mid - beg);
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
			int scoreMatrixStrand,
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
  params.setScoreMat(scoreMatrix, rowNames, colNames, scoreMatrixStrand);
  params.setSpliceSignals();
  if (!opts.genome.empty()) params.readGenome(opts.genome);
}

}
