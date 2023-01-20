// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

// Functions that find gapless X-drop alignments between a sequence
// and a PSSM.

// These functions are analogous to those described in
// gaplessXdrop.hh.  Here, "pssm" replaces both "seq2" and "scorer".

#ifndef GAPLESS_PSSM_XDROP_HH
#define GAPLESS_PSSM_XDROP_HH

#include "ScoreMatrixRow.hh"
#include "mcf_big_seq.hh"

#include <stdexcept>

namespace cbrc {

using namespace mcf;

static int forwardGaplessPssmXdropScore(BigPtr seq,
					const ScoreMatrixRow *pssm,
					int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*pssm++)[getNext(seq)];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in forward gapless extension with PSSM");
  return score;
}

static int reverseGaplessPssmXdropScore(BigPtr seq,
					const ScoreMatrixRow *pssm,
					int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*--pssm)[getPrev(seq)];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in reverse gapless extension with PSSM");
  return score;
}

static void gaplessPssmXdropScores(BigSeq seq, const ScoreMatrixRow *pssm,
				   int maxScoreDrop,
				   size_t pos1, size_t pos2,
				   int &fwdScore, int &revScore) {
  fwdScore = forwardGaplessPssmXdropScore(seq+pos1, pssm+pos2, maxScoreDrop);
  revScore = reverseGaplessPssmXdropScore(seq+pos1, pssm+pos2, maxScoreDrop);
}

static bool gaplessPssmXdropEnds(BigSeq seq, const ScoreMatrixRow *pssm,
				 int maxScoreDrop, int fwdScore, int revScore,
				 size_t &pos1, size_t &pos2, size_t &length) {
  size_t beg1 = pos1;
  size_t end1 = beg1;
  size_t beg2 = pos2;
  size_t end2 = beg2;
  while (fwdScore) fwdScore -= pssm[end2++][seq[end1++]];
  while (revScore) revScore -= pssm[--beg2][seq[--beg1]];
  pos1 = beg1;
  pos2 = beg2;
  length = end1 - beg1;

  int score = 0;
  int maxScore = 0;
  while (beg1 < end1) {
    score += pssm[beg2++][seq[beg1++]];
    if (score > maxScore) {
      maxScore = score;
    } else if (score <= 0 || beg1 == end1 || score < maxScore - maxScoreDrop) {
      return false;
    }
  }
  return true;
}

static int gaplessPssmXdropOverlap(BigPtr seq,
				   const ScoreMatrixRow *pssm,
				   int maxScoreDrop,
				   size_t &reverseLength,
				   size_t &forwardLength) {
  int minScore = 0;
  int maxScore = 0;
  int score = 0;

  BigPtr rs = seq;
  const ScoreMatrixRow *rp = pssm;
  while (true) {
    int s = (*--rp)[getPrev(rs)];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    else if (score < minScore) minScore = score;
  }

  maxScore = score - minScore;

  const ScoreMatrixRow *fp = pssm;
  while (true) {
    int s = (*fp++)[getNext(seq)];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
  }

  reverseLength = pssm - (rp + 1);
  forwardLength = (fp - 1) - pssm;
  return score;
}

static int gaplessPssmAlignmentScore(BigPtr seq,
				     const ScoreMatrixRow *pssm,
				     size_t length) {
  int score = 0;
  while (length--) score += (*pssm++)[getNext(seq)];
  return score;
}

}

#endif
