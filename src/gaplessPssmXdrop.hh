// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

// Functions that find gapless X-drop alignments between a sequence
// and a PSSM.

// These functions are analogous to those described in
// gaplessXdrop.hh.  Here, "pssm" replaces both "seq2" and "scorer".

#ifndef GAPLESS_PSSM_XDROP_HH
#define GAPLESS_PSSM_XDROP_HH

#include "ScoreMatrixRow.hh"

#include <stddef.h>
#include <stdexcept>

namespace cbrc {

typedef unsigned char uchar;

static int forwardGaplessPssmXdropScore(const uchar *seq,
					const ScoreMatrixRow *pssm,
					int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*pssm++)[*seq++];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in forward gapless extension with PSSM");
  return score;
}

static int reverseGaplessPssmXdropScore(const uchar *seq,
					const ScoreMatrixRow *pssm,
					int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*--pssm)[*--seq];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in reverse gapless extension with PSSM");
  return score;
}

static bool gaplessPssmXdropEnds(const uchar *seq, const ScoreMatrixRow *pssm,
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

static int gaplessPssmXdropOverlap(const uchar *seq,
				   const ScoreMatrixRow *pssm,
				   int maxScoreDrop,
				   size_t &reverseLength,
				   size_t &forwardLength) {
  int minScore = 0;
  int maxScore = 0;
  int score = 0;

  const uchar *rs = seq;
  const ScoreMatrixRow *rp = pssm;
  while (true) {
    --rs;  --rp;
    int s = (*rp)[*rs];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    else if (score < minScore) minScore = score;
  }

  maxScore = score - minScore;

  const uchar *fs = seq;
  const ScoreMatrixRow *fp = pssm;
  while (true) {
    int s = (*fp)[*fs];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    ++fs;  ++fp;
  }

  reverseLength = seq - (rs + 1);
  forwardLength = fs - seq;
  return score;
}

static int gaplessPssmAlignmentScore(const uchar *seq,
				     const uchar *seqEnd,
				     const ScoreMatrixRow *pssm) {
  int score = 0;
  while (seq < seqEnd) score += (*pssm++)[*seq++];
  return score;
}

}

#endif
