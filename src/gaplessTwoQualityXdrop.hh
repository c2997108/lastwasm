// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

// Functions that find gapless X-drop alignments between two sequences
// that both have quality scores.

// These functions are analogous to those described in
// gaplessXdrop.hh.

#ifndef GAPLESS_TWO_QUALITY_XDROP_HH
#define GAPLESS_TWO_QUALITY_XDROP_HH

#include "TwoQualityScoreMatrix.hh"

#include <stddef.h>
#include <stdexcept>

namespace cbrc {

typedef unsigned char uchar;

static int forwardGaplessTwoQualityXdropScore(const uchar *seq1,
					      const uchar *qual1,
					      const uchar *seq2,
					      const uchar *qual2,
					      const TwoQualityScoreMatrix &m,
					      int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += m(*seq1++, *seq2++, *qual1++, *qual2++);  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in forward gapless extension with qualities");
  return score;
}

static int reverseGaplessTwoQualityXdropScore(const uchar *seq1,
					      const uchar *qual1,
					      const uchar *seq2,
					      const uchar *qual2,
					      const TwoQualityScoreMatrix &m,
					      int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += m(*--seq1, *--seq2, *--qual1, *--qual2);  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in reverse gapless extension with qualities");
  return score;
}

static bool gaplessTwoQualityXdropEnds(const uchar *seq1, const uchar *qual1,
				       const uchar *seq2, const uchar *qual2,
				       const TwoQualityScoreMatrix &m,
				       int maxScoreDrop,
				       int fwdScore, int revScore,
				       size_t &pos1, size_t &pos2,
				       size_t &length) {
  size_t beg1 = pos1;
  size_t end1 = beg1;
  size_t beg2 = pos2;
  size_t end2 = beg2;
  while (fwdScore) {
    fwdScore -= m(seq1[end1], seq2[end2], qual1[end1], qual2[end2]);
    end1++;
    end2++;
  }
  while (revScore) {
    --beg1;
    --beg2;
    revScore -= m(seq1[beg1], seq2[beg2], qual1[beg1], qual2[beg2]);
  }
  pos1 = beg1;
  pos2 = beg2;
  length = end1 - beg1;

  int score = 0;
  int maxScore = 0;
  while (beg1 < end1) {
    score += m(seq1[beg1], seq2[beg2], qual1[beg1], qual2[beg2]);
    beg1++;
    beg2++;
    if (score > maxScore) {
      maxScore = score;
    } else if (score <= 0 || beg1 == end1 || score < maxScore - maxScoreDrop) {
      return false;
    }
  }
  return true;
}

static int gaplessTwoQualityXdropOverlap(const uchar *seq1,
					 const uchar *qual1,
					 const uchar *seq2,
					 const uchar *qual2,
					 const TwoQualityScoreMatrix &m,
					 int maxScoreDrop,
					 size_t &reverseLength,
					 size_t &forwardLength) {
  int minScore = 0;
  int maxScore = 0;
  int score = 0;

  const uchar *rs1 = seq1;
  const uchar *rq1 = qual1;
  const uchar *rs2 = seq2;
  const uchar *rq2 = qual2;
  while (true) {
    --rs1;  --rq1;  --rs2;  --rq2;
    int s = m(*rs1, *rs2, *rq1, *rq2);
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    else if (score < minScore) minScore = score;
  }

  maxScore = score - minScore;

  const uchar *fs1 = seq1;
  const uchar *fq1 = qual1;
  const uchar *fs2 = seq2;
  const uchar *fq2 = qual2;
  while (true) {
    int s = m(*fs1, *fs2, *fq1, *fq2);
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    ++fs1;  ++fq1;  ++fs2;  ++fq2;
  }

  reverseLength = seq1 - (rs1 + 1);
  forwardLength = fs1 - seq1;
  return score;
}

static int gaplessTwoQualityAlignmentScore(const uchar *seq1,
					   const uchar *seq1end,
					   const uchar *qual1,
					   const uchar *seq2,
					   const uchar *qual2,
					   const TwoQualityScoreMatrix &m) {
  int score = 0;
  while (seq1 < seq1end) score += m(*seq1++, *seq2++, *qual1++, *qual2++);
  return score;
}

}

#endif
