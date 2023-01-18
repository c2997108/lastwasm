// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

// Functions that find gapless X-drop alignments between two sequences.

#ifndef GAPLESS_XDROP_HH
#define GAPLESS_XDROP_HH

#include "ScoreMatrixRow.hh"

#include <stddef.h>
#include <stdexcept>

namespace cbrc {

typedef unsigned char uchar;

// Gets the maximum score for any gapless alignment starting at (seq1,
// seq2) and extending forwards.  The score is not allowed to drop by
// more than maxScoreDrop.  The sequences had better end with
// sentinels that have score < -maxScoreDrop.
// The score might suffer overflow, for huge sequences and/or huge
// scores.  If the function detects this (not guaranteed), it throws
// an exception.
static int forwardGaplessXdropScore(const uchar *seq1,
				    const uchar *seq2,
				    const ScoreMatrixRow *scorer,
				    int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += scorer[*seq1++][*seq2++];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in forward gapless extension");
  return score;
}

// As above, but extending backwards.
static int reverseGaplessXdropScore(const uchar *seq1,
				    const uchar *seq2,
				    const ScoreMatrixRow *scorer,
				    int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += scorer[*--seq1][*--seq2];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    throw std::overflow_error("score overflow in reverse gapless extension");
  return score;
}

// Find the shortest forward extension from (pos1, pos2) with score
// "fwdScore", and the shortest reverse extension with score
// "revScore".  Return the start coordinates and length of this alignment.
static bool gaplessXdropEnds(const uchar *seq1, const uchar *seq2,
			     const ScoreMatrixRow *scorer, int maxScoreDrop,
			     int fwdScore, int revScore,
			     size_t &pos1, size_t &pos2, size_t &length) {
  size_t beg1 = pos1;
  size_t end1 = beg1;
  size_t beg2 = pos2;
  size_t end2 = beg2;
  while (fwdScore) fwdScore -= scorer[seq1[end1++]][seq2[end2++]];
  while (revScore) revScore -= scorer[seq1[--beg1]][seq2[--beg2]];
  pos1 = beg1;
  pos2 = beg2;
  length = end1 - beg1;

  // Check whether the alignment has no prefix with score <= 0, no
  // suffix with score <= 0, and no region with score < -maxScoreDrop
  int score = 0;
  int maxScore = 0;
  while (beg1 < end1) {
    score += scorer[seq1[beg1++]][seq2[beg2++]];
    if (score > maxScore) {
      maxScore = score;
    } else if (score <= 0 || beg1 == end1 || score < maxScore - maxScoreDrop) {
      return false;
    }
  }
  return true;
}

// Returns the score, and sets the reverse and forward extension
// lengths, for a gapless "overlap" alignment starting at (seq1,
// seq2).  "Overlap" means that the alignment must extend, in each
// direction, until it hits a score <= -INF (presumably from a
// sentinel indicating a sequence end).  If the alignment would have
// any region with score < -maxScoreDrop, -INF is returned and the
// extension lengths are not set.
static int gaplessXdropOverlap(const uchar *seq1,
			       const uchar *seq2,
			       const ScoreMatrixRow *scorer,
			       int maxScoreDrop,
			       size_t &reverseLength,
			       size_t &forwardLength) {
  int minScore = 0;
  int maxScore = 0;
  int score = 0;

  const uchar *r1 = seq1;
  const uchar *r2 = seq2;
  while (true) {
    --r1;  --r2;
    int s = scorer[*r1][*r2];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    else if (score < minScore) minScore = score;
  }

  maxScore = score - minScore;

  const uchar *f1 = seq1;
  const uchar *f2 = seq2;
  while (true) {
    int s = scorer[*f1][*f2];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    ++f1;  ++f2;
  }

  reverseLength = seq1 - (r1 + 1);
  forwardLength = f1 - seq1;
  return score;
}

// Calculate the score of the gapless alignment starting at (seq1,
// seq2) and ending at seq1end.
static int gaplessAlignmentScore(const uchar *seq1,
				 const uchar *seq1end,
				 const uchar *seq2,
				 const ScoreMatrixRow *scorer) {
  int score = 0;
  while (seq1 < seq1end) score += scorer[*seq1++][*seq2++];
  return score;
}

}

#endif
