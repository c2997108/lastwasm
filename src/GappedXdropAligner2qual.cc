// Copyright 2011, 2012, 2013 Martin C. Frith

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"
#include "TwoQualityScoreMatrix.hh"

namespace cbrc {

static bool isDelimiter2qual(uchar c) {
  return c == 4;  // a bit ugly (hard-wired delimiter value)
}

int GappedXdropAligner::align2qual(const uchar *seq1,
                                   const uchar *qual1,
                                   const uchar *seq2,
                                   const uchar *qual2,
                                   bool isForward,
				   int globality,
                                   const TwoQualityScoreMatrix &scorer,
				   int delExistenceCost,
				   int delExtensionCost,
				   int insExistenceCost,
				   int insExtensionCost,
                                   int gapUnalignedCost,
				   bool isAffine,
                                   int maxScoreDrop,
                                   int maxMatchScore) {
  const int seqIncrement = isForward ? 1 : -1;

  size_t seq1beg = 0;
  size_t seq1end = 1;
  size_t diagPos = xdropPadLen - 1;
  size_t horiPos = xdropPadLen * 2 - 1;
  size_t thisPos = xdropPadLen * 2;

  int bestScore = 0;
  int bestEdgeScore = -INF;
  size_t bestEdgeAntidiagonal = 0;

  init();

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    int numCells = seq1end - seq1beg;
    int n = numCells - 1;

    size_t seq2pos = antidiagonal - seq1beg;

    const uchar *s1 = isForward ? seq1 + seq1beg : seq1 - seq1beg;
    const uchar *q1 = isForward ? qual1 + seq1beg : qual1 - seq1beg;
    const uchar *s2 = isForward ? seq2 + seq2pos : seq2 - seq2pos;
    const uchar *q2 = isForward ? qual2 + seq2pos : qual2 - seq2pos;

    initAntidiagonal(seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    Score *x0 = &xScores[thisPos];
    Score *y0 = &yScores[thisPos];
    Score *z0 = &zScores[thisPos];
    const Score *y1 = &yScores[horiPos];
    const Score *z1 = &zScores[horiPos + 1];
    const Score *x2 = &xScores[diagPos];

    if (!globality && (isDelimiter2qual(s1[n * seqIncrement]) ||
		       isDelimiter2qual(s2[0]))) {
      updateMaxScoreDrop(maxScoreDrop, n, maxMatchScore);
    }

    int minScore = bestScore - maxScoreDrop;

    if (globality && isDelimiter2qual(s2[0])) {
      const Score *z2 = &zScores[diagPos];
      int b = maxValue(x2[0], z1[0]-insExtensionCost, z2[0]-gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1beg);
    }

    if (globality && isDelimiter2qual(s1[n * seqIncrement])) {
      const Score *y2 = &yScores[diagPos];
      int b = maxValue(x2[n], y1[n]-delExtensionCost, y2[n]-gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1end-1);
    }

    const Score *x0last = x0 + n;

    if (isAffine) {
      if (isForward)
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
	    if (b > bestScore) {
	      bestScore = b;
	      bestAntidiagonal = antidiagonal;
	    }
            *x0 = b + scorer(*s1, *s2, *q1, *q2);
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  ++q1;  --s2;  --q2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
      else
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
	    if (b > bestScore) {
	      bestScore = b;
	      bestAntidiagonal = antidiagonal;
	    }
            *x0 = b + scorer(*s1, *s2, *q1, *q2);
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  --q1;  ++s2;  ++q2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
    } else {
      const Score *y2 = &yScores[diagPos];
      const Score *z2 = &zScores[diagPos];
      while (1) {
        int x = *x2;
        int y = maxValue(*y1 - delExtensionCost, *y2 - gapUnalignedCost);
        int z = maxValue(*z1 - insExtensionCost, *z2 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
	  if (b > bestScore) {
	    bestScore = b;
	    bestAntidiagonal = antidiagonal;
	  }
          *x0 = b + scorer(*s1, *s2, *q1, *q2);
          *y0 = maxValue(b - delExistenceCost, y);
          *z0 = maxValue(b - insExistenceCost, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;  ++y2;  ++z2;
	s1 += seqIncrement;  q1 += seqIncrement;
	s2 -= seqIncrement;  q2 -= seqIncrement;
      }
    }

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[0] > -INF / 2) {
      ++seq1end;
    }

    if (x0[-n] <= -INF / 2) {
      ++seq1beg;
      ++diagPos;
      ++horiPos;
      if (seq1beg >= seq1end) break;
    }
  }

  if (globality) {
    bestAntidiagonal = bestEdgeAntidiagonal;
    bestScore = bestEdgeScore;
  } else {
    calcBestSeq1position(bestScore);
  }
  return bestScore;
}

}
