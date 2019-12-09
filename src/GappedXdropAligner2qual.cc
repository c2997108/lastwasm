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
                                   int maxScoreDrop,
                                   int maxMatchScore) {
  const SimdInt mNegInf = simdSet1(-INF);
  const bool isAffine = isAffineGaps(delExistenceCost, delExtensionCost,
				     insExistenceCost, insExtensionCost,
				     gapUnalignedCost);

  size_t seq1beg = 0;
  size_t seq1end = 1;

  int bestScore = 0;
  int bestEdgeScore = -INF;
  size_t bestEdgeAntidiagonal = 0;

  init();

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    size_t scoreEnd = scoreEnds.back();
    size_t numCells = seq1end - seq1beg;

    size_t seq2pos = antidiagonal - seq1beg;

    const uchar *s1 = isForward ? seq1 + seq1beg : seq1 - seq1beg;
    const uchar *q1 = isForward ? qual1 + seq1beg : qual1 - seq1beg;
    const uchar *s2 = isForward ? seq2 + seq2pos : seq2 - seq2pos;
    const uchar *q2 = isForward ? qual2 + seq2pos : qual2 - seq2pos;

    initAntidiagonal(seq1end, scoreEnd + xdropPadLen + numCells);

    Score *x0 = &xScores[scoreEnd];
    Score *y0 = &yScores[scoreEnd];
    Score *z0 = &zScores[scoreEnd];
    const Score *y1 = &yScores[hori(antidiagonal, seq1beg)];
    const Score *z1 = &zScores[vert(antidiagonal, seq1beg)];
    const Score *x2 = &xScores[diag(antidiagonal, seq1beg)];

    simdStore(x0, mNegInf);  x0 += xdropPadLen;
    simdStore(y0, mNegInf);  y0 += xdropPadLen;
    simdStore(z0, mNegInf);  z0 += xdropPadLen;

    if (!globality && isDelimiter2qual(*s2))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    int minScore = bestScore - maxScoreDrop;

    const Score *x0last = x0 + numCells - 1;

    if (globality && isDelimiter2qual(*s2)) {
      const Score *z2 = &zScores[diag(antidiagonal, seq1beg)];
      int b = maxValue(*x2, *z1 - insExtensionCost, *z2 - gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1beg);
    }

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
      const Score *y2 = &yScores[diag(antidiagonal, seq1beg)];
      const Score *z2 = &zScores[diag(antidiagonal, seq1beg)];
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
        if (isForward) { ++s1;  ++q1;  --s2;  --q2; }
        else           { --s1;  --q1;  ++s2;  ++q2; }
      }
    }

    if (globality && isDelimiter2qual(*s1)) {
      const Score *y2 = &yScores[diag(antidiagonal, seq1end-1)];
      int b = maxValue(*x2, *y1 - delExtensionCost, *y2 - gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1end-1);
    }

    if (!globality && isDelimiter2qual(*s1))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    if (x0[0] > -INF / 2) {
      ++seq1end;
    }

    if (x0[1 - numCells] <= -INF / 2) {
      ++seq1beg;
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
