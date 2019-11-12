// Copyright 2011, 2012, 2013 Martin C. Frith

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

namespace cbrc {

int GappedXdropAligner::alignPssm(const uchar *seq,
                                  const ScoreMatrixRow *pssm,
                                  bool isForward,
				  int globality,
				  int delExistenceCost,
				  int delExtensionCost,
				  int insExistenceCost,
				  int insExtensionCost,
                                  int gapUnalignedCost,
                                  int maxScoreDrop,
                                  int maxMatchScore) {
  const SimdInt mNegInf = simdSet1(-INF);
  const int seqIncrement = isForward ? 1 : -1;
  const bool isAffine = isAffineGaps(delExistenceCost, delExtensionCost,
				     insExistenceCost, insExtensionCost,
				     gapUnalignedCost);

  size_t maxSeq1begs[] = { 0, 9 };
  size_t minSeq1ends[] = { 1, 0 };

  int bestScore = 0;
  int bestEdgeScore = -INF;
  size_t bestEdgeAntidiagonal = 0;

  init();

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    size_t seq1beg = std::min(maxSeq1begs[0], maxSeq1begs[1]);
    size_t seq1end = std::max(minSeq1ends[0], minSeq1ends[1]);

    if (seq1beg >= seq1end) break;

    size_t scoreEnd = scoreEnds.back();
    int numCells = seq1end - seq1beg;

    initAntidiagonal(seq1end, scoreEnd + xdropPadLen + numCells);

    size_t seq2pos = antidiagonal - seq1beg;

    const uchar *s1 = isForward ? seq + seq1beg : seq - seq1beg;
    const ScoreMatrixRow *s2 = isForward ? pssm + seq2pos : pssm - seq2pos;

    if (!globality && isDelimiter(0, *s2))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    int minScore = bestScore - maxScoreDrop;

    int *x0 = &xScores[scoreEnd];
    int *y0 = &yScores[scoreEnd];
    int *z0 = &zScores[scoreEnd];
    const int *y1 = &yScores[hori(antidiagonal, seq1beg)];
    const int *z1 = &zScores[vert(antidiagonal, seq1beg)];
    const int *x2 = &xScores[diag(antidiagonal, seq1beg)];

    simdStore(x0, mNegInf);  x0 += xdropPadLen;
    simdStore(y0, mNegInf);  y0 += xdropPadLen;
    simdStore(z0, mNegInf);  z0 += xdropPadLen;

    if (globality && isDelimiter(0, *s2)) {
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      int b = maxValue(x2[0], z1[0]-insExtensionCost, z2[0]-gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1beg);
    }

    if (isAffine) {
      if (isForward)
	for (int i = 0; i < numCells; ++i) {
          int x = x2[i];
          int y = y1[i] - delExtensionCost;
          int z = z1[i] - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
	    if (b > bestScore) {
	      bestScore = b;
	      bestAntidiagonal = antidiagonal;
	    }
            x0[i] = b + (*s2)[*s1];
            y0[i] = maxValue(b - delExistenceCost, y);
            z0[i] = maxValue(b - insExistenceCost, z);
          }
          else x0[i] = y0[i] = z0[i] = -INF;
          ++s1;
	  --s2;
        }
      else
	for (int i = 0; i < numCells; ++i) {
          int x = x2[i];
          int y = y1[i] - delExtensionCost;
          int z = z1[i] - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
	    if (b > bestScore) {
	      bestScore = b;
	      bestAntidiagonal = antidiagonal;
	    }
            x0[i] = b + (*s2)[*s1];
            y0[i] = maxValue(b - delExistenceCost, y);
            z0[i] = maxValue(b - insExistenceCost, z);
          }
          else x0[i] = y0[i] = z0[i] = -INF;
          --s1;
	  ++s2;
        }
    } else {
      const int *y2 = &yScores[diag(antidiagonal, seq1beg)];
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      for (int i = 0; i < numCells; ++i) {
        int x = x2[i];
        int y = maxValue(y1[i] - delExtensionCost, y2[i] - gapUnalignedCost);
        int z = maxValue(z1[i] - insExtensionCost, z2[i] - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
	  if (b > bestScore) {
	    bestScore = b;
	    bestAntidiagonal = antidiagonal;
	  }
          x0[i] = b + (*s2)[*s1];
          y0[i] = maxValue(b - delExistenceCost, y);
          z0[i] = maxValue(b - insExistenceCost, z);
        }
        else x0[i] = y0[i] = z0[i] = -INF;
	s1 += seqIncrement;
	s2 -= seqIncrement;
      }
    }

    s1 -= seqIncrement;

    if (globality && isDelimiter(*s1, *pssm)) {
      const int *y2 = &yScores[diag(antidiagonal, seq1beg)];
      int n = numCells - 1;
      int b = maxValue(x2[n], y1[n]-delExtensionCost, y2[n]-gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1end-1);
    }

    if (!globality && isDelimiter(*s1, *pssm))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    const int *x0base = x0 - seq1beg;
    updateFiniteEdges(maxSeq1begs, minSeq1ends, x0base, x0+numCells, numCells);
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
