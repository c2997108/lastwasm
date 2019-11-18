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
  const int *vectorOfMatchScores = *pssm;
  const SimdInt mNegInf = simdSet1(-INF);
  const SimdInt mDelOpenCost = simdSet1(delExistenceCost);
  const SimdInt mDelGrowCost = simdSet1(delExtensionCost);
  const SimdInt mInsOpenCost = simdSet1(insExistenceCost);
  const SimdInt mInsGrowCost = simdSet1(insExtensionCost);
  const int seqIncrement = isForward ? 1 : -1;
  const bool isAffine = isAffineGaps(delExistenceCost, delExtensionCost,
				     insExistenceCost, insExtensionCost,
				     gapUnalignedCost);

  size_t maxSeq1begs[] = { 0, 9 };
  size_t minSeq1ends[] = { 1, 0 };

  int bestScore = 0;
  SimdInt mBestScore = simdSet1(0);
  int bestEdgeScore = -INF;
  size_t bestEdgeAntidiagonal = 0;

  init();
  seq1queue.clear();
  seq2queue.clear();

  // xxx the queue names are flipped: seq1queue <=> seq2queue

  for (int i = 0; i < simdLen-1; ++i) {
    uchar c = *seq;
    seq2queue.push(c, 0);
    seq += seqIncrement * !isDelimiter(c, vectorOfMatchScores);
    seq1queue.push(vectorOfMatchScores, 0);
  }

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    size_t seq1beg = std::min(maxSeq1begs[0], maxSeq1begs[1]);
    size_t seq1end = std::max(minSeq1ends[0], minSeq1ends[1]);

    if (seq1beg >= seq1end) break;

    size_t scoreEnd = scoreEnds.back();
    int numCells = seq1end - seq1beg;

    initAntidiagonal(seq1end, scoreEnd + xdropPadLen + numCells);

    if (seq1end + (simdLen-1) > seq2queue.size()) {
      uchar c = *seq;
      seq2queue.push(c, seq1beg);
      seq += seqIncrement * !isDelimiter(c, vectorOfMatchScores);
    }
    const uchar *s1 = &seq2queue[seq1beg];

    size_t seq2pos = antidiagonal - seq1beg;
    if (seq2pos + simdLen > seq1queue.size()) {
      seq1queue.push(*pssm, seq2pos - numCells + 1);
      pssm += seqIncrement;
    }
    const const_int_ptr *s2 = &seq1queue[seq2pos + (simdLen-1)];

    if (!globality && isDelimiter(0, *s2))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    int minScore = bestScore - maxScoreDrop;
    SimdInt mMinScore = simdSet1(minScore);

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
      for (int i = 0; i < numCells; i += simdLen) {
	SimdInt s = simdSet(
#ifdef __SSE4_1__
//#ifdef __AVX2__
#ifdef WANT_AVX2
			    s2[-7][s1[7]],
			    s2[-6][s1[6]],
			    s2[-5][s1[5]],
			    s2[-4][s1[4]],
#endif
			    s2[-3][s1[3]],
			    s2[-2][s1[2]],
			    s2[-1][s1[1]],
#endif
			    s2[-0][s1[0]]);
	SimdInt x = simdLoad(x2+i);
	SimdInt y = simdSub(simdLoad(y1+i), mDelGrowCost);
	SimdInt z = simdSub(simdLoad(z1+i), mInsGrowCost);
	SimdInt b = simdMax(simdMax(x, y), z);
	SimdInt isDrop = simdGt(mMinScore, b);
	mBestScore = simdMax(b, mBestScore);
	simdStore(x0+i, simdBlend(simdAdd(b, s), mNegInf, isDrop));
	simdStore(y0+i, simdMax(simdSub(b, mDelOpenCost), y));
	simdStore(z0+i, simdMax(simdSub(b, mInsOpenCost), z));
	s1 += simdLen;
	s2 -= simdLen;
      }

      int newBestScore = simdHorizontalMax(mBestScore);
      if (newBestScore > bestScore) {
	bestScore = newBestScore;
	bestAntidiagonal = antidiagonal;
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
	++s1;
	--s2;
      }
    }

    uchar seq1back = seq2queue[seq1end - 1];

    if (globality && isDelimiter(seq1back, vectorOfMatchScores)) {
      const int *y2 = &yScores[diag(antidiagonal, seq1beg)];
      int n = numCells - 1;
      int b = maxValue(x2[n], y1[n]-delExtensionCost, y2[n]-gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		    b, antidiagonal, seq1end-1);
    }

    if (!globality && isDelimiter(seq1back, vectorOfMatchScores))
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
