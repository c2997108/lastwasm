// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

//#include <iostream>  // for debugging

namespace cbrc {

const int seqLoadLen = 16;

const int delimiter = 4;

int GappedXdropAligner::alignDna(const uchar *seq1,
				 const uchar *seq2,
				 bool isForward,
				 const ScoreMatrixRow *scorer,
				 int delOpenCost,
				 int delGrowCost,
				 int insOpenCost,
				 int insGrowCost,
				 int maxScoreDrop,
				 int maxMatchScore,
				 const uchar *toUnmasked) {
  // This guarantees that we avoid false-positive drops > maxScoreDrop:
  assert(maxScoreDrop < SCHAR_MIN - droppedTinyScore);
  // With non-saturating arithmetic, we also require:
  // -maxScoreDrop + shortDelimiterScore >= SHRT_MIN,
  // which is guaranteed by the definition of shortDelimiterScore.

  delGrowCost = std::min(delGrowCost, maxScoreDrop + 1);
  delOpenCost = std::min(delOpenCost, maxScoreDrop + 1 - delGrowCost);

  insGrowCost = std::min(insGrowCost, maxScoreDrop + 1);
  insOpenCost = std::min(insOpenCost, maxScoreDrop + 1 - insGrowCost);

  const SimdInt mNegInf = simdFill2(droppedTinyScore);
  const SimdInt mDelOpenCost = simdFill2(delOpenCost);
  const SimdInt mDelGrowCost = simdFill2(delGrowCost);
  const SimdInt mInsOpenCost = simdFill2(insOpenCost);
  const SimdInt mInsGrowCost = simdFill2(insGrowCost);
  const int seqIncrement = isForward ? 1 : -1;

  const __m128i reverse = _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);

  const __m128i scorer4x4 =
    _mm_set_epi8(scorer[3][3], scorer[3][2], scorer[3][1], scorer[3][0],
		 scorer[2][3], scorer[2][2], scorer[2][1], scorer[2][0],
		 scorer[1][3], scorer[1][2], scorer[1][1], scorer[1][0],
		 scorer[0][3], scorer[0][2], scorer[0][1], scorer[0][0]);

  size_t seq1beg = 0;
  size_t seq1end = 1;
  size_t diagPos = xdropPadLen - 1;
  size_t horiPos = xdropPadLen * 2 - 1;
  size_t thisPos = xdropPadLen * 2;

  int bestScore = 0;
  SimdInt mBestScore = simdZero();
  SimdInt mMinScore = simdFill2(-maxScoreDrop);
  SimdInt mScoreRise1 = simdZero();
  SimdInt mScoreRise2 = simdZero();

  initTiny();
  seq1queue.clear();
  seq2queue.clear();

  bool isDna = (toUnmasked[*seq1] < 4 && toUnmasked[*seq2] < 4);

  for (int i = 0; i < seqLoadLen; ++i) {
    uchar x = toUnmasked[*seq1];
    seq1queue.push(x, i);
    seq1 += seqIncrement * (x != delimiter);
    seq2queue.push(toUnmasked[*seq2], i);
  }

  seq2 += seqIncrement;

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    int numCells = seq1end - seq1beg;
    int n = numCells - 1;

    const uchar *s1 = &seq1queue.fromEnd(n + seqLoadLen);
    const uchar *s2 = &seq2queue.fromEnd(1);

    initAntidiagonalTiny(seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    TinyScore *x0 = &xTinyScores[thisPos];
    TinyScore *y0 = &yTinyScores[thisPos];
    TinyScore *z0 = &zTinyScores[thisPos];
    const TinyScore *y1 = &yTinyScores[horiPos];
    const TinyScore *z1 = &zTinyScores[horiPos + 1];
    const TinyScore *x2 = &xTinyScores[diagPos];

    const SimdInt mScoreRise12 = simdAdd2(mScoreRise1, mScoreRise2);
    const SimdInt mDelGrowCost1 = simdAdd2(mDelGrowCost, mScoreRise1);
    const SimdInt mInsGrowCost1 = simdAdd2(mInsGrowCost, mScoreRise1);

    if (isDna) {
      for (int i = 0; i < numCells; i += simdLen2) {
	__m128i fwd1 = simdLoad128(s1);
	__m128i fwd2 = simdLoad128(s2 - (seqLoadLen-1));
	__m128i rev2 = _mm_shuffle_epi8(fwd2, reverse);
	__m128i j = _mm_or_si128(_mm_slli_epi32(fwd1, 2), rev2);
	SimdInt s = simdBytesToInts2(_mm_shuffle_epi8(scorer4x4, j));
	SimdInt x = simdSub2(simdLoad(x2+i), mScoreRise12);
	SimdInt y = simdSub2(simdLoad(y1+i), mDelGrowCost1);
	SimdInt z = simdSub2(simdLoad(z1+i), mInsGrowCost1);
	SimdInt b = simdMax2(simdMax2(x, y), z);
	SimdInt isDrop = simdGt2(mMinScore, b);
	mBestScore = simdMax2(b, mBestScore);
	simdStore(x0+i, simdBlend(simdAdd2(b, s), mNegInf, isDrop));
	simdStore(y0+i, simdMax2(simdSub2(b, mDelOpenCost), y));
	simdStore(z0+i, simdMax2(simdSub2(b, mInsOpenCost), z));
	s1 += simdLen2;
	s2 -= simdLen2;
      }
    } else {
      if (s1[n] == delimiter || s2[0] == delimiter) {
	updateMaxScoreDrop(maxScoreDrop, n, maxMatchScore);
	mMinScore = simdFill2(-maxScoreDrop);
      }

      for (int i = 0; i < numCells; i += simdLen2) {
	SimdInt s = simdSet2(
#ifdef __SSE4_1__
#ifdef __AVX2__
			     scorer[s1[15]][s2[-15]],
			     scorer[s1[14]][s2[-14]],
			     scorer[s1[13]][s2[-13]],
			     scorer[s1[12]][s2[-12]],
			     scorer[s1[11]][s2[-11]],
			     scorer[s1[10]][s2[-10]],
			     scorer[s1[9]][s2[-9]],
			     scorer[s1[8]][s2[-8]],
#endif
			     scorer[s1[7]][s2[-7]],
			     scorer[s1[6]][s2[-6]],
			     scorer[s1[5]][s2[-5]],
			     scorer[s1[4]][s2[-4]],
			     scorer[s1[3]][s2[-3]],
			     scorer[s1[2]][s2[-2]],
			     scorer[s1[1]][s2[-1]],
#endif
			     scorer[s1[0]][s2[-0]]);

	SimdInt x = simdSub2(simdLoad(x2+i), mScoreRise12);
	SimdInt y = simdSub2(simdLoad(y1+i), mDelGrowCost1);
	SimdInt z = simdSub2(simdLoad(z1+i), mInsGrowCost1);
	SimdInt b = simdMax2(simdMax2(x, y), z);
	SimdInt isDrop = simdGt2(mMinScore, b);
	mBestScore = simdMax2(b, mBestScore);
	simdStore(x0+i, simdBlend(simdAdd2(b, s), mNegInf, isDrop));
	simdStore(y0+i, simdMax2(simdSub2(b, mDelOpenCost), y));
	simdStore(z0+i, simdMax2(simdSub2(b, mInsOpenCost), z));
	s1 += simdLen2;
	s2 -= simdLen2;
      }
    }

    mScoreRise2 = mScoreRise1;
    mScoreRise1 = simdZero();
    int newBestScore = simdHorizontalMax2(mBestScore);
    if (newBestScore > 0) {
      bestScore += newBestScore;
      bestAntidiagonal = antidiagonal;
      mBestScore = simdZero();
      mScoreRise1 = simdFill2(newBestScore);
    }
    scoreRises.push_back(newBestScore);

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[n] > droppedTinyScore) {
      ++seq1end;
      uchar x = toUnmasked[*seq1];
      seq1queue.push(x, n + seqLoadLen);
      seq1 += seqIncrement * (x != delimiter);
      uchar z = seq1queue.fromEnd(seqLoadLen);
      if (z >= 4) {
	isDna = false;
      }
    }

    if (x0[0] > droppedTinyScore) {
      uchar y = toUnmasked[*seq2];
      seq2queue.push(y, n + seqLoadLen);
      seq2 += seqIncrement;
      if (y >= 4) {
	isDna = false;
      }
    } else {
      ++seq1beg;
      ++diagPos;
      ++horiPos;
      if (seq1beg == seq1end) break;
    }
  }

  calcBestSeq1positionTiny();
  return bestScore;
}

bool GappedXdropAligner::getNextChunkDna(size_t &end1,
					 size_t &end2,
					 size_t &length,
					 int delOpenCost,
					 int delGrowCost,
					 int insOpenCost,
					 int insGrowCost) {
  if (bestAntidiagonal == 0) return false;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - bestSeq1position;

  int x, y, z;
  while (1) {
    size_t h = hori(bestAntidiagonal, bestSeq1position);
    size_t v = vert(bestAntidiagonal, bestSeq1position);
    size_t d = diag(bestAntidiagonal, bestSeq1position);
    x = xTinyScores[d] - scoreRises[bestAntidiagonal];
    y = yTinyScores[h] - delGrowCost;
    z = zTinyScores[v] - insGrowCost;
    if (x < y || x < z || bestAntidiagonal == 0) break;
    bestAntidiagonal -= 2;
    bestSeq1position -= 1;
  }

  length = end1 - bestSeq1position;
  if (bestAntidiagonal == 0) return true;

  while (1) {
    bool isDel = (y >= z);
    bestAntidiagonal -= 1;
    if (isDel) bestSeq1position -= 1;
    size_t h = hori(bestAntidiagonal, bestSeq1position);
    size_t v = vert(bestAntidiagonal, bestSeq1position);
    size_t d = diag(bestAntidiagonal, bestSeq1position);
    x = xTinyScores[d] - scoreRises[bestAntidiagonal];
    y = yTinyScores[h] - delGrowCost;
    z = zTinyScores[v] - insGrowCost;
    if (isDel) {
      y += delOpenCost;
    } else {
      z += insOpenCost;
    }
    if (x >= y && x >= z) return true;
  }
}

}
