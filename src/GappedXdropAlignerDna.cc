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
				 int delExistenceCost,
				 int delExtensionCost,
				 int insExistenceCost,
				 int insExtensionCost,
				 int maxScoreDrop,
				 int maxMatchScore,
				 const uchar *toUnmasked) {
  const SimdInt mNegInf = simdFill(-INF);
  const SimdInt mDelOpenCost = simdFill(delExistenceCost);
  const SimdInt mDelGrowCost = simdFill(delExtensionCost);
  const SimdInt mInsOpenCost = simdFill(insExistenceCost);
  const SimdInt mInsGrowCost = simdFill(insExtensionCost);
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

  init();
  seq1queue.clear();
  seq2queue.clear();

  bool isDelimiter1 = (*seq1 == delimiter);
  bool isDelimiter2 = (*seq2 == delimiter);
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

    initAntidiagonal(seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    Score *x0 = &xScores[thisPos];
    Score *y0 = &yScores[thisPos];
    Score *z0 = &zScores[thisPos];
    const Score *y1 = &yScores[horiPos];
    const Score *z1 = &zScores[horiPos + 1];
    const Score *x2 = &xScores[diagPos];

    if (isDelimiter1 || isDelimiter2) {
      updateMaxScoreDrop(maxScoreDrop, n, maxMatchScore);
    }

    int minScore = bestScore - maxScoreDrop;
    SimdInt mMinScore = simdFill(minScore);

    if (isDna) {
      for (int i = 0; i < numCells; i += simdLen) {
	__m128i fwd1 = simdLoad128(s1);
	__m128i fwd2 = simdLoad128(s2 - (seqLoadLen-1));
	__m128i rev2 = _mm_shuffle_epi8(fwd2, reverse);
	__m128i j = _mm_or_si128(_mm_slli_epi32(fwd1, 2), rev2);
	SimdInt s = simdBytesToInts(_mm_shuffle_epi8(scorer4x4, j));
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
    } else {
      for (int i = 0; i < numCells; i += simdLen) {
	SimdInt s = simdSet(
#ifdef __SSE4_1__
#ifdef __AVX2__
			    scorer[s1[7]][s2[-7]],
			    scorer[s1[6]][s2[-6]],
			    scorer[s1[5]][s2[-5]],
			    scorer[s1[4]][s2[-4]],
#endif
			    scorer[s1[3]][s2[-3]],
			    scorer[s1[2]][s2[-2]],
			    scorer[s1[1]][s2[-1]],
#endif
			    scorer[s1[0]][s2[-0]]);
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
    }

    int newBestScore = simdHorizontalMax(mBestScore);
    if (newBestScore > bestScore) {
      bestScore = newBestScore;
      bestAntidiagonal = antidiagonal;
    }

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[n] > -INF / 2) {
      ++seq1end;
      uchar x = toUnmasked[*seq1];
      seq1queue.push(x, n + seqLoadLen);
      seq1 += seqIncrement * (x != delimiter);
      uchar z = seq1queue.fromEnd(seqLoadLen);
      if (z >= 4) {
	isDna = false;
	isDelimiter1 = (z == delimiter);
      }
    }

    if (x0[0] > -INF / 2) {
      uchar y = toUnmasked[*seq2];
      seq2queue.push(y, n + seqLoadLen);
      seq2 += seqIncrement;
      if (y >= 4) {
	isDna = false;
	isDelimiter2 = (y == delimiter);
      }
    } else {
      ++seq1beg;
      ++diagPos;
      ++horiPos;
      if (seq1beg == seq1end) break;
    }
  }

  calcBestSeq1position(bestScore);
  return bestScore;
}

}
