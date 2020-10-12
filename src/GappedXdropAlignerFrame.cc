// Author: Martin C. Frith 2020
// SPDX-License-Identifier: GPL-3.0-or-later

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"
//#include <iostream>  // for debugging

namespace cbrc {

typedef const uchar *const_uchar_ptr;

// Puts 6 "dummy" antidiagonals at the start, so that we can safely
// look-back from subsequent antidiagonals.
void GappedXdropAligner::initFrame() {
  initAntidiagonal(0, 0, 0, 0);
  initAntidiagonal(1, 0, xdropPadLen, 0);
  initAntidiagonal(2, 0, xdropPadLen * 2, 0);
  initAntidiagonal(3, 0, xdropPadLen * 3, 0);
  initAntidiagonal(4, 0, xdropPadLen * 4, 0);
  initAntidiagonal(5, 0, xdropPadLen * 5, 0);
  xScores[xdropPadLen - 1] = 0;
  bestAntidiagonal = 0;
}

int GappedXdropAligner::alignFrame(const uchar *seq1,
				   const uchar *seq2frame0,
				   const uchar *seq2frame1,  // the +1 frame
				   const uchar *seq2frame2,  // the +2 frame
				   bool isForward,
				   const ScoreMatrixRow *scorer,
				   const GapCosts &gap,
				   int maxScoreDrop) {
  const const_uchar_ptr frames[] = {seq2frame0, seq2frame1, seq2frame2};
  const int seqIncrement = isForward ? 1 : -1;

  const int delOpenScore = -gap.delPieces[0].openCost;
  const int insOpenScore = -gap.insPieces[0].openCost;
  const int delScore1 = gap.delScore1;
  const int delScore2 = gap.delScore2;
  const int delScore3 = gap.delScore3;
  const int insScore1 = gap.insScore1;
  const int insScore2 = gap.insScore2;
  const int insScore3 = gap.insScore3;

  int runOfDrops = 2;
  int runOfEdges = 0;
  size_t seq1beg = 0;
  size_t seq1end = 1;
  size_t diagPos6 = xdropPadLen - 1;
  size_t horiPos5 = xdropPadLen * 2 - 1;
  size_t horiPos4 = xdropPadLen * 3 - 1;
  size_t horiPos3 = xdropPadLen * 4 - 1;
  size_t vertPos2 = xdropPadLen * 5;
  size_t vertPos1 = xdropPadLen * 6;
  size_t thisPos  = xdropPadLen * 6;

  int bestScore = 0;
  int newBestScore = 0;

  initFrame();

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    int numCells = seq1end - seq1beg;
    int n = numCells - 1;

    const uchar *seq2 = frames[antidiagonal % 3];
    size_t seq2pos = antidiagonal / 3 - seq1beg;
    const uchar *s1 = isForward ? seq1 + seq1beg : seq1 - seq1beg - 1;
    const uchar *s2 = isForward ? seq2 + seq2pos : seq2 - seq2pos - 1;

    initAntidiagonal(antidiagonal + 6, seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    Score *X0 = &xScores[thisPos];
    Score *Y0 = &yScores[thisPos];
    Score *Z0 = &zScores[thisPos];
    const Score *Z1 = &zScores[vertPos1];
    const Score *Z2 = &zScores[vertPos2];
    const Score *Y3 = &yScores[horiPos3];
    const Score *Z3 = &zScores[horiPos3 + 1];
    const Score *Y4 = &yScores[horiPos4];
    const Score *Y5 = &yScores[horiPos5];
    const Score *X6 = &xScores[diagPos6];

    int minScore = bestScore - maxScoreDrop;

    for (int i = 0; i < numCells; ++i) {
      int s = scorer[*s1][*s2];
      int y1 = Y5[i] + delScore1;
      int y2 = Y4[i] + delScore2;
      int y3 = Y3[i] + delScore3;
      int z1 = Z1[i] + insScore1;
      int z2 = Z2[i] + insScore2;
      int z3 = Z3[i] + insScore3;
      int b = maxValue(X6[i], y1, y2, y3, z1, z2, z3);
      bool isDrop = (b < minScore);
      newBestScore = maxValue(newBestScore, b);
      X0[i] = isDrop ? -INF : b + s;
      Y0[i] = maxValue(b + delOpenScore, y3);
      Z0[i] = maxValue(b + insOpenScore, z3);
      s1 += seqIncrement;
      s2 -= seqIncrement;
    }

    if (newBestScore > bestScore) {
      bestScore = newBestScore;
      bestAntidiagonal = antidiagonal;
    }

    diagPos6 = horiPos5;
    horiPos5 = horiPos4;
    horiPos4 = horiPos3;
    horiPos3 = vertPos2 - 1;
    vertPos2 = vertPos1;
    vertPos1 = thisPos;
    thisPos += numCells;

    if (X0[n] > -INF / 2 || runOfEdges) {
      ++runOfEdges;
      if (runOfEdges == 3) {
	++seq1end;
	runOfEdges = 0;
      }
    }

    if (X0[0] > -INF / 2) {
      runOfDrops = 0;
    } else {
      ++runOfDrops;
      if (runOfDrops == 3) {
	++seq1beg;
	if (seq1beg == seq1end) break;
	++diagPos6;
	++horiPos5;
	++horiPos4;
	++horiPos3;
	++vertPos2;
	++vertPos1;
	runOfDrops = 0;
      }
    }
  }

  calcBestSeq1position(bestScore, 6);
  return bestScore;
}

bool GappedXdropAligner::getNextChunkFrame(size_t &end1,
					   size_t &end2,
					   size_t &length,
					   int &gapCost,
					   const GapCosts &gap) {
  if (bestAntidiagonal == 0) return false;

  const int delOpenScore = -gap.delPieces[0].openCost;
  const int insOpenScore = -gap.insPieces[0].openCost;
  const int delScore1 = gap.delScore1;
  const int delScore2 = gap.delScore2;
  const int delScore3 = gap.delScore3;
  const int insScore1 = gap.insScore1;
  const int insScore2 = gap.insScore2;
  const int insScore3 = gap.insScore3;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - bestSeq1position * 3;

  int opt[7];
  int state = 0;

  while (1) {
    opt[0] = xScores[diag(bestAntidiagonal + 0, bestSeq1position)];
    opt[1] = yScores[hori(bestAntidiagonal + 0, bestSeq1position)] + delScore1;
    opt[2] = yScores[hori(bestAntidiagonal + 1, bestSeq1position)] + delScore2;
    opt[3] = yScores[hori(bestAntidiagonal + 2, bestSeq1position)] + delScore3;
    opt[4] = zScores[vert(bestAntidiagonal + 2, bestSeq1position)] + insScore3;
    opt[5] = zScores[vert(bestAntidiagonal + 3, bestSeq1position)] + insScore2;
    opt[6] = zScores[vert(bestAntidiagonal + 4, bestSeq1position)] + insScore1;
    state = std::max_element(opt, opt + 7) - opt;
    if (state != 0 || bestAntidiagonal == 0) break;
    bestAntidiagonal -= 6;
    bestSeq1position -= 1;
  }

  length = end1 - bestSeq1position;
  if (bestAntidiagonal == 0) {
    gapCost = 0;
    return true;
  }

  int scoreHere = opt[state];

  do {
    bool isDel = (state < 4);
    bestAntidiagonal -= 7 - state - isDel;
    if (isDel) bestSeq1position -= 1;
    opt[0] = xScores[diag(bestAntidiagonal + 0, bestSeq1position)];
    opt[1] = yScores[hori(bestAntidiagonal + 0, bestSeq1position)] + delScore1;
    opt[2] = yScores[hori(bestAntidiagonal + 1, bestSeq1position)] + delScore2;
    opt[3] = yScores[hori(bestAntidiagonal + 2, bestSeq1position)] + delScore3;
    opt[4] = zScores[vert(bestAntidiagonal + 2, bestSeq1position)] + insScore3;
    opt[5] = zScores[vert(bestAntidiagonal + 3, bestSeq1position)] + insScore2;
    opt[6] = zScores[vert(bestAntidiagonal + 4, bestSeq1position)] + insScore1;
    if (isDel) {
      opt[3] -= delOpenScore;
    } else {
      opt[4] -= insOpenScore;
    }
    state = std::max_element(opt, opt + 7) - opt;
  } while (state != 0);

  gapCost = opt[0] - scoreHere;
  return true;
}

}
