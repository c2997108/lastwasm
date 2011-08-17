// Copyright 2011 Martin C. Frith

// The algorithm is based on these recurrence formulas, for
// generalized affine gap costs.  For standard affine gap costs, set
// gup=infinity.
//
// gop = gapExistenceCost
// gep = gapExtensionCost
// gup = gapUnalignedCost
//
// The 1st sequence: s(1), s(2), s(3), ...
// The 2nd sequence: t(1), t(2), t(3), ...
//
// matchScore(i, j)  =  the score for aligning s(i) with t(j).
//
// Initialization:
// x(i, 0)  =  y(i, 0)  =  z(i, 0)  =  -INF  (for all i >= 0)
// x(0, j)  =  y(0, j)  =  z(0, j)  =  -INF  (for all j >= 0)
// x(0, 0)  =  0
//
// Recurrence (i > 0 and j > 0):
// X(i, j)  =  x(i-1, j-1)
// Y(i, j)  =  max[ y(i-1, j) - gep, y(i-1, j-1) - gup ]
// Z(i, j)  =  max[ z(i, j-1) - gep, z(i-1, j-1) - gup ]
// b(i, j)  =  max[ X(i, j), Y(i, j), Z(i, j) ]
// x(i, j)  =  b(i, j) + matchScore(i, j)
// y(i, j)  =  max[ b(i, j) - gop, Y(i, j) ]
// z(i, j)  =  max[ b(i, j) - gop, Z(i, j) ]
//
// Interpretations:
// X(i, j)  =  the best score for any alignment ending with s(i-1)
//             aligned to t(j-1).
// b(i, j)  =  the best score for any alignment ending at s(i-1) and
//             t(j-1).
// x(i, j)  =  the best score for any alignment ending with s(i)
//             aligned to t(j).

// The recurrences are calculated antidiagonal-by-antidiagonal, where:
// antidiagonal  =  i + j

// We store x(i, j), y(i, j), and z(i, j) in the following way.
// xScores: oxx2x33x444x5555x66666...
// yScores: xxx2x33x444x5555x66666...
// zScores: xxx2x33x444x5555x66666...
// "o" indicates a cell with score = 0.
// "x" indicates a pad cell with score = -INF.
// "2", "3", etc. indicate cells in antidiagonal 2, 3, etc.

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

namespace cbrc {

// Puts 2 "dummy" antidiagonals at the start, so that we can safely
// look-back from subsequent antidiagonals.
void GappedXdropAligner::init() {
  seq1starts.resize(0);
  scoreEnds.resize(1);

  initAntidiagonal(0, 0, 0);
  xScores[0] = 0;
  yScores[0] = -INF;
  zScores[0] = -INF;

  initAntidiagonal(0, 1, 0);
  xScores[1] = -INF;
  yScores[1] = -INF;
  zScores[1] = -INF;

  bestAntidiagonal = 2;
}

void GappedXdropAligner::initAntidiagonal(std::size_t seq1beg,
                                          std::size_t scoreEnd,
                                          std::size_t numCells) {
  std::size_t newEnd = scoreEnd + numCells + 1;  // + 1 pad cell

  if (xScores.size() < newEnd) {
    xScores.resize(newEnd);
    yScores.resize(newEnd);
    zScores.resize(newEnd);
  }

  seq1starts.push_back(seq1beg);
  scoreEnds.push_back(newEnd);
}

int GappedXdropAligner::align(const uchar *seq1,
                              const uchar *seq2,
                              bool isForward,
                              const ScoreMatrixRow *scorer,
                              int gapExistenceCost,
                              int gapExtensionCost,
                              int gapUnalignedCost,
                              int maxScoreDrop,
                              int maxMatchScore) {
  bool isAffine = gapUnalignedCost >= gapExistenceCost + 2 * gapExtensionCost;

  std::size_t maxSeq1begs[] = { 0, 9 };
  std::size_t minSeq1ends[] = { 1, 0 };

  int bestScore = 0;

  init();

  for (std::size_t antidiagonal = 2; /* noop */; ++antidiagonal) {
    std::size_t seq1beg = arrayMin(maxSeq1begs);
    std::size_t seq1end = arrayMax(minSeq1ends);

    if (seq1beg >= seq1end) break;

    std::size_t scoreEnd = scoreEnds.back();
    std::size_t numCells = seq1end - seq1beg;

    initAntidiagonal(seq1beg, scoreEnd, numCells);

    std::size_t seq2pos = antidiagonal - 2 - seq1beg;

    const uchar *s1 = isForward ? seq1 + seq1beg : seq1 - seq1beg - 1;
    const uchar *s2 = isForward ? seq2 + seq2pos : seq2 - seq2pos - 1;

    if (isDelimiter(*s2, *scorer))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    int minScore = bestScore - maxScoreDrop;

    int *x0 = &xScores[scoreEnd];
    int *y0 = &yScores[scoreEnd];
    int *z0 = &zScores[scoreEnd];
    const int *y1 = &yScores[hori(antidiagonal, seq1beg)];
    const int *z1 = &zScores[vert(antidiagonal, seq1beg)];
    const int *x2 = &xScores[diag(antidiagonal, seq1beg)];

    const int *x0last = x0 + numCells;

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0base = x0 - seq1beg;

    if (isAffine) {
      // We could avoid this code duplication, by using: transposed(scorer).
      if (isForward)
        while (1) {
          int x = *x2;
          int y = *y1 - gapExtensionCost;
          int z = *z1 - gapExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + scorer[*s1][*s2];
            int g = b - gapExistenceCost;
            *y0 = maxValue(g, y);
            *z0 = maxValue(g, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  --s2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
      else
        while (1) {
          int x = *x2;
          int y = *y1 - gapExtensionCost;
          int z = *z1 - gapExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + scorer[*s1][*s2];
            int g = b - gapExistenceCost;
            *y0 = maxValue(g, y);
            *z0 = maxValue(g, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  ++s2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
    } else {
      const int *y2 = &yScores[diag(antidiagonal, seq1beg)];
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      while (1) {
        int x = *x2;
        int y = maxValue(*y1 - gapExtensionCost, *y2 - gapUnalignedCost);
        int z = maxValue(*z1 - gapExtensionCost, *z2 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
          updateBest(bestScore, b, antidiagonal, x0, x0base);
          *x0 = b + scorer[*s1][*s2];
          int g = b - gapExistenceCost;
          *y0 = maxValue(g, y);
          *z0 = maxValue(g, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;  ++y2;  ++z2;
        if (isForward) { ++s1;  --s2; }
        else           { --s1;  ++s2; }
      }
    }

    if (isDelimiter(*s1, *scorer))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    updateFiniteEdges(maxSeq1begs, minSeq1ends, x0base, x0 + 1, numCells);
  }

  return bestScore;
}

bool GappedXdropAligner::getNextChunk(std::size_t &end1,
                                      std::size_t &end2,
                                      std::size_t &length,
                                      int gapExistenceCost,
                                      int gapExtensionCost,
                                      int gapUnalignedCost) {
  if (bestAntidiagonal == 2) return false;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - 2 - bestSeq1position;
  length = 0;

  int state = 0;

  while (1) {
    if (state < 1 || state > 2) bestAntidiagonal -= 2;
    else                        bestAntidiagonal -= 1;

    if (state != 2) bestSeq1position -= 1;

    assert(bestAntidiagonal >= 2);
    assert(bestSeq1position <= bestAntidiagonal - 2);

    std::size_t h = hori(bestAntidiagonal, bestSeq1position);
    std::size_t v = vert(bestAntidiagonal, bestSeq1position);
    std::size_t d = diag(bestAntidiagonal, bestSeq1position);

    int x = xScores[d];
    int y = yScores[h] - gapExtensionCost;
    int z = zScores[v] - gapExtensionCost;
    int a = yScores[d] - gapUnalignedCost;
    int b = zScores[d] - gapUnalignedCost;

    if (state == 1) {
      y += gapExistenceCost;
      a += gapExistenceCost;
    }

    if (state == 2) {
      z += gapExistenceCost;
      b += gapExistenceCost;
    }

    state = maxIndex(x, y, z, a, b);

    if (length == 0 && (state > 0 || bestAntidiagonal == 2))
      length = end1 - bestSeq1position;

    if (length > 0 && state == 0) return true;
  }
}

}
