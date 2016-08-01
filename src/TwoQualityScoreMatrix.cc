// Copyright 2011 Martin C. Frith

#include "TwoQualityScoreMatrix.hh"

#include "qualityScoreUtil.hh"

#include <algorithm>  // min
#include <cassert>
#include <cmath>

namespace cbrc {

void TwoQualityMatrixIndexer::init(const uchar *toUnmasked) {
  indexMap.resize(qualityCapacity * scoreMatrixRowSize);

  for (int quality = 0; quality < qualityCapacity; ++quality) {
    int normalStart = quality * numQualityLetters;
    int maskedStart = normalStart + numNormalLetters;
    int abnormalPos = 0;

    for (int letter = 0; letter < scoreMatrixRowSize; ++letter) {
      int unmasked = toUnmasked[letter];
      int i = subindex(letter, quality);

      if (letter < numNormalLetters)
        indexMap[i] = normalStart + letter;
      else if (unmasked < numNormalLetters)
        indexMap[i] = maskedStart + unmasked;
      else
        indexMap[i] = abnormalPos++;
    }

    assert(abnormalPos <= minQuality * numQualityLetters);
  }
}

static int qualityEnd(const TwoQualityMatrixIndexer &indexer, int letter) {
  if (indexer(0, letter, 0, 0) == indexer(0, letter, 0, 1))
    return indexer.minQuality + 1;
  else
    return indexer.qualityCapacity;
}

void TwoQualityScoreMatrix::init(const ScoreMatrixRow *scoreMatrix,
                                 double lambda,
                                 const double *letterProbs1,
                                 const double *letterProbs2,
                                 bool isPhred1,
                                 int qualityOffset1,
                                 bool isPhred2,
                                 int qualityOffset2,
                                 const uchar *toUnmasked,
                                 bool isMask) {
  typedef TwoQualityMatrixIndexer Indexer;

  indexer.init(toUnmasked);
  data.resize(Indexer::numSymbols * Indexer::numSymbols);

  double certainties1[Indexer::numNormalLetters][Indexer::qualityCapacity];
  double certainties2[Indexer::numNormalLetters][Indexer::qualityCapacity];

  for (int q = Indexer::minQuality; q < Indexer::qualityCapacity; ++q) {
    double e1 = errorProbFromQual(q, qualityOffset1, isPhred1);
    double e2 = errorProbFromQual(q, qualityOffset2, isPhred2);
    for (int x = 0; x < Indexer::numNormalLetters; ++x) {
      certainties1[x][q] = qualityCertainty(e1, letterProbs1[x]);
      certainties2[x][q] = qualityCertainty(e2, letterProbs2[x]);
    }
  }

  for (int x1 = 0; x1 < scoreMatrixRowSize; ++x1) {
    for (int x2 = 0; x2 < scoreMatrixRowSize; ++x2) {
      bool isNormal1 = (x1 < Indexer::numNormalLetters);
      bool isNormal2 = (x2 < Indexer::numNormalLetters);
      bool isNormal = (isNormal1 && isNormal2);

      int unmasked1 = toUnmasked[x1];
      int unmasked2 = toUnmasked[x2];

      bool isMasked = (unmasked1 != x1 || unmasked2 != x2);

      int score = scoreMatrix[unmasked1][unmasked2];
      double expScore = std::exp(lambda * score);

      int end1 = qualityEnd(indexer, x1);
      int end2 = qualityEnd(indexer, x2);

      for (int q1 = Indexer::minQuality; q1 < end1; ++q1) {
        for (int q2 = Indexer::minQuality; q2 < end2; ++q2) {
	  if (isMasked) {
	    score = data[indexer(unmasked1, unmasked2, q1, q2)];
	    if (isMask) score = std::min(score, 0);
	  } else if (isNormal) {
            double c1 = certainties1[x1][q1];
            double c2 = certainties2[x2][q2];
            score = qualityPairScore(expScore, c1, c2, lambda);
          }
          data[indexer(x1, x2, q1, q2)] = score;
        }
      }
    }
  }
}

void TwoQualityExpMatrix::init(const TwoQualityScoreMatrix &m,
                               double temperature) {
  assert(temperature > 0);
  indexer = m.indexer;
  data.resize(indexer.numSymbols * indexer.numSymbols);

  for (int i1 = 0; i1 < scoreMatrixRowSize; ++i1) {
    for (int i2 = 0; i2 < scoreMatrixRowSize; ++i2) {
      int end1 = qualityEnd(indexer, i1);
      int end2 = qualityEnd(indexer, i2);

      for (int q1 = indexer.minQuality; q1 < end1; ++q1)
        for (int q2 = indexer.minQuality; q2 < end2; ++q2)
          data[indexer(i1, i2, q1, q2)] =
              std::exp(m(i1, i2, q1, q2) / temperature);
    }
  }
}

}
