// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_substitution_matrix_stats.hh"
#include "LambdaCalculator.hh"
#include "cbrc_linalg.hh"
#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <math.h>

static double checkedExp(double lambda, int score) {
  double y = exp(lambda * score);
  if (!(y < HUGE_VAL)) throw std::overflow_error("exp overflow");
  return y;
}

static void permuteComplement(std::vector<double> &v,
			      const unsigned char *complement) {
  for (unsigned i = 0; i < v.size(); ++i) {
    unsigned j = complement[i];
    if (j < i) std::swap(v[i], v[j]);
  }
}

static double calcLetterProbs(std::vector<double> &probs, unsigned size,
			      double **expMat) {
  probs.assign(size, 1.0);
  cbrc::linalgSolve(expMat, &probs[0], size);
  double sum = 0;
  for (unsigned i = 0; i < size; ++i) {
    if (probs[i] < 0) {
      throw std::runtime_error("got a probability < 0 for the "
			       "substitution matrix");
    }
    sum += probs[i];
  }
  assert(sum > 0);
  double bias = 1 / sum;
  for (unsigned i = 0; i < size; ++i) {
    probs[i] *= bias;
  }
  return bias;
}

namespace mcf {

void SubstitutionMatrixStats::calcFromScale(const const_int_ptr *scoreMatrix,
					    unsigned size, double scale) {
  assert(size > 0);
  assert(scale > 0);
  mLambda = 1 / scale;

  std::vector<double> expVec(size * size);
  std::vector<double *> expMat(size);
  for (unsigned i = 0; i < size; ++i)
    expMat[i] = &expVec[i * size];

  for (unsigned i = 0; i < size; ++i)
    for (unsigned j = 0; j < size; ++j)
      expMat[i][j] = checkedExp(mLambda, scoreMatrix[j][i]);
  double bias1 = calcLetterProbs(mLetterProbs1, size, &expMat[0]);

  for (unsigned i = 0; i < size; ++i)
    for (unsigned j = 0; j < size; ++j)
      expMat[i][j] = checkedExp(mLambda, scoreMatrix[i][j]);
  double bias2 = calcLetterProbs(mLetterProbs2, size, &expMat[0]);

  // bias1 and bias2 should be equal
  mBias = (bias1 + bias2) / 2;
}

void SubstitutionMatrixStats::calcUnbiased(const const_int_ptr *scoreMatrix,
					   unsigned size) {
  cbrc::LambdaCalculator c;
  c.calculate(scoreMatrix, size);
  mLambda = c.lambda();
  if (!isBad()) {
    mBias = 1;
    mLetterProbs1.assign(c.letterProbs1(), c.letterProbs1() + size);
    mLetterProbs2.assign(c.letterProbs2(), c.letterProbs2() + size);
  }
}

void SubstitutionMatrixStats::flipDnaStrands(const unsigned char *complement) {
  permuteComplement(mLetterProbs1, complement);
  permuteComplement(mLetterProbs2, complement);
}

}
