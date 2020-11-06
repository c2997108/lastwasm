// Author: Martin C. Frith 2020
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef FRAMESHIFT_XDROP_ALIGNER_HH
#define FRAMESHIFT_XDROP_ALIGNER_HH

#include "mcf_gap_costs.hh"

#include <stddef.h>

#include <vector>

namespace mcf {

typedef const double *const_dbl_ptr;
typedef double *dbl_ptr;

class FrameshiftXdropAligner {
  typedef unsigned char uchar;
  typedef const uchar *const_uchar_ptr;
  typedef double Prob;

  enum { padSize = 1 };
  enum { rescaleStep = 32 };

public:
  double forward(const uchar *protein,
		 const uchar *frame0,  // translated DNA,  0 frame
		 const uchar *frame1,  // translated DNA, +1 frame
		 const uchar *frame2,  // translated DNA, +2 frame
		 bool isRightwardExtension,
		 const const_dbl_ptr *substitutionProbs,
		 const GapCosts &gapCosts,
		 double probDropLimit);

  void backward(bool isRightwardExtension,
		const const_dbl_ptr *substitutionProbs,
		const GapCosts &gapCosts);

  void count(bool isRightwardExtension,
	     const GapCosts &gapCosts,
	     const dbl_ptr *substitutionCounts,
	     double *transitionCounts);

  double matchProb(size_t proteinCoordinate, size_t dnaCoordinate) const {
    size_t antidiagonal = proteinCoordinate * 3 + dnaCoordinate;
    if (antidiagonal + 6 >= numOfAntidiagonals) return 0;
    size_t a = antidiagonal * 2;
    size_t i = xdropShape[a + 2] - xdropShape[a + 1] + proteinCoordinate;
    if (i < xdropShape[a] + padSize || i >= xdropShape[a + 2]) return 0;
    size_t b = (antidiagonal + 6) * 2;
    size_t j = xdropShape[b + 2] - xdropShape[b + 1] + proteinCoordinate + 1;
    size_t totalNumOfCells = xdropShape[numOfAntidiagonals * 2];
    return xFwdProbs[i] * xBckProbs[totalNumOfCells - j + padSize*7 - 1];
  }

private:
  std::vector<Prob> xFwdProbs;
  std::vector<Prob> yFwdProbs;
  std::vector<Prob> zFwdProbs;

  std::vector<Prob> xBckProbs;
  std::vector<Prob> yBckProbs;
  std::vector<Prob> zBckProbs;

  std::vector<size_t> xdropShape;

  std::vector<double> rescales;

  double rescaledSumOfProbRatios;

  size_t numOfAntidiagonals;

  const uchar *proteinPtr;
  const_uchar_ptr frames[3];

  size_t begInProbs(size_t antidiagonal) const
  { return xdropShape[antidiagonal * 2]; }

  size_t endInProtein(size_t antidiagonal) const
  { return xdropShape[antidiagonal * 2 + 1]; }

  int numOfCellsInAntidiagonal(size_t antidiagonal) const
  { return begInProbs(antidiagonal + 1) - begInProbs(antidiagonal) - padSize; }

  void resizeFwdProbsIfSmaller(size_t size) {
    if (xFwdProbs.size() < size) {
      xFwdProbs.resize(size);
      yFwdProbs.resize(size);
      zFwdProbs.resize(size);
    }
  }

  void resizeBckProbsIfSmaller(size_t size) {
    if (xBckProbs.size() < size) {
      xBckProbs.resize(size);
      yBckProbs.resize(size);
      zBckProbs.resize(size);
    }
  }

  void rescaleFwdProbs(size_t beg, size_t end, double scale) {
    for (size_t i = beg; i < end; ++i) {
      xFwdProbs[i] *= scale;
      yFwdProbs[i] *= scale;
      zFwdProbs[i] *= scale;
    }
  }

  void rescaleBckProbs(size_t beg, size_t end, double scale) {
    for (size_t i = beg; i < end; ++i) {
      xBckProbs[i] *= scale;
      yBckProbs[i] *= scale;
      zBckProbs[i] *= scale;
    }
  }
};

}

#endif
