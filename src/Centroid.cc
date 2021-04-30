// Copyright 2008, 2009, 2010, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#include "Centroid.hh"
#include "GappedXdropAlignerInl.hh"

#include <algorithm>
#include <cmath> // for exp
#include <cfloat>   // for DBL_MAX

static const double DINF = DBL_MAX / 2;

namespace{
  double EXP ( double x ) {
    return std::exp (x);
  }
}

namespace cbrc{

  ExpectedCount::ExpectedCount ()
  {
    double d0 = 0;
    toMatch = d0;
    delInit = d0;
    delNext = d0;
    insInit = d0;
    insNext = d0;

    for (int n=0; n<scoreMatrixRowSize; n++)
      for (int m=0; m<scoreMatrixRowSize; m++) emit[n][m] = d0;
  }

  void Centroid::setPssm( const ScoreMatrixRow* pssm, size_t qsize, double T,
			  const OneQualityExpMatrix& oqem,
			  const uchar* sequenceBeg, const uchar* qualityBeg ) {
    pssmExp.resize( qsize * scoreMatrixRowSize );
    pssmExp2 = reinterpret_cast<ExpMatrixRow*> ( &pssmExp[0] );

    if( oqem ){  // fast special case
      makePositionSpecificExpMatrix( oqem, sequenceBeg, sequenceBeg + qsize,
                                     qualityBeg, &pssmExp[0] );
    }
    else{  // slow general case
      for ( size_t i=0; i<qsize; ++i ) {
        for ( unsigned j=0; j<scoreMatrixRowSize; ++j ) {
          pssmExp2[ i ][ j ] = EXP ( pssm[ i ][ j ] / T );
        }
      }
    }
  }

  void Centroid::setLetterProbsPerPosition(unsigned alphabetSize,
					   size_t sequenceLength,
					   const uchar *sequence,
					   const uchar *qualityCodes,
					   bool isFastq,
					   const double *qualToProbCorrect,
					   const double *letterProbs,
					   const uchar *toUnmasked) {
    letterProbsPerPosition.resize(sequenceLength * alphabetSize);
    for (size_t i = 0; i < sequenceLength; ++i) {
      size_t j = i * alphabetSize;
      double *out = &letterProbsPerPosition[j];
      if (isFastq) {
	unsigned letter = toUnmasked[sequence[i]];
	if (letter < alphabetSize) {
	  double p = qualToProbCorrect[qualityCodes[i]];
	  for (unsigned k = 0; k < alphabetSize; ++k) {
	    out[k] = (1 - p) * letterProbs[k];
	  }
	  out[letter] = p * (1 - letterProbs[letter]);
	  // it's OK to scale the "out" values by a constant, per "i"
	} else {
	  std::fill_n(out, alphabetSize, 0.0);
	}
      } else {
	for (unsigned k = 0; k < alphabetSize; ++k) {
	  out[k] = qualToProbCorrect[qualityCodes[j + k]];
	}
      }
    }
  }

  double Centroid::forward(const uchar *seq1, const uchar *seq2,
			   size_t start2, bool isExtendFwd,
			   const const_dbl_ptr *substitutionProbs,
			   const GapCosts &gapCosts, int globality) {
    const ExpMatrixRow *pssm = pssmExp.empty() ? 0 : pssmExp2 + start2;
    const int seqIncrement = isExtendFwd ? 1 : -1;

    const double delInit = gapCosts.delProbPieces[0].openProb;
    const double delNext = gapCosts.delProbPieces[0].growProb;
    const double insInit = gapCosts.insProbPieces[0].openProb;
    const double insNext = gapCosts.insProbPieces[0].growProb;

    initForward();

    size_t thisPos = xdropPadLen * 2;

    double Z = 0.0;  // partion function of forward values

    for (size_t antidiagonal = 0; antidiagonal < numAntidiagonals; ++antidiagonal) {
      const size_t seq1beg = xa.seq1start(antidiagonal);
      const size_t seq2pos = antidiagonal - seq1beg;
      const double scale2 = rescales[antidiagonal];
      const double scale1 = rescales[antidiagonal + 1];
      double sum_f = 0.0; // sum of forward values

      const double scaledDelNext = scale1 * delNext;
      const double scaledInsNext = scale1 * insNext;

      double *fM0 = &fM[thisPos];
      double *fD0 = &fD[thisPos];
      double *fI0 = &fI[thisPos];

      const size_t horiPos = xa.hori(antidiagonal, seq1beg);
      const size_t vertPos = xa.vert(antidiagonal, seq1beg);
      const size_t diagPos = xa.diag(antidiagonal, seq1beg);
      const double *fD1 = &fD[horiPos];
      const double *fI1 = &fI[vertPos];
      const double *fM2 = &fM[diagPos];

      const size_t nextPos = xa.scoreEndIndex(antidiagonal + 1);
      const int numCells = nextPos - thisPos - xdropPadLen;
      const uchar* s1 = isExtendFwd ? seq1 + seq1beg : seq1 - seq1beg;

      for (int i = 0; i < xdropPadLen; ++i) {
	*fM0++ = *fD0++ = *fI0++ = 0.0;
      }

      if (!pssm) {
	const uchar* s2 = isExtendFwd ? seq2 + seq2pos : seq2 - seq2pos;
	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = substitutionProbs[*s1][*s2];
	  const double xM = fM2[i];
	  const double xD = fD1[i];
	  const double xI = fI1[i];
	  const double xSum = (xM * scale2 + xD + xI) * scale1;
	  fD0[i] = xSum * delInit + xD * scaledDelNext;
	  fI0[i] = xSum * insInit + xI * scaledInsNext;
	  fM0[i] = xSum * matchProb;
	  sum_f += xSum;
	  if (globality && matchProb <= 0) Z += xSum;  // xxx
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}
      } else {
	const ExpMatrixRow* p2 = isExtendFwd ? pssm + seq2pos : pssm - seq2pos;
	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = (*p2)[*s1];
	  const double xM = fM2[i];
	  const double xD = fD1[i];
	  const double xI = fI1[i];
	  const double xSum = (xM * scale2 + xD + xI) * scale1;
	  fD0[i] = xSum * delInit + xD * scaledDelNext;
	  fI0[i] = xSum * insInit + xI * scaledInsNext;
	  fM0[i] = xSum * matchProb;
	  sum_f += xSum;
	  if (globality && matchProb <= 0) Z += xSum;  // xxx
	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	}
      }

      thisPos = nextPos;

      if( !globality ) Z += sum_f;
      rescales[antidiagonal + 2] = 1.0 / (sum_f + 1.0);  // seems ugly
      Z *= rescales[antidiagonal + 2];
    }

    assert( Z > 0.0 );
    rescales[numAntidiagonals + 1] /= Z;  // this causes scaled Z to equal 1
    return logPartitionFunction();
  }

  // added by M. Hamada
  // compute posterior probabilities while executing backward algorithm
  void Centroid::backward(const uchar* seq1, const uchar* seq2,
			  size_t start2, bool isExtendFwd,
			  const const_dbl_ptr *substitutionProbs,
			  const GapCosts& gapCosts, int globality) {
    const ExpMatrixRow *pssm = pssmExp.empty() ? 0 : pssmExp2 + start2;
    const int seqIncrement = isExtendFwd ? 1 : -1;

    const double delInit = gapCosts.delProbPieces[0].openProb;
    const double delNext = gapCosts.delProbPieces[0].growProb;
    const double insInit = gapCosts.insProbPieces[0].openProb;
    const double insNext = gapCosts.insProbPieces[0].growProb;

    size_t oldPos = xa.scoreEndIndex(numAntidiagonals);
    initBackward(oldPos);

    double scaledUnit = 1.0;

    for (size_t antidiagonal = numAntidiagonals; antidiagonal-- > 0;) {
      const size_t seq1beg = xa.seq1start(antidiagonal);
      const size_t seq2pos = antidiagonal - seq1beg;
      const double scale2 = rescales[antidiagonal];
      const double scale1 = rescales[antidiagonal + 1];
      scaledUnit *= rescales[antidiagonal + 2];

      const double scaledDelNext = scale1 * delNext;
      const double scaledInsNext = scale1 * insNext;

      const size_t newPos = xa.scoreEndIndex(antidiagonal);
      const double *bM0 = &bM[newPos + xdropPadLen];
      const double *bD0 = &bD[newPos + xdropPadLen];
      const double *bI0 = &bI[newPos + xdropPadLen];

      const size_t horiPos = xa.hori(antidiagonal, seq1beg);
      const size_t vertPos = xa.vert(antidiagonal, seq1beg);
      const size_t diagPos = xa.diag(antidiagonal, seq1beg);
      double *bD1 = &bD[horiPos];
      double *bI1 = &bI[vertPos];
      double *bM2 = &bM[diagPos];

      const double *fD1 = &fD[horiPos];
      const double *fI1 = &fI[vertPos];

      double* mDout = &mD[ seq1beg ];
      double* mIout = &mI[ seq2pos ];

      const int numCells = oldPos - newPos - xdropPadLen;
      const uchar *s1 = isExtendFwd ? seq1 + seq1beg : seq1 - seq1beg;

      if (!pssm) {
	const uchar *s2 = isExtendFwd ? seq2 + seq2pos : seq2 - seq2pos;

	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = substitutionProbs[*s1][*s2];
	  const double yM = bM0[i];
	  const double yD = bD0[i];
	  const double yI = bI0[i];
	  double ySum = yM * matchProb + yD * delInit + yI * insInit;

	  if (!globality || matchProb <= 0) ySum += scaledUnit;
	  // xxx matchProb should be 0 only at delimiters, but will be
	  // 0 for non-delimiters with severe mismatch scores

	  ySum *= scale1;

	  bM2[i] = ySum * scale2;
	  bD1[i] = ySum + yD * scaledDelNext;
	  bI1[i] = ySum + yI * scaledInsNext;

	  *mDout += fD1[i] * bD1[i];
	  *mIout += fI1[i] * bI1[i];

	  mDout++; mIout--;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}
      } else {
	const ExpMatrixRow *p2 = isExtendFwd ? pssm + seq2pos : pssm - seq2pos;

	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = (*p2)[*s1];
	  const double yM = bM0[i];
	  const double yD = bD0[i];
	  const double yI = bI0[i];
	  double ySum = yM * matchProb + yD * delInit + yI * insInit;

	  if (!globality || matchProb <= 0) ySum += scaledUnit;  // xxx

	  ySum *= scale1;

	  bM2[i] = ySum * scale2;
	  bD1[i] = ySum + yD * scaledDelNext;
	  bI1[i] = ySum + yI * scaledInsNext;

	  *mDout += fD1[i] * bD1[i];
	  *mIout += fI1[i] * bI1[i];

	  mDout++; mIout--;
	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	}
      }

      oldPos = newPos;
    }
  }

  double Centroid::dp_centroid( double gamma ){
    for( size_t k = 1; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      size_t seq1pos = xa.seq1start( k );

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const size_t h = xa.hori( k, seq1pos );
      const size_t d = xa.diag( k, seq1pos );
      const double* X1 = &X[h];
      const double* X2 = &X[d];
      const double* fM2 = &fM[d];
      const double* bM2 = &bM[d];

      for (int i = 0; i < xdropPadLen; ++i) {
	*X0++ = -DINF;
      }

      do{
	const double matchProb = (*fM2++) * (*bM2++);
	const double s = ( gamma + 1 ) * matchProb - 1;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max( std::max( oldX1, *X1 ), *X2++ + s );
	//assert ( score >= 0 );
	updateScore ( score, k, seq1pos );
	*X0++ = score;
	seq1pos++;
      }while( X0 != x0end );
    }
    return bestScore;
  }

  void Centroid::traceback_centroid( std::vector< SegmentPair >& chunks,
				     double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const size_t h = xa.hori( k, i );
      const size_t v = xa.vert( k, i );
      const size_t d = xa.diag( k, i );
      const double matchProb = fM[d] * bM[d];
      const int m = maxIndex( X[d] + (gamma + 1) * matchProb - 1, X[h], X[v] );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  double Centroid::dp_ama( double gamma ){
    mX1.assign ( numAntidiagonals + 2, 1.0 );
    mX2.assign ( numAntidiagonals + 2, 1.0 );

    for (size_t k = 0; k < numAntidiagonals; ++k) {
      size_t seq1pos = xa.seq1start(k);
      size_t seq2pos = k - seq1pos;
      size_t loopBeg = xa.diag(k, seq1pos);
      size_t loopEnd = loopBeg + xa.numCellsAndPads(k) - xdropPadLen;
      for (size_t i = loopBeg; i < loopEnd; ++i) {
	const double matchProb = fM[i] * bM[i];
	mX1[seq1pos++] -= matchProb;
	mX2[seq2pos--] -= matchProb;
      }
    }

    for( size_t k = 1; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      size_t seq1pos = xa.seq1start( k );
      size_t seq2pos = k - seq1pos;

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const size_t h = xa.hori( k, seq1pos );
      const size_t d = xa.diag( k, seq1pos );
      const double* X1 = &X[h];
      const double* X2 = &X[d];
      const double* fM2 = &fM[d];
      const double* bM2 = &bM[d];

      for (int i = 0; i < xdropPadLen; ++i) {
	*X0++ = -DINF;
      }

      do{
	const double matchProb = (*fM2++) * (*bM2++);
	const double thisD = mD[seq1pos];
	const double thisI = mI[seq2pos];
	const double thisXD = mX1[seq1pos] - thisD;
	const double thisXI = mX2[seq2pos] - thisI;
	const double s = 2 * gamma * matchProb - (thisXD + thisXI);
	const double u = gamma * thisD - thisXD;
	const double t = gamma * thisI - thisXI;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max(std::max(oldX1 + u, *X1 + t), *X2++ + s);
	updateScore ( score, k, seq1pos );
	*X0++ = score;
	seq1pos++;
	seq2pos--;
      }while( X0 != x0end );
    }

    return bestScore;
  }

  void Centroid::traceback_ama( std::vector< SegmentPair >& chunks,
			    double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const size_t j = k - i;
      const size_t h = xa.hori( k, i );
      const size_t v = xa.vert( k, i );
      const size_t d = xa.diag( k, i );
      const double matchProb = fM[d] * bM[d];
      const double thisD = mD[i];
      const double thisI = mI[j];
      const double thisXD = mX1[i] - thisD;
      const double thisXI = mX2[j] - thisI;
      const double s = 2 * gamma * matchProb - (thisXD + thisXI);
      const double t = gamma * thisI - thisXI;
      const double u = gamma * thisD - thisXD;
      const int m = maxIndex( X[d] + s, X[h] + u, X[v] + t );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  void Centroid::getMatchAmbiguities(std::vector<char>& ambiguityCodes,
				     size_t seq1end, size_t seq2end,
				     size_t length) const {
    while (length) {
      size_t d = xa.diag(seq1end + seq2end, seq1end);
      double p = fM[d] * bM[d];
      ambiguityCodes.push_back(asciiProbability(p));
      --seq1end;  --seq2end;  --length;
    }
  }

  void Centroid::getDeleteAmbiguities(std::vector<char>& ambiguityCodes,
				      size_t seq1end, size_t seq1beg) const {
    for (size_t i = seq1end; i > seq1beg; --i)
      ambiguityCodes.push_back(asciiProbability(mD[i]));
  }

  void Centroid::getInsertAmbiguities(std::vector<char>& ambiguityCodes,
				      size_t seq2end, size_t seq2beg) const {
    for (size_t i = seq2end; i > seq2beg; --i)
      ambiguityCodes.push_back(asciiProbability(mI[i]));
  }

  double Centroid::logPartitionFunction() const{
    double x = 0.0;
    for( size_t k = 0; k < numAntidiagonals; ++k ){
      x -= std::log(rescales[k + 2]);
    }
    return x;
  }

  static void countUncertainLetters(double *counts, double alignProb,
				    unsigned alphabetSize,
				    const double *probRatios,
				    const double *letterProbs) {
    double ratioParts[scoreMatrixRowSize];
    double sum = 0;
    for (unsigned letter = 0; letter < alphabetSize; ++letter) {
      double r = probRatios[letter] * letterProbs[letter];
      ratioParts[letter] = r;
      sum += r;
    }
    if (sum > 0) {
      const double mul = alignProb / sum;
      for (unsigned letter = 0; letter < alphabetSize; ++letter) {
	counts[letter] += mul * ratioParts[letter];
      }
    }
  }

  void Centroid::computeExpectedCounts(const uchar *seq1, const uchar *seq2,
				       size_t start2, bool isExtendFwd,
				       const const_dbl_ptr *substitutionProbs,
				       const GapCosts &gapCosts,
				       unsigned alphabetSize,
				       ExpectedCount &c) const {
    const ExpMatrixRow *pssm = pssmExp.empty() ? 0 : pssmExp2 + start2;

    const double *letterProbs = 0;
    if (!letterProbsPerPosition.empty()) {
      letterProbs = &letterProbsPerPosition[0] + start2 * alphabetSize;
    }

    const int seqIncrement = isExtendFwd ? 1 : -1;
    int alphabetSizeIncrement = alphabetSize;
    if (!isExtendFwd) alphabetSizeIncrement *= -1;

    size_t thisPos = xdropPadLen * 3;

    double delInitCount = 0; double delNextCount = 0;
    double insInitCount = 0; double insNextCount = 0;

    for (size_t antidiagonal = 0; antidiagonal < numAntidiagonals;) {
      const size_t seq1beg = xa.seq1start(antidiagonal);
      const size_t seq2pos = antidiagonal - seq1beg;
      const double scale2 = rescales[antidiagonal];
      const double scale1 = rescales[antidiagonal + 1];

      const double *bM0 = &bM[thisPos];
      const double *bD0 = &bD[thisPos];
      const double *bI0 = &bI[thisPos];

      const size_t horiPos = xa.hori(antidiagonal, seq1beg);
      const size_t vertPos = xa.vert(antidiagonal, seq1beg);
      const size_t diagPos = xa.diag(antidiagonal, seq1beg);
      const double *fD1 = &fD[horiPos];
      const double *fI1 = &fI[vertPos];
      const double *fM2 = &fM[diagPos];

      double dNextCount = 0;
      double iNextCount = 0;

      ++antidiagonal;
      const size_t nextPos = xa.scoreEndIndex(antidiagonal);
      const int numCells = nextPos - thisPos;
      const uchar* s1 = isExtendFwd ? seq1 + seq1beg : seq1 - seq1beg;

      if (!pssm) {
	const uchar* s2 = isExtendFwd ? seq2 + seq2pos : seq2 - seq2pos;

	for (int i = 0; i < numCells; ++i) {
	  const unsigned letter1 = *s1;
	  const unsigned letter2 = *s2;
	  const double matchProb = substitutionProbs[letter1][letter2];

	  const double yM = bM0[i];
	  const double yD = bD0[i];
	  const double yI = bI0[i];

	  const double xM = fM2[i];
	  const double xD = fD1[i];
	  const double xI = fI1[i];
	  const double xSum = (xM * scale2 + xD + xI) * scale1;

	  const double alignProb = xSum * yM * matchProb;
	  c.emit[letter1][letter2] += alignProb;
	  c.toMatch += alignProb;
	  delInitCount += xSum * yD;
	  dNextCount += xD * yD;
	  insInitCount += xSum * yI;
	  iNextCount += xI * yI;

	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}
      } else {
	const ExpMatrixRow* p2 = isExtendFwd ? pssm + seq2pos : pssm - seq2pos;
	const size_t a2 = seq2pos * alphabetSize;
	const double* lp2 = isExtendFwd ? letterProbs + a2 : letterProbs - a2;

	for (int i = 0; i < numCells; ++i) {
	  const unsigned letter1 = *s1;
	  const double matchProb = (*p2)[letter1];

	  const double yM = bM0[i];
	  const double yD = bD0[i];
	  const double yI = bI0[i];

	  const double xM = fM2[i];
	  const double xD = fD1[i];
	  const double xI = fI1[i];
	  const double xSum = (xM * scale2 + xD + xI) * scale1;

	  const double alignProb = xSum * yM * matchProb;
	  countUncertainLetters(c.emit[letter1], alignProb,
				alphabetSize, substitutionProbs[letter1], lp2);
	  c.toMatch += alignProb;
	  delInitCount += xSum * yD;
	  dNextCount += xD * yD;
	  insInitCount += xSum * yI;
	  iNextCount += xI * yI;

	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	  lp2 -= alphabetSizeIncrement;
	}
      }

      delNextCount += dNextCount * scale1;
      insNextCount += iNextCount * scale1;

      thisPos = nextPos + xdropPadLen;
    }

    const double delInit = gapCosts.delProbPieces[0].openProb;
    const double delNext = gapCosts.delProbPieces[0].growProb;
    const double insInit = gapCosts.insProbPieces[0].openProb;
    const double insNext = gapCosts.insProbPieces[0].growProb;

    c.delInit += delInitCount * delInit;
    c.delNext += delNextCount * delNext;
    c.insInit += insInitCount * insInit;
    c.insNext += insNextCount * insNext;
  }

}  // end namespace cbrc
