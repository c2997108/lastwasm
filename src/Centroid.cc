// Copyright 2008 Michiaki Hamada

#include "Centroid.hh"
#include <algorithm>  // max_element
#include <cassert>
#include <cmath> // for exp
//#include <iostream>  // for debugging
#include <string>
#include <cfloat>   // for DBL_MAX

#define CI(type) std::vector<type>::const_iterator  // added by MCF

// We should use "logadd" for avoiding overflow althogh
// it is slower than "exp"
static const std::string stype = "logadd"; // or exp
// log (0)
static const double LOG_ZERO = - DBL_MAX / 2;
static const double DINF = DBL_MAX / 2;

namespace{

  double EXP ( double x ) {
    return std::exp (x);
  }

  double LOG_ADD( double x, double y ){
    if ( y < x ) 
      return x + std::log( 1 + std::exp( y - x ) );
    return y + std::log( 1 + std::exp ( x - y ) );
  }

  double LOG_ADD( double x, double y, double z ){
    return LOG_ADD( x, LOG_ADD( y, z ) );
  }

  double LOG_ADD( double x, double y, double z, double w ){
    return LOG_ADD( x, LOG_ADD( y, LOG_ADD( z, w ) ) );
  }

}


namespace cbrc{

  // get DP matrix value "diagonal from" the given position
  // return 0 if there does not exist diagonal value
  double Centroid::diag( const dmatrix_t& matrix,
			 size_t antiDiagonal, size_t seq1pos ) const{
    if( antiDiagonal > 1 &&
	seq1pos > xa.fillBeg( antiDiagonal-2 ) &&
	seq1pos <= xa.fillEnd( antiDiagonal-2 ) ){
      return xa.cell( matrix, antiDiagonal-2, seq1pos-1 );
    }else{
      if( stype == "exp" )
	return 0.0;
      else if( stype == "logadd" )
	return LOG_ZERO;
      else assert( false );
      return 0.0;
    }
  }

  void Centroid::initForwardMatrix(){
    fM.resize( lastAntiDiagonal + 1 );
    fD.resize( lastAntiDiagonal + 1 );
    fI.resize( lastAntiDiagonal + 1 );
    fP.resize( lastAntiDiagonal + 1 );
    for( size_t k=0; k < fM.size(); ++k ){
      fM[k].resize( xa.x[k].size() ); 
      fD[k].resize( xa.x[k].size() ); 
      fI[k].resize( xa.x[k].size() ); 
      fP[k].resize( xa.x[k].size() ); 
    }
//     if( stype=="exp" ) Z = 0.0;
//     else if( stype=="logadd" ) Z = LOG_ZERO;
//     else assert( false );
    if( stype == "exp" ) {
      fM[0][0] = 1;
      fD[0][0] = 0;
      fI[0][0] = 0;
      fP[0][0] = 0;
    }
    else if( stype == "logadd" ){
      fM[0][0] = 0;
      fD[0][0] = LOG_ZERO;
      fI[0][0] = LOG_ZERO;
      fP[0][0] = LOG_ZERO;
    } else assert( false );
    Z = fM[0][0];
  }
  
  void Centroid::initBackwardMatrix(){
    bM.resize( fM.size() );
    bD.resize( fM.size() );
    bI.resize( fM.size() );
    bP.resize( fM.size() );
    pp.resize( fM.size() );
    for( size_t k = 0; k < fM.size(); ++k ){
      double d1 = 1.0;
      if( stype == "logadd" ) d1 = 0.0;
      double d0 = 0.0;
      if( stype == "logadd" ) d0 = LOG_ZERO;
      bM[ k ].resize( fM[ k ].size(), d1 ); 
      bD[ k ].resize( fM[ k ].size(), d0 ); 
      bI[ k ].resize( fM[ k ].size(), d0 ); 
      bP[ k ].resize( fM[ k ].size(), d0 ); 
      pp[ k ].resize( fM[ k ].size() );
    }
    //bM[ 0 ][ 0 ] = 0;
  }

  void Centroid::initDecodingMatrix(){
    X.resize( fM.size() );
    for( size_t k = 0; k < fM.size(); ++k ){
      double d1 = 1.0;
      if( stype == "logadd" ) d1 = 0.0;
      double d0 = 0.0;
      if( stype == "logadd" ) d0 = LOG_ZERO;
      X[ k ].resize( fM[ k ].size(), d1 ); 
    }
  }

  void Centroid::updateScore( double score, size_t antiDiagonal, size_t cur ){
    if( bestScore < score ){
      bestScore = score;
      bestAntiDiagonal = antiDiagonal;
      bestPos1 = cur;
    }
  }

  double Centroid::forward( const uchar* seq1, const uchar* seq2,
			    size_t start1, size_t start2, XdropAligner::direction dir,
			    const int sm[64][64], 
			    const GeneralizedAffineGapCosts& gap ){

    const int seqIncrement = (dir == XdropAligner::FORWARD) ? 1 : -1;

    initForwardMatrix();

    const int E = gap.extend;
    const int F = gap.first;
    const int P = gap.extendPair;
    const int Q = gap.firstPair;
    const double tE = E / T;
    const double tF = F / T;
    const double tP = P / T;
    const double tQ = Q / T;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );

    for( size_t k = 1; k <= lastAntiDiagonal; ++k ){  // loop over antidiagonals
      const size_t k1 = k - 1;
      const size_t k2 = k - 2;  // might wrap around
      const size_t off1 = xa.offsets[ k1 ];
      const size_t end1 = off1 + xa.x[ k1 ].size();
      const size_t off0 = xa.offsets[ k ];
      const size_t end0 = off0 + xa.x[ k ].size();
      const size_t loopBeg = off0 + ( off0 == off1 );
      const size_t loopEnd = end0 - ( end0 > end1 );
      const int* x0 = &xa.x[ k ][ 0 ];
      double* fM0 = &fM[ k ][ 0 ];
      double* fD0 = &fD[ k ][ 0 ];
      double* fI0 = &fI[ k ][ 0 ];
      double* fP0 = &fP[ k ][ 0 ];

      if( off0 == off1 ){  // do first cell on boundary
	int score = 0;
	double eS = 0;
	if ( off0 > 0 && k - off0 > 0 ) {
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, off0 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - off0 );
	  assert ( *s1 < 64 && *s2 < 64 );
	  score = sm[ *s1 ][ *s2 ];
	  eS = EXP( score / T );
	}
	*x0++;
	const double fM1 = fM[ k1 ].front();
	const double fD1 = fD[ k1 ].front();
	const double fI1 = fI[ k1 ].front();
	const double fP1 = fP[ k1 ].front();
	const double fM2 = diag( fM, k, off0 );
	const double fD2 = diag( fD, k, off0 );
	const double fI2 = diag( fI, k, off0 );
	const double fP2 = diag( fP, k, off0 );

	if( stype == "exp" ){
	  *fM0 = ( fM2 + fD2 + fI2 + fP2 ) * eS;
	  *fD0 = 0;
	  *fI0 = ( fM1 + fD1 ) * eF + ( fI1 + fP1 ) * eE;
	  *fP0 = ( fM2 ) * eQ + ( fP2 ) * eP;
	  Z += *fM0;
	  assert ( *fM0 >= 0 && *fD0 >= 0 && *fI0 >= 0 && *fP0 >=0);
	} else if( stype == "logadd" ) {
	  *fM0 = LOG_ADD( fM2 + score / T, fD2 + score / T, fI2 + score / T, fP2 + score / T);
	  *fD0 = LOG_ZERO;
	  *fI0 = LOG_ADD( fM1 - tF, fD1 - tF, fI1 - tE, fP1 - tE );
	  *fP0 = LOG_ADD( fM2 - Q / T, fP2 - P / T );
	  Z = LOG_ADD( Z, *fM0 );
	} else assert( false );
	fM0++; fD0++; fI0++; fP0++;
      }

      if( loopBeg < loopEnd ){
	assert( k > 1 );
	const int* const x0end = x0 + (loopEnd - loopBeg);
	const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, loopBeg );
	const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - loopBeg );
	assert ( *s1 < 64 && *s2 < 64 );
	const size_t horiBeg = loopBeg - 1 - off1;
	const size_t diagBeg = loopBeg - 1 - xa.offsets[ k2 ];
	const double* fM1 = &fM[ k1 ][ horiBeg ];
	const double* fD1 = &fD[ k1 ][ horiBeg ];
	const double* fI1 = &fI[ k1 ][ horiBeg ];
	const double* fP1 = &fP[ k1 ][ horiBeg ];
	const double* fM2 = &fM[ k2 ][ diagBeg ];
	const double* fD2 = &fD[ k2 ][ diagBeg ];
	const double* fI2 = &fI[ k2 ][ diagBeg ];
	const double* fP2 = &fP[ k2 ][ diagBeg ];
	do{
	  *x0++;
	  if( stype == "exp" ){
	    const double eS = EXP( sm[ *s1 ][ *s2 ] / T );
	    *fM0 = ( *fM2 + *fD2 + *fI2 + *fP2 ) * eS;
	    *fD0 = ( *fM1++ ) * eF + ( *fD1++ + *fP1++ ) * eE; fI1++;
	    *fI0 = ( *fM1 + *fD1 ) * eF + ( *fI1 + *fP1 ) * eE;
	    *fP0 = ( *fM2 ) * eQ + ( *fP2 ) * eP;
	    fM2++; fD2++; fI2++; fP2++;
	    Z += *fM0;
	    assert ( *fM0 >= 0 && *fD0 >= 0 && *fI0 >= 0 && *fP0 >=0);
	  }else if( stype == "logadd" ) {
	    const double tS = sm[ *s1 ][ *s2 ] / T;
	    *fM0 = LOG_ADD( *fM2, *fD2, *fI2, *fP2) + tS;
	    *fD0 = LOG_ADD( *fM1++ - tF, *fD1++ - tE, *fP1++ - tE); fI1++;
 	    *fI0 = LOG_ADD( *fM1 - tF, *fD1 - tF, *fI1 - tE, *fP1 - tE );
	    *fP0 = LOG_ADD( *fM2 - tQ, *fP2 - tP );
	    fM2++; fD2++; fI2++; fP2++;
	    Z = LOG_ADD( Z, *fM0 );
	  }else assert( false );
	  fM0++; fD0++; fI0++; fP0++;

	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( x0 != x0end );
      }
     
      if( end0 > end1 ){  // do last cell on boundary
	assert ( end0 == end1 + 1 );
	const double fM1 = fM[ k1 ].back();
	const double fD1 = fD[ k1 ].back();
	const double fP1 = fP[ k1 ].back();
	const double fM2 = diag( fM, k, end0 - 1 );
	const double fD2 = diag( fD, k, end0 - 1 );
	const double fI2 = diag( fI, k, end0 - 1 );
	const double fP2 = diag( fP, k, end0 - 1 );

	int score = 0;
	double eS = 0;
	if ( end0 > 1 && k - end0 + 1 > 0 ) {
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, end0 - 1 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - end0 + 1 );
	  assert ( *s1 < 64 && *s2 < 64 );
	  score = sm[ *s1 ][ *s2 ];
	  eS = EXP ( score / T );
	}
	*x0++;
	if( stype == "exp" ){
	  *fM0 = ( fM2 + fD2 + fI2 + fP2 ) * eS;
	  *fD0 = fM1 * eF + ( fD1 + fP1 ) * eE;
	  *fI0 = 0;
	  *fP0 = fM2 * eQ + fP2 * eP;
	  Z += *fM0;
	  assert ( *fM0 >= 0 && *fD0 >= 0 && *fI0 >= 0 && *fP0 >=0 );
	}else if( stype == "logadd" ){
	  *fM0 = LOG_ADD( fM2 + score / T, fD2 + score / T, fI2 + score / T, fP2 + score / T);
	  *fD0 = LOG_ADD( fM1 - tF, fD1 - tE, fP1 - tE);
	  *fI0 = LOG_ZERO;
	  *fP0 = LOG_ADD( fM2 - tQ, fP2 - tP);
	  Z = LOG_ADD( Z, *fM0 );
	} else assert( false );
	fM0++; fD0++; fI0++; fP0++;
      }
    }
    if( stype=="exp") return log(Z);
    return Z;
  }

  // added by M. Hamada
  // compute posterior probabilities while executing backward algorithm 
  // posterior probabilities are stored in pp
  double Centroid::backward( const uchar* seq1, const uchar* seq2,
			     size_t start1, size_t start2, XdropAligner::direction dir,
			     const int sm[64][64], 
			     const GeneralizedAffineGapCosts& gap ){

    const int seqIncrement = (dir == XdropAligner::FORWARD) ? 1 : -1;

    initBackwardMatrix();

    const int E = gap.extend;
    const int F = gap.first;
    const int P = gap.extendPair;
    const int Q = gap.firstPair;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );
    const double tE = E / T;
    const double tF = F / T;
    const double tP = P / T;
    const double tQ = Q / T;
    for( size_t k = lastAntiDiagonal; k > 0; --k ){  // loop over antidiagonals
      const size_t k1 = k - 1;  
      const size_t k2 = k - 2;  // might wrap around
      const size_t off1 = xa.offsets[ k1 ];
      const size_t end1 = off1 + xa.x[ k1 ].size();
      const size_t off0 = xa.offsets[ k ];
      const size_t end0 = off0 + xa.x[ k ].size();
      const size_t loopBeg = off0 + ( off0 == off1 );
      const size_t loopEnd = end0 - ( end0 > end1 );

      const int* x0 = &xa.x[ k ][ 0 ];
      const double* bM0 = &bM[ k ][ 0 ];
      const double* bD0 = &bD[ k ][ 0 ];
      const double* bI0 = &bI[ k ][ 0 ];
      const double* bP0 = &bP[ k ][ 0 ];
      double* pp0 = &pp[ k ][ 0 ];

      const double* fM0 = &fM[ k ][ 0 ];

      if( off0 == off1 ){  // do first cell on boundary
	double* bM1 = &bM[ k1 ][ 0 ];
	double* bD1 = &bD[ k1 ][ 0 ];
	double* bI1 = &bI[ k1 ][ 0 ];
	double* bP1 = &bP[ k1 ][ 0 ];
	if( stype == "exp" ){
	  *bM1 += *bI0 * eF;
	  *bD1 += *bI0 * eF;
	  *bI1 += *bI0 * eE;
	  *bP1 += *bI0 * eE;
	}
	else if( stype == "logadd" ) {
	  *bM1 = LOG_ADD( *bM1, *bI0 - tF );
	  *bD1 = LOG_ADD( *bD1, *bI0 - tF );
	  *bI1 = LOG_ADD( *bI1, *bI0 - tE );
	  *bP1 = LOG_ADD( *bP1, *bI0 - tE );
	} else assert( false );
	if( k2 < k && xa.offsets[ k2 ] + 1 <= off0 ) { // there exists diagonal values
	  double eS = 1;
	  int score = 0;
	  //if ( off0 > 0 && k - off0 > 0 ) {
	    const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, off0 );
	    const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - off0 );
	    assert ( *s1 < 64 && *s2 < 64 );
	    score = sm[ *s1 ][ *s2 ];
	    eS = EXP( score / T );
	    //}
	  const size_t dig = off0 - 1 - xa.offsets[ k2 ]; //
	  double* bM2 = &bM[ k2 ][ dig ];
	  double* bD2 = &bD[ k2 ][ dig ];
	  double* bI2 = &bI[ k2 ][ dig ];
	  double* bP2 = &bP[ k2 ][ dig ];
	  if( stype == "exp" ){
	    *bM2 += *bM0 * eS + *bP0 * eQ;
	    *bD2 += *bM0 * eS;
	    *bI2 += *bM0 * eS;
	    *bP2 += *bM0 * eS + *bP0 * eP;
	  }else if( stype == "logadd" ){
	    double tS = score / T;
	    *bM2 = LOG_ADD( *bM2, *bM0 + tS, *bP0 - tQ );
	    *bD2 = LOG_ADD( *bD2, *bM0 + tS );
	    *bI2 = LOG_ADD( *bI2, *bM0 + tS );
	    *bP2 = LOG_ADD( *bP2, *bM0 + tS, *bP0 - tP );
	  }else assert (false);
	}
	double prob = 0;
	if( stype == "exp" ) prob = *fM0 * *bM0 / Z; 
	else if( stype == "logadd" ) prob = std::exp( *fM0 + *bM0 - Z );
	assert( 0 <= prob && prob <= 1 );
	*pp0 = prob;
	// iteration
	bM0++; bD0++; bI0++; bP0++; 
	fM0++; 
	pp0++;
	x0++;
      }

      if( loopBeg < loopEnd ){
	assert( k > 1 );
	const int* const x0end = x0 + (loopEnd - loopBeg);
	const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, loopBeg );
	const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - loopBeg );
	assert ( *s1 < 64 && *s2 < 64 );
	const size_t horiBeg = loopBeg - 1 - off1;
	const size_t diagBeg = loopBeg - 1 - xa.offsets[ k2 ];
	double* bM1 = &bM[ k1 ][ horiBeg ];
	double* bD1 = &bD[ k1 ][ horiBeg ];
	double* bI1 = &bI[ k1 ][ horiBeg ];
	double* bP1 = &bP[ k1 ][ horiBeg ];
	double* bM2 = &bM[ k2 ][ diagBeg ];
	double* bD2 = &bD[ k2 ][ diagBeg ];
	double* bI2 = &bI[ k2 ][ diagBeg ];
	double* bP2 = &bP[ k2 ][ diagBeg ];

	do{
	  if( stype == "exp" ) {
	    const double eS = EXP( sm[ *s1 ][ *s2 ] / T );
	    *bM2 += *bM0 * eS + *bP0 * eQ;
	    *bD2 += *bM0 * eS;
	    *bI2 += *bM0 * eS;
	    *bP2 += *bM0 * eS + *bP0 * eP;
	    *bM1++ += *bD0 * eF;
	    *bD1++ += *bD0 * eE;
	    bI1++;
	    *bP1++ += *bD0 * eE;
	    *bM1 += *bI0 * eF;
	    *bD1 += *bI0 * eF;
	    *bI1 += *bI0 * eE;
	    *bP1 += *bI0 * eE;
	  } else if( stype == "logadd" ){
	    const double tS = ( sm[ *s1 ][ *s2 ] / T );
	    *bM2 = LOG_ADD( *bM2, *bM0 + tS, *bP0 - tQ );
	    *bD2 = LOG_ADD( *bD2, *bM0 + tS );
	    *bI2 = LOG_ADD( *bI2, *bM0 + tS );
	    *bP2 = LOG_ADD( *bP2, *bM0 + tS, *bP0 - tP );
	    *bM1 = LOG_ADD( *bM1, *bD0 - tF );
	    *bD1 = LOG_ADD( *bD1, *bD0 - tE );
	    *bP1 = LOG_ADD( *bP1, *bD0 - tE );
	    bM1++; bD1++; bI1++; bP1++;
	    *bM1 = LOG_ADD( *bM1, *bI0 - tF );
	    *bD1 = LOG_ADD( *bD1, *bI0 - tF );
	    *bI1 = LOG_ADD( *bI1, *bI0 - tE );
	    *bP1 = LOG_ADD( *bP1, *bI0 - tE );
	  }
	  double prob = 0;
	  if( stype == "exp" ) prob = *fM0 * *bM0 / Z; 
	  else if( stype == "logadd" ) prob = std::exp( *fM0 + *bM0 - Z );
	  assert( 0 <= prob && prob <= 1 );
	  *pp0 = prob;
	  // iteration
	  bM2++; bD2++; bI2++; bP2++;
	  bM0++; bD0++; bI0++; bP0++;
	  fM0++; 
	  pp0++;
	  x0++;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( x0 != x0end );
      }

      if( end0 > end1 ){  // do last cell on boundary
	assert ( end0 == end1 + 1 );
	double* bM1 = &bM[ k1 ][ bM[ k1 ].size() - 1 ];
	double* bD1 = &bD[ k1 ][ bD[ k1 ].size() - 1 ];
	double* bP1 = &bP[ k1 ][ bP[ k1 ].size() - 1 ];
	if( stype == "exp" ){
	  *bM1 += *bD0 * eF;
	  *bD1 += *bD0 * eE;
	  *bP1 += *bD0 * eE;
	}else if( stype == "logadd" ){
	  *bM1 = LOG_ADD( *bM1, *bD0 - tF );
	  *bD1 = LOG_ADD( *bD1, *bD0 - tE );
	  *bP1 = LOG_ADD( *bP1, *bD0 - tE );	  
	}
	if( k2 < k ) {
	  const size_t off2 = xa.offsets[ k2 ];
	  const size_t end2 = off2 + xa.x[ k2 ].size();
	  if( end2 + 1 >= end0 ) { // there exists diagonal
	    const size_t dig = end0 - 2 - off2; // diagonal 
	    double* bM2 = &bM[ k2 ][ dig ];
	    double* bD2 = &bD[ k2 ][ dig ];
	    double* bI2 = &bI[ k2 ][ dig ];
	    double* bP2 = &bP[ k2 ][ dig ];
	    int score = 0;
	    double eS = 1;
	    double tS = 0;
	    if ( end0 > 1 && k - end0 + 1 > 0 ) { 
	      const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, end0 - 1 );
	      const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - end0 + 1);
	      assert ( *s1 < 64 && *s2 < 64 );
	      score = sm[ *s1 ][ *s2 ];
	      eS = EXP ( score / T );
	      tS = score / T;
	    }
	    if( stype == "exp" ) {
	      const double tmp1 = *bM0 * eS;
	      *bM2 += ( tmp1 + *bP0 * eQ );
	      *bD2 += tmp1;
	      *bI2 += tmp1;
	      *bP2 += tmp1 + *bP0 * eP;
	    } else if( stype == "logadd" ) {
	      const double tmp1 = *bM0 + tS;
	      *bM2 = LOG_ADD( *bM2, tmp1, *bP0 - tQ );
	      *bD2 = LOG_ADD( *bD2, tmp1 );
	      *bI2 = LOG_ADD( *bI2, tmp1 );
	      *bP2 = LOG_ADD( *bP2, tmp1,  *bP0 - tP );	      
	    } else assert (false);
	  }
	} // if
	double prob = 0;
	if( stype == "exp" ) prob = *fM0 * *bM0 / Z; 
	else if( stype == "logadd" ) prob = std::exp( *fM0 + *bM0 - Z );
	assert( 0 <= prob && prob <= 1 );
	*pp0 = prob;
	bM0++; bD0++; bI0++; bP0++;
	fM0++; 
	pp0++;
	x0++;
      }
    }
    if( stype=="exp") return log(bM[0][0]);
    return bM[0][0];
  }

  double Centroid::dp( double gamma ){

    initDecodingMatrix();

    for( size_t k = 1; k <= lastAntiDiagonal; ++k ){  // loop over antidiagonals
      const size_t k1 = k - 1;
      const size_t k2 = k - 2;  // might wrap around
      const size_t off1 = xa.offsets[ k1 ];
      const size_t end1 = off1 + xa.x[ k1 ].size();
      const size_t off0 = xa.offsets[ k ];
      const size_t end0 = off0 + xa.x[ k ].size();
      const size_t loopBeg = off0 + ( off0 == off1 );
      const size_t loopEnd = end0 - ( end0 > end1 );
      const double* p0 = &pp[ k ][ 0 ]; // 

      double* X0 = &X[ k ][ 0 ];
      const double* P0 = &pp[ k ][ 0 ];
      size_t cur = off0;

      if( off0 == off1 ){  // do first cell on boundary
	p0++;
	const double X1 = X[ k1 ].front();
	const double X2 = xa.diag( X, k, off0 );
	const double s = ( gamma + 1 ) * ( *P0++ ) - 1;
	const double score = std::max( X1, X2 + s );
	assert ( score >= 0 );
	updateScore ( score, k, cur );
	*X0++ = score;
	cur++;
      }

      if( loopBeg < loopEnd ){
	assert( k > 1 );
	const double* const p0end = p0 + (loopEnd - loopBeg);
	const size_t horiBeg = loopBeg - 1 - off1;
	const size_t diagBeg = loopBeg - 1 - xa.offsets[ k2 ];
	const double* X1 = &X[ k1 ][ horiBeg ];
	const double* X2 = &X[ k2 ][ diagBeg ];
	do{
	  p0++;
	  const double s = ( gamma + 1 ) * ( *P0++ ) - 1;
	  const double oldX1 = *X1++;  // Added by MCF
	  const double score = std::max( std::max( oldX1, *X1 ), *X2++ + s );
	  assert ( score >= 0 );
	  updateScore ( score, k, cur );
	  *X0++ = score;
	  cur++;
	}while( p0 != p0end );
      }
 
      if( end0 > end1 ){  // do last cell on boundary
	assert ( end0 == end1 + 1 );
	const double X1 = X[ k1 ].back();
	const double X2 = xa.diag( X, k, end0 - 1 );
	const double s = ( gamma + 1 ) * ( *P0++ ) - 1;
	const double score = std::max( X1, X2 + s );
	assert ( score >= 0 );
	updateScore ( score, k, cur );
	*X0++ = score;
	cur++;
      }
    }

    return bestScore;
  }

  void Centroid::traceback( std::vector< SegmentPair >& chunks,
			    double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const int m =
	maxIndex3( xa.diag( X, k, i ) + ( gamma + 1 ) * xa.cell( pp, k, i ) - 1,
		   xa.hori( X, k, i ),
		   xa.vert( X, k, i ) );
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

  // Added by MCF:
  void Centroid::chunkProbabilities( std::vector<double>& probs,
				     const std::vector<SegmentPair>& chunks ){
    for( CI(SegmentPair) i = chunks.begin(); i < chunks.end(); ++i ){
      size_t seq1pos = i->end1();
      size_t seq2pos = i->end2();

      for( size_t j = 0; j < i->size; ++j ){
	probs.push_back( xa.cell( pp, seq1pos + seq2pos, seq1pos ) );
	--seq1pos;
	--seq2pos;
      }
    }
  }

  //
  // for debug
  //
  void Centroid::print_forward_matrix( std::ostream& os ) const{
    os << "fM" << std::endl;
    print_matrix ( os, fM, xa.offsets );
    os << "fD" << std::endl;
    print_matrix ( os, fD, xa.offsets );
    os << "fI" << std::endl;
    print_matrix ( os, fI, xa.offsets );
    os << "fP" << std::endl;
    print_matrix ( os, fP, xa.offsets );
  }

  void Centroid::print_backward_matrix( std::ostream& os ) const{
    os << "bM" << std::endl;
    print_matrix ( os, bM, xa.offsets );
    os << "bD" << std::endl;
    print_matrix ( os, bD, xa.offsets );
    os << "bI" << std::endl;
    print_matrix ( os, bI, xa.offsets );
    os << "bP" << std::endl;
    print_matrix ( os, bP, xa.offsets );
  }

  void Centroid::print_prob_matrix( std::ostream& os ) const{
    for( size_t k = 0; k < pp.size(); ++k ){
      const size_t off = xa.offsets[ k ];
      os << "(" << k << ")" << " [" << off << "] ";
      for( size_t p = 0; p < pp[ k ].size(); ++p ){
	if( ! (0 <= pp[ k ][ p ] && pp[ k ][ p ] <= 1) ) {
	  os << "pp[" << k << "][" << p << "]=" <<  pp[k][p] << std::endl;
	}
	assert (0 <= pp[ k ][ p ] && pp[ k ][ p ] <= 1);
	os << pp[ k ][ p ] << " ";
      }
      os << "\n";
    }    
  }

}  // end namespace cbrc

