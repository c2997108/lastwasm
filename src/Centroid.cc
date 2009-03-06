// Copyright 2008, 2009 Michiaki Hamada

#include "Centroid.hh"
#include <algorithm>  // max_element
#include <cassert>
#include <cmath> // for exp
#include <string>
#include <cfloat>   // for DBL_MAX

#define CI(type) std::vector<type>::const_iterator  // added by MCF

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
    MM = d0; MD = d0; MP = d0; MI = d0; MQ = d0;
    DD = d0; DM = d0; DI = d0;
    PP = d0; PM = d0; PD = d0; PI = d0;
    II = d0; IM = d0;
    SM = d0; SD = d0; SP = d0; SI = d0; SQ = d0;

    for (int n=0; n<MAT; n++)
      for (int m=0; m<MAT; m++) emit[n][m] = d0;
  }
  
  std::ostream& ExpectedCount::write (std::ostream& os, double Z) const
  {
    for (int n=0; n<MAT; ++n) {
      for (int m=0; m<MAT; ++m) {
	double prob = emit[n][m] / Z;
	if (prob > 0) 
	  os << "emit[" << n << "][" << m << "]=" << emit[n][m] / Z << std::endl;
      }
    }
    os << "M->M=" << MM / Z << std::endl;
    os << "M->D=" << MD / Z << std::endl;
    os << "M->P=" << MP / Z << std::endl;
    os << "M->I=" << MI / Z << std::endl;
    os << "M->Q=" << MQ / Z << std::endl;

    os << "D->D=" << DD / Z << std::endl;
    os << "D->M=" << DM / Z << std::endl;
    os << "D->I=" << DI / Z << std::endl;

    os << "P->P=" << PP / Z << std::endl;
    os << "P->M=" << PM / Z << std::endl;
    os << "P->D=" << PD / Z << std::endl;
    os << "P->I=" << PI / Z << std::endl;
  
    os << "I->I=" << II / Z << std::endl;
    os << "I->M=" << IM / Z << std::endl;

    os << "S->Q=" << SQ / Z << std::endl;

    //os << ( MQ + SQ ) / Z << std::endl; // must be equal to 1
    assert ( std::abs ((MQ + SQ) / Z - 1) < 1e-10);
    return os;
  }

  // get DP matrix value "diagonal from" the given position
  // return 0 if there does not exist diagonal value
  double Centroid::diag( const dmatrix_t& matrix,
			 size_t antiDiagonal, size_t seq1pos ) const{
    if( antiDiagonal > 1 &&
	seq1pos > xa.fillBeg( antiDiagonal-2 ) &&
	seq1pos <= xa.fillEnd( antiDiagonal-2 ) ){
      return xa.cell( matrix, antiDiagonal-2, seq1pos-1 );
    }else{
      return 0.0;
    }
  }

  Centroid::Centroid( const XdropAligner& xa_, const int sm[MAT][MAT], double T_ ) 
    : xa( xa_ ), T( T_ ), lastAntiDiagonal ( xa_.offsets.size () - 1 ), bestScore ( 0 ),
      bestAntiDiagonal (0), bestPos1 (0) {
    for ( int n=0; n<MAT; ++n )
      for ( int m=0; m<MAT; ++m ) {
	match_score[n][m] = EXP ( sm[ n ][ m ] / T );
      }
  }


  void Centroid::initForwardMatrix(){
    fM.resize( lastAntiDiagonal + 1 );
    fD.resize( lastAntiDiagonal + 1 );
    fI.resize( lastAntiDiagonal + 1 );
    fP.resize( lastAntiDiagonal + 1 );
    scale.resize ( lastAntiDiagonal + 1, 1.0 ); // scaling

    for( size_t k=0; k < fM.size(); ++k ){
      fM[k].resize( xa.x[k].size() ); 
      fD[k].resize( xa.x[k].size() ); 
      fI[k].resize( xa.x[k].size() ); 
      fP[k].resize( xa.x[k].size() ); 
    }
    fM[0][0] = 1;
    fD[0][0] = 0;
    fI[0][0] = 0;
    fP[0][0] = 0;

    Z = fM[0][0];
  }
  
  void Centroid::initBackwardMatrix(){
    bM.resize( fM.size() );
    bD.resize( fM.size() );
    bI.resize( fM.size() );
    bP.resize( fM.size() );
    pp.resize( fM.size() );
    for( size_t k = fM.size() - 1; k != (size_t) -1; --k ){
      double d1 = 1.0;
      if ( k != fM.size() - 1 ) d1 = bM[ k + 1 ][0] / (scale[k]);
      else d1 = 1.0 / scale[k];

      double d0 = 0.0;
      bM[ k ].assign( fM[ k ].size(), d1 );
      bD[ k ].assign( fM[ k ].size(), d0 );
      bI[ k ].assign( fM[ k ].size(), d0 );
      bP[ k ].assign( fM[ k ].size(), d0 );
      pp[ k ].resize( fM[ k ].size() );
    }
  }

  void Centroid::initDecodingMatrix(){
    X.resize( fM.size() );
    for( size_t k = 0; k < fM.size(); ++k ){
      X[ k ].resize( fM[ k ].size(), 0.0 );
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
			    const GeneralizedAffineGapCosts& gap ){

    const int seqIncrement = (dir == XdropAligner::FORWARD) ? 1 : -1;

    initForwardMatrix();

    const int E = gap.extend;
    const int F = gap.first;
    const int P = gap.extendPair;
    const int Q = gap.firstPair;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );

    for( size_t k = 1; k <= lastAntiDiagonal; ++k ){  // loop over antidiagonals
      double sum_f = 0.0; // sum of forward values
      const size_t k1 = k - 1;
      const size_t k2 = k - 2;  // might wrap around
      const size_t off1 = xa.offsets[ k1 ];
      const size_t end1 = off1 + xa.x[ k1 ].size();
      const size_t off0 = xa.offsets[ k ];
      const size_t end0 = off0 + xa.x[ k ].size();
      const size_t loopBeg = off0 + ( off0 == off1 );
      const size_t loopEnd = end0 - ( end0 > end1 );

      const double scale12 = ( k > 1 ) ? 1.0 / ( scale[k1] * scale[k2] ) : 1.0; // scaling factor
      const double scale1  = 1.0 / scale[k1];

      const double seF = eF * scale1;
      const double seE = eE * scale1;
      const double seQ = eQ * scale12;
      const double seP = eP * scale12;

      double* fM0 = &fM[ k ][ 0 ];
      double* fD0 = &fD[ k ][ 0 ];
      double* fI0 = &fI[ k ][ 0 ];
      double* fP0 = &fP[ k ][ 0 ];

      if( off0 == off1 ){  // do first cell on boundary
	double eS = 0;
	if ( off0 > 0 && k - off0 > 0 ) {
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, off0 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - off0 );
	  assert ( *s1 < MAT && *s2 < MAT );
	  eS = match_score[ *s1 ][ *s2 ];
	}
	const double fM1 = fM[ k1 ].front();
	const double fD1 = fD[ k1 ].front();
	const double fI1 = fI[ k1 ].front();
	const double fP1 = fP[ k1 ].front();
	const double fM2 = diag( fM, k, off0 );
	const double fD2 = diag( fD, k, off0 );
	const double fI2 = diag( fI, k, off0 );
	const double fP2 = diag( fP, k, off0 );

	if ( k > 1 )
	  *fM0 = ( fM2 + fD2 + fI2 + fP2 ) * eS * scale12;
	*fD0 = 0;
	*fI0 = ( ( fM1 + fD1 ) * eF + ( fI1 + fP1 ) * eE ) * scale1;
	if ( k > 1 )
	  *fP0 = ( fM2 * eQ +  fP2 * eP ) * scale12;
	sum_f += *fM0;
	assert ( *fM0 < DINF && *fD0 < DINF && *fI0 < DINF && *fP0 < DINF );

	fM0++; fD0++; fI0++; fP0++;
      }

      if( loopBeg < loopEnd ){
	assert( k > 1 );
	const double* const fM0end = fM0 + (loopEnd - loopBeg);
	const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, loopBeg );
	const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - loopBeg );
	assert ( *s1 < MAT && *s2 < MAT );
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

	double xM1 = *fM1, xD1 = *fD1, xI1 = *fI1, xP1 = *fP1;
	do{ // start: inner most loop
	  const double S = match_score[ *s1 ][ *s2 ] * scale12;
	  const double xM2 = *fM2, xD2 = *fD2, xI2 = *fI2, xP2 = *fP2;
	  *fM0 = ( xM2 + xD2 + xI2 + xP2 ) * S;
	  *fD0 = ( xM1 ) * seF + ( xD1 + xP1 ) * seE;
	  xM1 = *++fM1; xD1 = *++fD1; xI1 = *++fI1; xP1 = *++fP1;
	  *fI0 = ( xM1 + xD1 ) * seF + ( xI1 + xP1 ) * seE;
	  *fP0 = xM2 * seQ + xP2 * seP;
	  fM2++; fD2++; fI2++; fP2++;
	  sum_f += *fM0;
	  //assert ( *fM0 < DINF && *fD0 < DINF && *fI0 < DINF && *fP0 < DINF );
	  fM0++; fD0++; fI0++; fP0++;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( fM0 != fM0end ); // end: inner most loop
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

	double eS = 0;
	if ( end0 > 1 && k - end0 + 1 > 0 ) {
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, end0 - 1 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - end0 + 1 );
	  assert ( *s1 < MAT && *s2 < MAT );
	  eS = match_score[ *s1 ][ *s2 ];
	}
	if ( k > 1 ) *fM0 = ( fM2 + fD2 + fI2 + fP2 ) * eS * scale12;
	*fD0 = ( fM1 * eF + ( fD1 + fP1 ) * eE ) * scale1;
	*fI0 = 0;
	if ( k > 1 ) *fP0 = ( fM2 * eQ + fP2 * eP ) * scale12;
	sum_f += *fM0;
	assert ( *fM0 < DINF && *fD0 < DINF && *fI0 < DINF && *fP0 < DINF );
	fM0++; fD0++; fI0++; fP0++;
      }

      Z += sum_f;
      scale[k] = sum_f + 1.0;  // seems ugly
      Z /= scale[k]; // scaling
    } // k
    return log(Z);
  }

  // added by M. Hamada
  // compute posterior probabilities while executing backward algorithm 
  // posterior probabilities are stored in pp
  double Centroid::backward( const uchar* seq1, const uchar* seq2,
			     size_t start1, size_t start2, XdropAligner::direction dir,
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

    for( size_t k = lastAntiDiagonal; k > 0; --k ){  // loop over antidiagonals
      const size_t k1 = k - 1;  
      const size_t k2 = k - 2;  // might wrap around
      const size_t off2 = (k2 < k) ? xa.offsets[ k2 ] : 0;
      const size_t end2 = (k2 < k) ? off2 + xa.x[ k2 ].size () : 0;
      const size_t off1 = xa.offsets[ k1 ];
      const size_t end1 = off1 + xa.x[ k1 ].size();
      const size_t off0 = xa.offsets[ k ];
      const size_t end0 = off0 + xa.x[ k ].size();
      const size_t loopBeg = off0 + ( off0 == off1 );
      const size_t loopEnd = end0 - ( end0 > end1 );

      const double scale12 = ( k > 1 ) ? 1.0 / ( scale[k1] * scale[k2] ) : 1.0; // scaling factor
      const double scale1  = 1.0 / scale[k1];

      const double seF = eF * scale1;
      const double seE = eE * scale1;
      const double seQ = eQ * scale12;
      const double seP = eP * scale12;

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
	*bM1 += *bI0 * eF * scale1;
	*bD1 += *bI0 * eF * scale1;
	*bI1 += *bI0 * eE * scale1;
	*bP1 += *bI0 * eE * scale1;

	if( k2 < k && xa.offsets[ k2 ] + 1 <= off0 ) { // there exists diagonal values
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, off0 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - off0 );
	  assert ( *s1 < MAT && *s2 < MAT );
	  double eS = match_score[ *s1 ][ *s2 ];

	  const size_t dig = off0 - 1 - xa.offsets[ k2 ];
	  double* bM2 = &bM[ k2 ][ dig ];
	  double* bD2 = &bD[ k2 ][ dig ];
	  double* bI2 = &bI[ k2 ][ dig ];
	  double* bP2 = &bP[ k2 ][ dig ];
	  *bM2 += ( *bM0 * eS + *bP0 * eQ ) * scale12;
	  *bD2 += *bM0 * eS * scale12;
	  *bI2 += *bM0 * eS * scale12;
	  *bP2 += ( *bM0 * eS + *bP0 * eP ) * scale12;
	}
	double prob = *fM0 * *bM0 / Z; 
	//assert( 0 <= prob && prob <= 1 );
	*pp0 = prob;
	// iteration
	bM0++; bD0++; bI0++; bP0++; 
	fM0++; 
	pp0++;
      }

      if( loopBeg < loopEnd ){
	assert( k > 1 );
	const double* const fM0end = fM0 + (loopEnd - loopBeg);
	const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, loopBeg );
	const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - loopBeg );
	assert ( *s1 < MAT && *s2 < MAT );
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

	do{ // inner most loop
	  const double S = match_score[ *s1 ][ *s2 ];
	  const double tmp1 = *bM0 * S * scale12;
	  const double tmp2 = *bP0;
	  *bM2 += tmp1 + tmp2 * seQ; 
	  *bD2 += tmp1; 
	  *bI2 += tmp1; 
	  *bP2 += tmp1 + tmp2 * seP; 
	  const double tmp3 = *bD0;
	  *bM1++ += tmp3 * seF; 
	  *bD1++ += tmp3 * seE; 
	  bI1++;
	  *bP1++ += tmp3 * seE; 
	  const double tmp4 = *bI0;
	  const double tmp5 = tmp4 * seF; 
	  const double tmp6 = tmp4 * seE;
	  *bM1 += tmp5; 
	  *bD1 += tmp5; 
	  *bI1 += tmp6; 
	  *bP1 += tmp6; 
	  double prob = *fM0 * *bM0 / Z; 
	  //assert( 0 <= prob && prob <= 1 );
	  *pp0 = prob;
	  // iteration
	  bM2++; bD2++; bI2++; bP2++;
	  bM0++; bD0++; bI0++; bP0++;
	  fM0++; 
	  pp0++;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( fM0 != fM0end ); // inner most loop end;

      }

      if( end0 > end1 ){  // do last cell on boundary
	assert ( end0 == end1 + 1 );
	double* bM1 = &bM[ k1 ][ bM[ k1 ].size() - 1 ];
	double* bD1 = &bD[ k1 ][ bD[ k1 ].size() - 1 ];
	double* bP1 = &bP[ k1 ][ bP[ k1 ].size() - 1 ];
	*bM1 += *bD0 * eF * scale1;
	*bD1 += *bD0 * eE * scale1;
	*bP1 += *bD0 * eE * scale1;

	if( k2 < k && end2 + 1 >= end0 ) { // there exists diagonal
	  const size_t dig = end0 - 2 - off2; // diagonal 
	  double* bM2 = &bM[ k2 ][ dig ];
	  double* bD2 = &bD[ k2 ][ dig ];
	  double* bI2 = &bI[ k2 ][ dig ];
	  double* bP2 = &bP[ k2 ][ dig ];
	  //double eS = 1;
	  //if ( end0 > 1 && k - end0 + 1 > 0 ) { // comment out because not need
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, end0 - 1 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - end0 + 1);
	  assert ( *s1 < MAT && *s2 < MAT );
	  double eS = match_score[ *s1 ][ *s2 ];
	  //} 
	  const double tmp1 = *bM0 * eS;
	  *bM2 += ( tmp1 + *bP0 * eQ ) * scale12;
	  *bD2 += tmp1 * scale12;
	  *bI2 += tmp1 * scale12;
	  *bP2 += ( tmp1 + *bP0 * eP ) * scale12;
	}

	double prob = *fM0 * *bM0 / Z; 
	//assert( 0 <= prob && prob <= 1 );
	*pp0 = prob;
	bM0++; bD0++; bI0++; bP0++;
	fM0++; 
	pp0++;
      }
    }
    // This is a test for computeExpectedCounts (...)
    //ExpectedCount c;
    //computeExpectedCounts (seq1, seq2, start1, start2, dir, gap, c); 
    //c.write (std::cout, Z);
    //
    return log( bM[0][0] );
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

  
  void Centroid::computeExpectedCounts ( const uchar* seq1, const uchar* seq2,
					 size_t start1, size_t start2, XdropAligner::direction dir,
					 const GeneralizedAffineGapCosts& gap, 
					 ExpectedCount& c ) const{

    const int seqIncrement = (dir == XdropAligner::FORWARD) ? 1 : -1;

    const int E = gap.extend;
    const int F = gap.first;
    const int P = gap.extendPair;
    const int Q = gap.firstPair;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );

    c.SQ = 1; 

    for( size_t k = 1; k <= lastAntiDiagonal; ++k ){  // loop over antidiagonals
      const size_t k1 = k - 1;
      const size_t k2 = k - 2;  // might wrap around
      const size_t off1 = xa.offsets[ k1 ];
      const size_t end1 = off1 + xa.x[ k1 ].size();
      const size_t off0 = xa.offsets[ k ];
      const size_t end0 = off0 + xa.x[ k ].size();
      const size_t loopBeg = off0 + ( off0 == off1 );
      const size_t loopEnd = end0 - ( end0 > end1 );

      const double* fM0 = &fM[ k ][ 0 ];
      const double* fD0 = &fD[ k ][ 0 ];
      const double* fI0 = &fI[ k ][ 0 ];
      const double* fP0 = &fP[ k ][ 0 ];
      const double* bM0 = &bM[ k ][ 0 ];
      const double* bD0 = &bD[ k ][ 0 ];
      const double* bI0 = &bI[ k ][ 0 ];
      const double* bP0 = &bP[ k ][ 0 ];

      const double scale12 = ( k > 1 ) ? 1.0 / ( scale[k1] * scale[k2] ) : 1.0; 
      const double scale1  = 1.0 / scale[k1];
      const double scale0  = 1.0 / scale[k];

      if( off0 == off1 ){  // do first cell on boundary
	double S = 0;
	if ( off0 > 0 && k - off0 > 0 ) {
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, off0 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - off0 );
	  assert ( *s1 < MAT && *s2 < MAT );
	  S = match_score[ *s1 ][ *s2 ];
	  c.emit[ *s1 ][ *s2 ] += ( *fM0 * *bM0 );
	}
	const double fM1 = fM[ k1 ].front();
	const double fD1 = fD[ k1 ].front();
	const double fI1 = fI[ k1 ].front();
	const double fP1 = fP[ k1 ].front();
	const double fM2 = diag( fM, k, off0 );
	const double fD2 = diag( fD, k, off0 );
	const double fI2 = diag( fI, k, off0 );
	const double fP2 = diag( fP, k, off0 );

	c.MM += ( fM2 * S ) * *bM0 * scale12;
	c.PM += ( fP2 * S ) * *bM0 * scale12;
	c.DM += ( fD2 * S ) * *bM0 * scale12;
	c.IM += ( fI2 * S ) * *bM0 * scale12;
	c.MP += ( fM2 * eQ ) * *bP0 * scale12;
	c.PP += ( fP2 * eP ) * *bP0 * scale12;
	c.MQ += *fM0;

	c.MI += ( fM1 * eF ) * *bI0 * scale1;
	c.DI += ( fD1 * eF ) * *bI0 * scale1;
	c.PI += ( fP1 * eE ) * *bI0 * scale1;
	c.II += ( fI1 * eE ) * *bI0 * scale1;

	fM0++; fD0++; fI0++; fP0++;
	bM0++; bD0++; bI0++; bP0++;
      }

      if( loopBeg < loopEnd ){
	assert( k > 1 );
	const double* const fM0end = fM0 + (loopEnd - loopBeg);
	const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, loopBeg );
	const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - loopBeg );
	assert ( *s1 < MAT && *s2 < MAT );
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
	do{ // inner most loop
	  const double S = match_score[ *s1 ][ *s2 ];
	  c.emit[*s1][*s2] += ( *fM0 * *bM0 ) ;

	  const double tmp1 = S * *bM0 * scale12;
	  const double tmp2 = *bP0 * scale12;

	  c.MM += *fM2 * tmp1; 
	  c.PM += *fP2 * tmp1; 
	  c.DM += *fD2 * tmp1; 
	  c.IM += *fI2 * tmp1; 
	  c.MP += *fM2 * eQ * tmp2; 
	  c.PP += *fP2 * eP * tmp2;
	  c.MQ += *fM0;

	  const double tmp3 = *bD0 * scale1;
	  c.MD += ( *fM1 * eF ) * tmp3;
	  c.DD += ( *fD1 * eE ) * tmp3;
	  c.PD += ( *fP1 * eE ) * tmp3;

	  fM1++; fD1++; fP1++; fI1++;

	  const double tmp4 = *bI0 * scale1;
	  c.MI += ( *fM1 * eF )  * tmp4;
	  c.DI += ( *fD1 * eF )  * tmp4;
	  c.PI += ( *fP1 * eE )  * tmp4;
	  c.II += ( *fI1 * eE )  * tmp4;

	  fM2++; fD2++; fI2++; fP2++;
	  fM0++; fD0++; fI0++; fP0++;
	  bM0++; bD0++; bI0++; bP0++;
	    
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	} while( fM0 != fM0end );
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

	double S = 0;
	if ( end0 > 1 && k - end0 + 1 > 0 ) {
	  const uchar* s1 = XdropAligner::seqPtr( seq1, start1, dir, end0 - 1 );
	  const uchar* s2 = XdropAligner::seqPtr( seq2, start2, dir, k - end0 + 1 );
	  assert ( *s1 < MAT && *s2 < MAT );
	  S = match_score[ *s1 ][ *s2 ];
	  c.emit[ *s1 ][ *s2 ] += ( *fM0 * *bM0 ) ;
	}

	c.MM += ( fM2 * S ) * *bM0 * scale12;
	c.PM += ( fP2 * S ) * *bM0 * scale12;
	c.DM += ( fD2 * S ) * *bM0 * scale12;
	c.IM += ( fI2 * S ) * *bM0 * scale12;
	c.MP += ( fM2 * eQ ) * *bP0 * scale12;
	c.PP += ( fP2 * eP ) * *bP0 * scale12;
	c.MQ += *fM0;

	c.MD += ( fM1 * eF ) * *bD0 * scale1;
	c.DD += ( fD1 * eE ) * *bD0 * scale1;
	c.PD += ( fP1 * eE ) * *bD0 * scale1;

	fM0++; fD0++; fI0++; fP0++;
	bM0++; bD0++; bI0++; bP0++;
      }
      c.MQ *= scale0;
      c.SQ *= scale0;
    }
  }
}  // end namespace cbrc

