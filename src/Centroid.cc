// Copyright 2008, 2009, 2010, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#include "Centroid.hh"
#include "GappedXdropAlignerInl.hh"
#include <algorithm>
#include <cassert>
#include <cmath> // for exp
#include <cfloat>   // for DBL_MAX
#include <cstdlib>  // for abs
#include <iomanip> 

#define CI(type) std::vector<type>::const_iterator  // added by MCF

static const double DINF = DBL_MAX / 2;

namespace{
  double EXP ( double x ) {
    return std::exp (x);
  }
}


namespace cbrc{

  typedef double ScoreMatrixRowDouble[scoreMatrixRowSize];

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

  Centroid::Centroid( const GappedXdropAligner& xa_ ) 
    : xa( xa_ ), numAntidiagonals ( xa_.numAntidiagonals () ), bestScore ( 0 ), bestAntiDiagonal (0), bestPos1 (0)  {
  }

  void Centroid::setScoreMatrix( const int sm[MAT][MAT], double T ) {
    this -> T = T;
    this -> isPssm = false;
    for ( int n=0; n<MAT; ++n )
      for ( int m=0; m<MAT; ++m ) {
	match_score[n][m] = EXP ( sm[ n ][ m ] / T );
      }
  }

  void Centroid::setPssm( const int pssm[][MAT], const unsigned int qsize, double T, const OneQualityExpMatrix& oqem, const uchar* sequenceBeg, const uchar* qualityBeg ) {
    this->T = T;
    this -> isPssm = true;
    pssmExp.resize( qsize * MAT );
    pssmExp2 = reinterpret_cast<ScoreMatrixRowDouble*> ( &pssmExp[0] );

    if( oqem ){  // fast special case
      makePositionSpecificExpMatrix( oqem, sequenceBeg, sequenceBeg + qsize,
                                     qualityBeg, &pssmExp[0] );
    }
    else{  // slow general case
      for ( unsigned i=0; i<qsize; ++i ) {
        for ( unsigned j=0; j<MAT; ++j ) {
          pssmExp2[ i ][ j ] = EXP ( pssm[ i ][ j ] / T );
        }
      }
    }
  }

  void Centroid::initForwardMatrix(){
    scale.assign ( numAntidiagonals, 1.0 ); // scaling
    std::size_t n = xa.scoreEndIndex( numAntidiagonals );

    if ( fM.size() < n ) {
      fM.resize( n );
      fD.resize( n );
      fI.resize( n );
      fP.resize( n );
    }

    fM[3] = 1;
  }
  
  void Centroid::initBackwardMatrix(){
    pp.resize( fM.size() );
    mD.assign( numAntidiagonals, 0.0 );
    mI.assign( numAntidiagonals, 0.0 );
    mX1.assign ( numAntidiagonals, 1.0 );
    mX2.assign ( numAntidiagonals, 1.0 );

    std::size_t n = xa.scoreEndIndex( numAntidiagonals );
    bM.assign( n, 0.0 );
    bD.assign( n, 0.0 );
    bI.assign( n, 0.0 );
    bP.assign( n, 0.0 );
  }

  void Centroid::initDecodingMatrix(){
    X.resize( fM.size() );
  }

  void Centroid::updateScore( double score, size_t antiDiagonal, size_t cur ){
    if( bestScore < score ){
      bestScore = score;
      bestAntiDiagonal = antiDiagonal;
      bestPos1 = cur;
    }
  }

  // xxx this will go wrong for non-delimiters with severe mismatch scores
  static bool isDelimiter(uchar c, const double *expScores) {
    return expScores[c] <= 0.0;
  }

  double Centroid::forward( const uchar* seq1, const uchar* seq2,
			    size_t start1, size_t start2,
			    bool isForward, int globality,
			    const GeneralizedAffineGapCosts& gap ){

    //std::cout << "[forward] start1=" << start1 << "," << "start2=" << start2 << "," << "isForward=" << isForward << std::endl;
    const int seqIncrement = isForward ? 1 : -1;

    initForwardMatrix();

    if( globality ) {
      Z = 0.0;
      const uchar* s1 = seqPtr( seq1, start1, isForward, 1 );
      if (! isPssm) {
	const uchar* s2 = seqPtr( seq2, start2, isForward, 1 );
	if( isDelimiter(*s1, *match_score) || isDelimiter(*s2, *match_score) ){
	  Z = fM[3];
	}
      }else{
	const double (*p2)[MAT] = seqPtr( pssmExp2, start2, isForward, 1 );
	if( isDelimiter(*s1, *(pssmExp2+start2)) || isDelimiter(0, *p2) ){
	  Z = fM[3];
	}
      }
    }else{
      Z = fM[3];
    }

    const bool isAffine = gap.isAffine();
    const int E = gap.delExtend;
    const int F = gap.delExist + gap.delExtend;
    const int EI = gap.insExtend;
    const int FI = gap.insExist + gap.insExtend;
    const int P = gap.pairExtend;
    const int Q = gap.delExist + gap.pairExtend;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eEI = EXP ( - EI / T );
    const double eFI = EXP ( - FI / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );

    assert( gap.insExist == gap.delExist || eQ <= 0.0 );

    for( size_t k = 3; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      double sum_f = 0.0; // sum of forward values
      const size_t k1 = k - 1;
      const size_t k2 = k - 2;
      const size_t loopBeg = xa.seq1start( k );
      const std::size_t seq2pos = k2 - loopBeg;
      const double scale12 = 1.0 / ( scale[k1] * scale[k2] ); // scaling factor
      const double scale1  = 1.0 / scale[k1];

      const double seF = eF * scale1;
      const double seE = eE * scale1;
      const double seFI = eFI * scale1;
      const double seEI = eEI * scale1;
      const double seQ = eQ * scale12;
      const double seP = eP * scale12;

      const std::size_t scoreEnd = xa.scoreEndIndex( k );
      double* fM0 = &fM[ scoreEnd ];
      double* fD0 = &fD[ scoreEnd ];
      double* fI0 = &fI[ scoreEnd ];
      double* fP0 = &fP[ scoreEnd ];

      const double* const fM0end = fM0 + xa.numCellsAndPads( k );
      const std::size_t horiBeg = xa.hori( k, loopBeg );
      const std::size_t diagBeg = xa.diag( k, loopBeg );
      const double* fM1 = &fM[ horiBeg ];
      const double* fD1 = &fD[ horiBeg ];
      const double* fI1 = &fI[ horiBeg ];
      const double* fP1 = &fP[ horiBeg ];
      const double* fM2 = &fM[ diagBeg ];
      const double* fD2 = &fD[ diagBeg ];
      const double* fI2 = &fI[ diagBeg ];
      const double* fP2 = &fP[ diagBeg ];
      double xM1 = *fM1, xD1 = *fD1, xI1 = *fI1, xP1 = *fP1;

      const uchar* s1 = seqPtr( seq1, start1, isForward, loopBeg );

      *fM0++ = *fD0++ = *fI0++ = *fP0++ = 0.0;  // add one pad cell

      if (! isPssm) {	  
	const uchar* s2 = seqPtr( seq2, start2, isForward, seq2pos );
	while (1) {	// start: inner most loop
	  const double S = match_score[ *s1 ][ *s2 ] * scale12; // (loopBeg, k - loopBeg)
	  const double xM2 = *fM2, xD2 = *fD2, xI2 = *fI2, xP2 = *fP2;
	  *fD0 = ( xM1 ) * seF + ( xD1 + xP1 ) * seE;
	  xM1 = *++fM1; xD1 = *++fD1; xI1 = *++fI1; xP1 = *++fP1;
	  *fI0 = ( xM1 + xD1 ) * seFI + ( xI1 + xP1 ) * seEI;
	  *fM0 = ( xM2 + xD2 + xI2 + xP2 ) * S;
	  *fP0 = xM2 * seQ + xP2 * seP;
	  fM2++; fD2++; fI2++; fP2++;
	  sum_f += *fM0;
	  if( globality && (isDelimiter(*(s2+seqIncrement), *match_score) ||
			    isDelimiter(*(s1+seqIncrement), *match_score)) ){
	    Z += *fM0 + *fD0 + *fI0 + *fP0;
	  }
	  fM0++; fD0++; fI0++; fP0++;
	  if (fM0 == fM0end) break;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}	// end: inner most loop
      } // end: if (! isPssm)
      else { // when pssm2 is given
	const double (*p2)[MAT] = seqPtr( pssmExp2, start2, isForward, seq2pos );
	  
	if (isAffine) {
	  while (1) { // start: inner most loop
	    const double S = (*p2)[ *s1 ]  * scale12;
	    const double xM2 = *fM2, xD2 = *fD2, xI2 = *fI2;
	    *fD0 = ( xM1 ) * seF + ( xD1 ) * seE;
	    xM1 = *++fM1; xD1 = *++fD1; xI1 = *++fI1; 
	    *fI0 = ( xM1 + xD1 ) * seF + ( xI1 ) * seE;
	    *fM0 = ( xM2 + xD2 + xI2 ) * S;
	    fM2++; fD2++; fI2++; 
	    sum_f += *fM0;
	    if ( globality && (isDelimiter(0, *(p2+seqIncrement)) ||
			       isDelimiter(*(s1+seqIncrement), *(pssmExp2+start2))) ){
	      Z += *fM0 + *fD0 + *fI0;
	    }
	    fM0++; fD0++; fI0++; 
	    if (fM0 == fM0end) break;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }	// end: inner most loop
	}else{
	  while (1) { // start: inner most loop
	    const double S = (*p2)[ *s1 ]  * scale12;
	    const double xM2 = *fM2, xD2 = *fD2, xI2 = *fI2, xP2 = *fP2;
	    *fD0 = ( xM1 ) * seF + ( xD1 + xP1 ) * seE;
	    xM1 = *++fM1; xD1 = *++fD1; xI1 = *++fI1; xP1 = *++fP1;
	    *fI0 = ( xM1 + xD1 ) * seFI + ( xI1 + xP1 ) * seEI;
	    *fM0 = ( xM2 + xD2 + xI2 + xP2 ) * S;
	    *fP0 = xM2 * seQ + xP2 * seP;
	    fM2++; fD2++; fI2++; fP2++;
	    sum_f += *fM0;
	    if ( globality && (isDelimiter(0, *(p2+seqIncrement)) ||
			       isDelimiter(*(s1+seqIncrement), *(pssmExp2+start2))) ){
	      Z += *fM0 + *fD0 + *fI0 + *fP0;
	    }
	    fM0++; fD0++; fI0++; fP0++;
	    if (fM0 == fM0end) break;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }	// end: inner most loop
	}
      }
      if( !globality ) Z += sum_f;
      scale[k] = sum_f + 1.0;  // seems ugly
      Z /= scale[k]; // scaling
    } // k
    //std::cout << "Z=" << Z << std::endl;
    assert( Z > 0.0 );
    return log(Z);
  }

  // added by M. Hamada
  // compute posterior probabilities while executing backward algorithm 
  // posterior probabilities are stored in pp
  double Centroid::backward( const uchar* seq1, const uchar* seq2,
			     size_t start1, size_t start2,
			     bool isForward, int globality,
			     const GeneralizedAffineGapCosts& gap ){

    //std::cout << "[backward] start1=" << start1 << "," << "start2=" << start2 << "," << "isForward=" << isForward << std::endl;
    const int seqIncrement = isForward ? 1 : -1;

    initBackwardMatrix();

    const bool isAffine = gap.isAffine();
    const int E = gap.delExtend;
    const int F = gap.delExist + gap.delExtend;
    const int EI = gap.insExtend;
    const int FI = gap.insExist + gap.insExtend;
    const int P = gap.pairExtend;
    const int Q = gap.delExist + gap.pairExtend;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eEI = EXP ( - EI / T );
    const double eFI = EXP ( - FI / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );

    assert( gap.insExist == gap.delExist || eQ <= 0.0 );

    double d1 = 1.0;

    for( size_t k = numAntidiagonals-1; k > 2; --k ){  // loop over antidiagonals
      const size_t k1 = k - 1;  
      const size_t k2 = k - 2;
      const size_t loopBeg = xa.seq1start( k );
      const std::size_t seq2pos = k2 - loopBeg;
      const double scale12 = 1.0 / ( scale[k1] * scale[k2] ); // scaling factor
      const double scale1  = 1.0 / scale[k1];

      const double seF = eF * scale1;
      const double seE = eE * scale1;
      const double seFI = eFI * scale1;
      const double seEI = eEI * scale1;
      const double seQ = eQ * scale12;
      const double seP = eP * scale12;

      const std::size_t scoreEnd = xa.scoreEndIndex( k );
      double* bM0 = &bM[ scoreEnd + 1 ];
      double* bD0 = &bD[ scoreEnd + 1 ];
      double* bI0 = &bI[ scoreEnd + 1 ];
      double* bP0 = &bP[ scoreEnd + 1 ];

      double* pp0 = &pp[ scoreEnd ];

      const double* fM0 = &fM[ scoreEnd + 1 ];
      const double* fD0 = &fD[ scoreEnd + 1 ];
      const double* fI0 = &fI[ scoreEnd + 1 ];
      const double* fP0 = &fP[ scoreEnd + 1 ];

      const double* const fM0end = fM0 + xa.numCellsAndPads( k ) - 1;
      const std::size_t horiBeg = xa.hori( k, loopBeg );
      const std::size_t diagBeg = xa.diag( k, loopBeg );
      double* bM1 = &bM[ horiBeg ];
      double* bD1 = &bD[ horiBeg ];
      double* bI1 = &bI[ horiBeg ];
      double* bP1 = &bP[ horiBeg ];
      double* bM2 = &bM[ diagBeg ];
      double* bD2 = &bD[ diagBeg ];
      double* bI2 = &bI[ diagBeg ];
      double* bP2 = &bP[ diagBeg ];

      int i = loopBeg; int j = seq2pos;

      d1 /= scale[k];

      const uchar* s1 = seqPtr( seq1, start1, isForward, loopBeg );

      if (! isPssm ) {
	const uchar* s2 = seqPtr( seq2, start2, isForward, seq2pos );

	do{ // inner most loop
	  if( globality ){
	    if( isDelimiter(*(s2+seqIncrement), *match_score) ||
		isDelimiter(*(s1+seqIncrement), *match_score) ){
	      *bM0 += d1;  *bD0 += d1;  *bI0 += d1;  *bP0 += d1;
	    }
	  }else{
	    *bM0 += d1;
	  }
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
	  const double tmp5 = tmp4 * seFI; 
	  const double tmp6 = tmp4 * seEI;
	  *bM1 += tmp5; 
	  *bD1 += tmp5; 
	  *bI1 += tmp6; 
	  *bP1 += tmp6; 

	  *pp0 = *fM0 * *bM0 / Z; 
	  double probd = *fD0 * *bD0 / Z;
	  double probi = *fI0 * *bI0 / Z;
	  double probp = *fP0 * *bP0 / Z;
	  mD[ i ] += probd + probp; 
	  mI[ j ] += probi + probp; 
	  mX1 [ i ] -= ( *pp0 + probd + probp );
	  mX2 [ j ] -= ( *pp0 + probi + probp );
	  i++; j--;
	  // iteration
	  bM2++; bD2++; bI2++; bP2++;
	  bM0++; bD0++; bI0++; bP0++;
	  fM0++; fD0++; fI0++; fP0++; 
	  pp0++; 
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( fM0 != fM0end ); // inner most loop end;
      } // if (!ppsm2)
      else {
	const double (*p2)[MAT] = seqPtr( pssmExp2, start2, isForward, seq2pos );
	  
	if (isAffine) {
	  do{ // inner most loop
	    if( globality ){
	      if( isDelimiter(0, *(p2+seqIncrement)) ||
		  isDelimiter(*(s1+seqIncrement), *(pssmExp2+start2)) ){
		*bM0 += d1;  *bD0 += d1;  *bI0 += d1;
	      }
	    }else{
	      *bM0 += d1;
	    }
	    const double S = ( *p2 )[ *s1 ];
	    const double tmp1 = *bM0 * S * scale12;
	    *bM2 += tmp1; 
	    *bD2 += tmp1; 
	    *bI2 += tmp1; 
	    const double tmp3 = *bD0;
	    *bM1++ += tmp3 * seF; 
	    *bD1++ += tmp3 * seE; 
	    bI1++;
	    const double tmp4 = *bI0;
	    const double tmp5 = tmp4 * seF; 
	    const double tmp6 = tmp4 * seE;
	    *bM1 += tmp5; 
	    *bD1 += tmp5; 
	    *bI1 += tmp6; 

	    double prob = *fM0 * *bM0 / Z; 
	    *pp0 = prob;
	    double probd = *fD0 * *bD0 / Z;
	    double probi = *fI0 * *bI0 / Z;
	    mD[ i ] += probd; 
	    mI[ j ] += probi; 
	    mX1 [ i ] -= ( prob + probd );
	    mX2 [ j ] -= ( prob + probi );
	    i++; j--;
	    // iteration
	    bM2++; bD2++; bI2++;
	    bM0++; bD0++; bI0++;
	    fM0++; fD0++; fI0++; 
	    pp0++;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }while( fM0 != fM0end ); // inner most loop end;
	}else{
	  do{
	    if( globality ){
	      if( isDelimiter(0, *(p2+seqIncrement)) ||
		  isDelimiter(*(s1+seqIncrement), *(pssmExp2+start2)) ){
		*bM0 += d1;  *bD0 += d1;  *bI0 += d1;  *bP0 += d1;
	      }
	    }else{
	      *bM0 += d1;
	    }
	    const double S = ( *p2 )[ *s1 ];
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
	    const double tmp5 = tmp4 * seFI; 
	    const double tmp6 = tmp4 * seEI;
	    *bM1 += tmp5; 
	    *bD1 += tmp5; 
	    *bI1 += tmp6; 
	    *bP1 += tmp6; 

	    double prob = *fM0 * *bM0 / Z; 
	    *pp0 = prob;
	    double probd = *fD0 * *bD0 / Z;
	    double probi = *fI0 * *bI0 / Z;
	    double probp = *fP0 * *bP0 / Z;
	    mD[ i ] += probd + probp; 
	    mI[ j ] += probi + probp; 
	    mX1 [ i ] -= ( prob + probd + probp );
	    mX2 [ j ] -= ( prob + probi + probp );
	    i++; j--;
	    // iteration
	    bM2++; bD2++; bI2++; bP2++;
	    bM0++; bD0++; bI0++; bP0++;
	    fM0++; fD0++; fI0++; fP0++; 
	    pp0++;
	    s1 += seqIncrement;
	    p2 -= seqIncrement;
	  }while( fM0 != fM0end ); // inner most loop end;
	}
      }
    }

    d1 /= scale[2];
    if( globality ){
      // something is missing (but I think it doesn't matter)
    }else{
      bM[3] += d1;
    }

    // Modify for test start
    //ExpectedCount ec;
    //computeExpectedCounts ( seq1, seq2, start1, start2, isForward, gap, ec );
    //ec.write (std::cerr, Z);
    // Modify for test end
    return log( bM[3] );
  }

  double Centroid::dp( double gamma ){
    if (outputType == 5 ) return dp_centroid( gamma );
    else if (outputType == 6 ) return dp_ama (gamma);
    return 0;
  }

  void Centroid::traceback( std::vector< SegmentPair >& chunks,
			    double gamma ) const{
    if (outputType==5) traceback_centroid( chunks, gamma );
    else if (outputType==6) traceback_ama( chunks, gamma);
  }

  double Centroid::dp_centroid( double gamma ){

    initDecodingMatrix();

    for( size_t k = 3; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const std::size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      const double* P0 = &pp[ scoreEnd ];
      size_t cur = xa.seq1start( k );

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const double* X1 = &X[xa.hori(k, cur)];
      const double* X2 = &X[xa.diag(k, cur)];

      *X0++ = -DINF;		// add one pad cell

      do{
	const double s = ( gamma + 1 ) * ( *P0++ ) - 1;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max( std::max( oldX1, *X1 ), *X2++ + s );
	//assert ( score >= 0 );
	updateScore ( score, k, cur );
	*X0++ = score;
	cur++;
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

    while( k > 2 ){
      const int m =
	maxIndex( diagx( X, k, i ) + ( gamma + 1 ) * cellx( pp, k, i ) - 1,
                  horix( X, k, i ),
                  vertx( X, k, i ) );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 2 ){
	chunks.push_back( SegmentPair( i, k - i - 2, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  double Centroid::dp_ama( double gamma ){

    initDecodingMatrix();

    for( size_t k = 3; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const std::size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      const double* P0 = &pp[ scoreEnd ];
      size_t cur = xa.seq1start( k );
      size_t seq2pos = k - 2 - cur;

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const double* X1 = &X[ xa.hori(k, cur) ];
      const double* X2 = &X[ xa.diag(k, cur) ];

      *X0++ = -DINF;		// add one pad cell

      do{
	const double s = 2 * gamma * *P0++ - ( mX1[ cur ] + mX2[ seq2pos ] );
	const double oldX1 = *X1++;  // Added by MCF
	const double u = gamma * mD[ cur ] - mX1[ cur ];
	const double t = gamma * mI[ seq2pos ] - mX2[ seq2pos ];
	const double score = std::max( std::max( oldX1 + u, *X1 + t), *X2++ + s );
	updateScore ( score, k, cur );
	*X0++ = score;
	cur++;
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

    while( k > 2 ){
      const size_t j = k - i - 2;
      const double s = 2 * gamma * cellx( pp, k, i ) - ( mX1[ i ] + mX2[ j ] );
      const double t = gamma * mI[ j ] - mX2[ j ];
      const double u = gamma * mD[ i ] - mX1[ i ];
      const int m =
	maxIndex( diagx( X, k, i ) + s,
                  horix( X, k, i ) + u,
                  vertx( X, k, i ) + t );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 2 ){
	chunks.push_back( SegmentPair( i, k - i - 2, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  // Return an ASCII code representing an error probability.  The
  // printable codes are 33--126, but here we use a maximum of 125, so
  // that 126 is reserved for special cases.
  static uchar asciiProbability( double probCorrect ){
    assert( probCorrect >= 0 );
    assert( probCorrect <= 1 );
    double e = 1 - probCorrect;
    double f = std::max( e, 1e-10 );  // avoid overflow errors
    double g = -10 * std::log10(f);
    int i = static_cast<int>(g);  // round fractions down
    int j = i + 33;
    int k = std::min( j, 125 );
    return static_cast<uchar>(k);
  }

  static void getGapAmbiguities( std::vector<uchar>& ambiguityCodes,
                                    const std::vector<double>& probs,
                                    size_t rbeg, size_t rend ){
    for( size_t i = rbeg; i > rend; --i ){
      ambiguityCodes.push_back( asciiProbability( probs[ i ] ) );
    }
  }

  // Added by MCF:
  void Centroid::getColumnAmbiguities( std::vector<uchar>& ambiguityCodes,
                                       const std::vector<SegmentPair>& chunks,
                                       bool isForward ){
    for( CI(SegmentPair) i = chunks.begin(); i < chunks.end(); ++i ){
      size_t seq1pos = i->end1();
      size_t seq2pos = i->end2();

      for( size_t j = 0; j < i->size; ++j ){
	double p = cellx( pp, seq1pos + seq2pos + 2, seq1pos );
	ambiguityCodes.push_back( asciiProbability(p) );
	--seq1pos;
	--seq2pos;
      }

      CI(SegmentPair) j = i + 1;
      size_t end1 = (j < chunks.end()) ? j->end1() : 0;
      size_t end2 = (j < chunks.end()) ? j->end2() : 0;

      // ASSUMPTION: if there is an insertion adjacent to a deletion,
      // the deletion will get printed first.
      if( isForward ){
        getGapAmbiguities( ambiguityCodes, mI, seq2pos, end2 );
        getGapAmbiguities( ambiguityCodes, mD, seq1pos, end1 );
      }
      else{
        getGapAmbiguities( ambiguityCodes, mD, seq1pos, end1 );
        getGapAmbiguities( ambiguityCodes, mI, seq2pos, end2 );
      }
    }
  }

  double Centroid::logPartitionFunction() const{
    assert( Z > 0 );
    double x = std::log(Z);
    for( std::size_t k = 3; k < numAntidiagonals; ++k ){
      x += std::log( scale[k] );
    }
    return T * x;
  }

  void Centroid::computeExpectedCounts ( const uchar* seq1, const uchar* seq2,
					 size_t start1, size_t start2, bool isForward,
					 const GeneralizedAffineGapCosts& gap, 
					 ExpectedCount& c ) const{

    const int seqIncrement = isForward ? 1 : -1;

    const bool isAffine = gap.isAffine();
    const int E = gap.delExtend;
    const int F = gap.delExist + gap.delExtend;
    const int EI = gap.insExtend;
    const int FI = gap.insExist + gap.insExtend;
    const int P = gap.pairExtend;
    const int Q = gap.delExist + gap.pairExtend;
    const double eE = EXP ( - E / T );
    const double eF = EXP ( - F / T );
    const double eEI = EXP ( - EI / T );
    const double eFI = EXP ( - FI / T );
    const double eP = EXP ( - P / T );
    const double eQ = EXP ( - Q / T );

    assert( gap.insExist == gap.delExist || eQ <= 0.0 );

    c.SQ = 1; 

    for( size_t k = 3; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t k1 = k - 1;
      const size_t k2 = k - 2;
      const size_t loopBeg = xa.seq1start( k );
      const std::size_t seq2pos = k2 - loopBeg;

      const std::size_t scoreEnd = xa.scoreEndIndex( k );
      const double* fM0 = &fM[ scoreEnd + 1 ];
      const double* fD0 = &fD[ scoreEnd + 1 ];
      const double* fI0 = &fI[ scoreEnd + 1 ];
      const double* fP0 = &fP[ scoreEnd + 1 ];
      const double* bM0 = &bM[ scoreEnd + 1 ];
      const double* bD0 = &bD[ scoreEnd + 1 ];
      const double* bI0 = &bI[ scoreEnd + 1 ];
      const double* bP0 = &bP[ scoreEnd + 1 ];

      const double scale12 = 1.0 / ( scale[k1] * scale[k2] );
      const double scale1  = 1.0 / scale[k1];
      const double scale0  = 1.0 / scale[k];

      const double* const fM0end = fM0 + xa.numCellsAndPads( k ) - 1;
      const uchar* s1 = seqPtr( seq1, start1, isForward, loopBeg );
      const uchar* s2 = seqPtr( seq2, start2, isForward, seq2pos );
      //assert ( *s1 < MAT && *s2 < MAT );
      const std::size_t horiBeg = xa.hori( k, loopBeg );
      const std::size_t diagBeg = xa.diag( k, loopBeg );
      const double* fM1 = &fM[ horiBeg ];
      const double* fD1 = &fD[ horiBeg ];
      const double* fI1 = &fI[ horiBeg ];
      const double* fP1 = &fP[ horiBeg ];
      const double* fM2 = &fM[ diagBeg ];
      const double* fD2 = &fD[ diagBeg ];
      const double* fI2 = &fI[ diagBeg ];
      const double* fP2 = &fP[ diagBeg ];

      // Modify start for PSSM
      if (! isPssm ) {
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
	  c.MI += ( *fM1 * eFI )  * tmp4;
	  c.DI += ( *fD1 * eFI )  * tmp4;
	  c.PI += ( *fP1 * eEI )  * tmp4;
	  c.II += ( *fI1 * eEI )  * tmp4;

	  fM2++; fD2++; fI2++; fP2++;
	  fM0++; fD0++; fI0++; fP0++;
	  bM0++; bD0++; bI0++; bP0++;
	    
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	} while( fM0 != fM0end ); // inner most loop end;
      } // if (!ppsm2)
      else {
	const double (*p2)[MAT] = seqPtr( pssmExp2, start2, isForward, seq2pos );

	if (isAffine) {
	  do{ // inner most loop
	    const double S = ( *p2 )[ *s1 ];
	    c.emit[*s1][*s2] += ( *fM0 * *bM0 ) ;

	    const double tmp1 = S * *bM0 * scale12;

	    c.MM += *fM2 * tmp1; 
	    c.DM += *fD2 * tmp1; 
	    c.IM += *fI2 * tmp1; 
	    c.MQ += *fM0;

	    const double tmp3 = *bD0 * scale1;
	    c.MD += ( *fM1 * eF ) * tmp3;
	    c.DD += ( *fD1 * eE ) * tmp3;

	    fM1++; fD1++; fI1++;

	    const double tmp4 = *bI0 * scale1;
	    c.MI += ( *fM1 * eF )  * tmp4;
	    c.DI += ( *fD1 * eF )  * tmp4;
	    c.II += ( *fI1 * eE )  * tmp4;

	    fM2++; fD2++; fI2++;
	    fM0++; fD0++; fI0++;
	    bM0++; bD0++; bI0++;
	    
	    s1 += seqIncrement;
	    s2 -= seqIncrement;  // xxx p2 ???
	  } while( fM0 != fM0end ); // inner most loop end;
	}else{
	  do{ // inner most loop
	    const double S = ( *p2 )[ *s1 ];
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
	    c.MI += ( *fM1 * eFI )  * tmp4;
	    c.DI += ( *fD1 * eFI )  * tmp4;
	    c.PI += ( *fP1 * eEI )  * tmp4;
	    c.II += ( *fI1 * eEI )  * tmp4;

	    fM2++; fD2++; fI2++; fP2++;
	    fM0++; fD0++; fI0++; fP0++;
	    bM0++; bD0++; bI0++; bP0++;
	    
	    s1 += seqIncrement;
	    s2 -= seqIncrement;  // xxx p2 ???
	  } while( fM0 != fM0end ); // inner most loop end;
	}
      }
      // Modify end for PSSM

      c.MQ *= scale0;
      c.SQ *= scale0;
    }
  }
}  // end namespace cbrc
