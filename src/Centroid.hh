// Copyright 2008, 2009 Michiaki Hamada

#ifndef CENTROID_HH
#define CENTROID_HH
#include "XdropAligner.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include <vector>
#include <iostream> // for debug

namespace cbrc{

  struct ExpectedCount{
  public:
    enum { MAT = 64 };
    double emit[MAT][MAT];
    double MM, MD, MP, MI, MQ;
    double DD, DM, DI;
    double PP, PM, PD, PI;
    double II, IM;
    double SM, SD, SP, SI, SQ;
  public:
    ExpectedCount ();
    std::ostream& write (std::ostream& os, double Z) const;
  };
  /**
   * (1) Forward and backward algorithm on the DP region given by Xdrop algorithm
   * (2) \gamma-centroid decoding
   */
  class Centroid{
  public:
    enum { MAT = 64 };
    Centroid( const XdropAligner& xa_, const int sm[MAT][MAT], double T_ = 1.0 );
    void reset( ) { 
      lastAntiDiagonal = xa.offsets.size () - 1;
      bestScore = 0;
      bestAntiDiagonal = 0;
      bestPos1 =0;
    }

    typedef unsigned char uchar;
    double forward( const uchar* seq1, const uchar* seq2, 
		    size_t start1, size_t start2, XdropAligner::direction dir,
		    const GeneralizedAffineGapCosts& gap );
    
    double backward( const uchar* seq1, const uchar* seq2, 
		     size_t start1, size_t start2, XdropAligner::direction dir,
		     const GeneralizedAffineGapCosts& gap );
    double dp( double gamma );
    void traceback( std::vector< SegmentPair >& chunks, double gamma ) const;

    // Added by MCF: get the probabilities of each match in each chunk:
    void chunkProbabilities( std::vector< double >& probs,
			     const std::vector< SegmentPair >& chunks );

    // Added by MH (2008/10/10) : compute expected counts for transitions and emissions
    void computeExpectedCounts ( const uchar* seq1, const uchar* seq2,
				 size_t start1, size_t start2, XdropAligner::direction dir,
				 const GeneralizedAffineGapCosts& gap,
				 ExpectedCount& count ) const;

  private:
    const XdropAligner& xa;
    double T; // temperature
    size_t lastAntiDiagonal;
    double match_score[ MAT ][ MAT ]; // pre-computed match score
    typedef std::vector< std::vector< double > > dmatrix_t;
    typedef std::vector< double > dvec_t;

    dmatrix_t fM; // f^M(i,j)
    dmatrix_t fD; // f^D(i,j)
    dmatrix_t fI; // f^I(i,j)
    dmatrix_t fP; // f^P(i,j)

    double    Z; // partion function of forward values

    dmatrix_t bM; // b^M(i,j)
    dmatrix_t bD; // b^D(i,j)
    dmatrix_t bI; // b^I(i,j)
    dmatrix_t bP; // b^P(i,j)

    dmatrix_t pp; // posterior match probabilities

    dmatrix_t X; // DP tables for $gamma$-decoding

    dvec_t scale; // scale[n] is a scaling factor for the n-th anti-diagonal

    double bestScore;
    size_t bestAntiDiagonal;
    size_t bestPos1;

    void initForwardMatrix();
    void initBackwardMatrix();
    void initDecodingMatrix();

    void updateScore( double score, size_t antiDiagonal, size_t cur );

    double cell( const dmatrix_t& matrix,
		 size_t antiDiagonal, size_t seq1pos ) const;
    double diag( const dmatrix_t& matrix,
		 size_t antiDiagonal, size_t seq1pos ) const;
  };

}  // end namespace cbrc
#endif  // CENTROID_HH
