// Copyright 2008, 2009, 2011 Michiaki Hamada

#ifndef CENTROID_HH
#define CENTROID_HH
#include "XdropAligner.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include "OneQualityScoreMatrix.hh"
#include <vector>
#include <iostream> // for debug

namespace cbrc{

  typedef unsigned char uchar;

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

    Centroid( const XdropAligner& xa_ );

    // Setters
    void setScoreMatrix( const int sm[MAT][MAT], double T );
    void setPssm ( const int pssm[][MAT], const unsigned int qsize, double T,
                   const OneQualityExpMatrix& oqem,
                   const uchar* sequenceBeg, const uchar* qualityBeg );
    void setOutputType( int m ) { outputType = m; }

    void reset( ) { 
      lastAntiDiagonal = xa.offsets.size () - 1;
      bestScore = 0;
      bestAntiDiagonal = 0;
      bestPos1 =0;
    }

    double forward( const uchar* seq1, const uchar* seq2, 
		    size_t start1, size_t start2, XdropAligner::direction dir,
		    const GeneralizedAffineGapCosts& gap );
    
    double backward( const uchar* seq1, const uchar* seq2, 
		     size_t start1, size_t start2, XdropAligner::direction dir,
		     const GeneralizedAffineGapCosts& gap );

    double dp( double gamma );
    void traceback( std::vector< SegmentPair >& chunks, double gamma ) const;

    double dp_centroid( double gamma );
    void traceback_centroid( std::vector< SegmentPair >& chunks, double gamma ) const;

    double dp_ama( double gamma );
    void traceback_ama( std::vector< SegmentPair >& chunks, double gamma ) const;

    // Added by MCF: get the probability of each column in the alignment:
    void getColumnAmbiguities( std::vector< uchar >& ambiguityCodes,
                               const std::vector< SegmentPair >& chunks,
                               XdropAligner::direction dir );

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
    //const int (*pssm2)[MAT];
    bool isPssm;
    std::vector<double> pssmExp; //
    /* const */ double (*pssmExp2)[MAT]; // pre-computed pssm for prob align
    int outputType;

    typedef std::vector< std::vector< double > > dmatrix_t;
    typedef std::vector< double > dvec_t;

    dmatrix_t fM; // f^M(i,j)
    dmatrix_t fD; // f^D(i,j) Ix 
    dmatrix_t fI; // f^I(i,j) Iy
    dmatrix_t fP; // f^P(i,j)

    double    Z; // partion function of forward values

    dmatrix_t bM; // b^M(i,j)
    dmatrix_t bD; // b^D(i,j)
    dmatrix_t bI; // b^I(i,j)
    dmatrix_t bP; // b^P(i,j)

    dmatrix_t pp; // posterior match probabilities

    dvec_t mD;
    dvec_t mI;
    dvec_t mX1; 
    dvec_t mX2;

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
