// Copyright 2008, 2009, 2011 Michiaki Hamada

#ifndef CENTROID_HH
#define CENTROID_HH
#include "GappedXdropAligner.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include "OneQualityScoreMatrix.hh"
#include <cassert>
#include <stddef.h>  // size_t
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

    Centroid( const GappedXdropAligner& xa_ );

    // Setters
    void setScoreMatrix( const int sm[MAT][MAT], double T );
    void setPssm ( const int pssm[][MAT], const unsigned int qsize, double T,
                   const OneQualityExpMatrix& oqem,
                   const uchar* sequenceBeg, const uchar* qualityBeg );
    void setOutputType( int m ) { outputType = m; }

    void reset( ) { 
      lastAntiDiagonal = xa.numAntidiagonals () - 3;
      bestScore = 0;
      bestAntiDiagonal = 0;
      bestPos1 =0;
    }

    double forward( const uchar* seq1, const uchar* seq2, 
		    size_t start1, size_t start2, bool isForward,
		    const GeneralizedAffineGapCosts& gap );
    
    double backward( const uchar* seq1, const uchar* seq2, 
		     size_t start1, size_t start2, bool isForward,
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
                               bool isForward );

    // Added by MH (2008/10/10) : compute expected counts for transitions and emissions
    void computeExpectedCounts ( const uchar* seq1, const uchar* seq2,
				 size_t start1, size_t start2, bool isForward,
				 const GeneralizedAffineGapCosts& gap,
				 ExpectedCount& count ) const;

  private:
    const GappedXdropAligner& xa;
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

    // The next two functions use "antidiagonal + 2".  This is because
    // GappedXdropAligner uses 2 extra padding antidiagonals at the
    // beginning, whereas Centroid does not.  (Centroid should probably
    // be changed so that it does: that makes it easier to look-back by
    // up to 2 antidiagonals.)

    size_t seq1start( size_t antiDiagonal ) const{
      return xa.seq1start( antiDiagonal + 2 );
    }

    // The next function has "- 1", because GappedXdropAligner uses one
    // padding cell per antidiagonal, whereas Centroid does not.
    // (Centroid should probably be changed so that it does: then we
    // would not need special code for boundary cases.)

    size_t numCells( size_t antiDiagonal ) const{
      return xa.numCellsAndPads( antiDiagonal + 2 ) - 1;
    }

    size_t seq1end( size_t antiDiagonal ) const{
      return seq1start( antiDiagonal ) + numCells( antiDiagonal );
    }

    // get DP matrix value at the given position
    double cellx( const dmatrix_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      assert( seq1pos >= seq1start( antiDiagonal ) );
      assert( seq1pos < seq1end( antiDiagonal ) );
      return matrix[ antiDiagonal ][ seq1pos - seq1start( antiDiagonal ) ];
    }

    // get DP matrix value "left of" the given position
    double horix( const dmatrix_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      assert( antiDiagonal > 0 );
      if( seq1pos > seq1start( antiDiagonal-1 ) )
        return cellx( matrix, antiDiagonal-1, seq1pos-1 );
      else return -INF;
    }

    // get DP matrix value "above" the given position
    double vertx( const dmatrix_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      assert( antiDiagonal > 0 );
      if( seq1pos < seq1end( antiDiagonal-1 ) )
        return cellx( matrix, antiDiagonal-1, seq1pos );
      else return -INF;
    }

    // get DP matrix value "diagonal from" the given position
    double diagx( const dmatrix_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      if( antiDiagonal > 1 &&
          seq1pos > seq1start( antiDiagonal-2 ) &&
          seq1pos <= seq1end( antiDiagonal-2 ) )
        return cellx( matrix, antiDiagonal-2, seq1pos-1 );
      else return -INF;
    }

    double diag( const dmatrix_t& matrix,
		 size_t antiDiagonal, size_t seq1pos ) const;

    // get a pointer into a sequence, taking start and direction into account
    template < typename T >
    static const T* seqPtr( const T* seq, size_t start,
                            bool isForward, size_t pos ){
      if( isForward ) return &seq[ start + pos - 1 ];
      else            return &seq[ start - pos ];
    }
  };

}  // end namespace cbrc
#endif  // CENTROID_HH
