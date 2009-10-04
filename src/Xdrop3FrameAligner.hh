// Copyright 2009 Toshiyuki Sato

// X-drop algorithm for gapped extension of alignments, with
// frameshifts.  The underlying algorithm is "three-frame alignment",
// as implemented in FASTX and described in: "Aligning a DNA sequence
// with a protein sequence", Z Zhang, WR Pearson & W Miller, J Comput
// Biol. 1997 4(3):339-49.  So this code often gives the same results
// as FASTX, which is a good way of testing it.

#ifndef XDROP3FRAMEALIGNER_HH
#define XDROP3FRAMEALIGNER_HH
#include "XdropAligner.hh"

namespace cbrc{

class Xdrop3FrameAligner : public XdropAligner{

public:
  // Extend an alignment, from the given start point, in the given direction
  // Return the score (bestScore)
  int fillThreeFrame( const uchar* seq1, const uchar* seq2,
		      size_t start1, size_t start2, direction dir,
		      const int sm[MAT][MAT], int maxDrop,
		      const GeneralizedAffineGapCosts& gap,
		      int frameshiftCost, size_t frameSize );

  // Get the ungapped segments of the extension, putting them in "chunks"
  // The extension might be empty!
  void traceThreeFrame( std::vector< SegmentPair >& chunks,
			const uchar* seq1, const uchar* seq2,
			size_t start1, size_t start2, direction dir,
			const int sm[MAT][MAT],
			const GeneralizedAffineGapCosts& gap,
			int frameshiftCost, size_t frameSize );

protected:
  int match2( const uchar* seq1, const uchar* seq2,
	      size_t start1, size_t start2, direction dir,
	      const int sm[MAT][MAT], size_t frameSize,
	      size_t antiDiagonal, size_t seq1pos );

  void reset3( const GeneralizedAffineGapCosts& gap );
  size_t finiteBeg( size_t kn );
  size_t finiteEnd( size_t kn );

  static const uchar* seqPtr2( const uchar* seq2, size_t start2,
			       direction dir, size_t frameSize, size_t pos );

  // get DP matrix value "left of" the given position
  template < typename T >
  T hori3( const std::vector< std::vector< T > >& matrix,
	   size_t antiDiagonal, size_t seq1pos ) const{
    assert( antiDiagonal > 2 );
    if( seq1pos >   fillBeg( antiDiagonal - 3 ) &&
	seq1pos <=  fillEnd( antiDiagonal - 3 ) ){
      return cell( matrix, antiDiagonal - 3, seq1pos - 1 );
    }
    else return -INF;
  }

  // get DP matrix value "above" the given position
  template < typename T >
  T vert3( const std::vector< std::vector< T > >& matrix,
	   size_t antiDiagonal, size_t seq1pos ) const{
    assert( antiDiagonal > 2 );
    if( seq1pos >= fillBeg( antiDiagonal - 3 ) &&
	seq1pos <  fillEnd( antiDiagonal - 3 ) ){
      return cell( matrix, antiDiagonal - 3, seq1pos );
    }
    else return -INF;
  }

  // get DP matrix value "diagonal from" the given position
  template < typename T >
  T diag3( const std::vector< std::vector< T > >& matrix,
	   size_t antiDiagonal, size_t seq1pos ) const{
    if( antiDiagonal > 5 &&
	seq1pos >  fillBeg( antiDiagonal - 6 ) &&
	seq1pos <= fillEnd( antiDiagonal - 6 ) ){
      return cell( matrix, antiDiagonal - 6, seq1pos - 1 );
    }
    else return -INF;
  }
};

template< typename T > T max5( T a, T b, T c, T d, T e ){
  return a > b ? max4(a, c, d, e) : max4(b, c, d, e);
}

template< typename T > int maxIndex5( T a, T b, T c, T d, T e ){
  return e > a ? maxIndex4( b, c, d, e ) + 1 : maxIndex4( a, b, c, d );
}

template< typename T > T min3( T a, T b, T c ){
  return a < b ? std::min(a, c) : std::min(b, c);
}
template< typename T > T min4( T a, T b, T c, T d ){
  return a < b ? min3(a, c, d) : min3(b, c, d);
}

template< typename T > T min5( T a, T b, T c, T d, T e ){
  return a < b ? min4(a, c, d, e) : min4(b, c, d, e);
}

}  // end namespace cbrc
#endif  // XDROP3FRAMEALIGNER_HH
