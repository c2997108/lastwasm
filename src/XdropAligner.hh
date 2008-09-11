// Copyright 2008 Martin C. Frith

// X-drop algorithm for gapped extension of alignments.

#ifndef XDROPALIGNER_HH
#define XDROPALIGNER_HH
#include <vector>
#include <climits>  // INT_MAX

namespace cbrc{

class GeneralizedAffineGapCosts;
class SegmentPair;

class XdropAligner{
public:
  typedef unsigned char uchar;

  enum direction{ FORWARD, REVERSE };

  // Extend an alignment, from the given start point, in the given direction
  // Return the score (bestScore)
  int fill( const uchar* seq1, const uchar* seq2,
	    size_t start1, size_t start2, direction dir,
	    const int sm[64][64], int maxDrop,
	    const GeneralizedAffineGapCosts& gap );

  // Get the ungapped segments of the extension, putting them in "chunks"
  // The extension might be empty!
  void traceback( std::vector< SegmentPair >& chunks,
		  const uchar* seq1, const uchar* seq2,
		  size_t start1, size_t start2, direction dir,
		  const int sm[64][64],
		  const GeneralizedAffineGapCosts& gap ) const;

private:
  typedef std::vector< std::vector< int > > matrix_t;

  enum { INF = INT_MAX / 2 };  // big, but try to avoid overflow

  int bestScore;
  size_t bestAntiDiagonal;
  size_t bestPos1;

  // maybe use std::vector< std::vector< int[3] > > instead of x, y, z???
  matrix_t x;
  matrix_t y;
  matrix_t z;
  std::vector< size_t > offsets;

  static int drop( int score, int minScore );
  static bool isDelimiter( uchar c, const int sm[64][64] );
  static const uchar* seqPtr( const uchar* seq, size_t start,
			      direction dir, size_t pos );
  static int match( const uchar* seq1, const uchar* seq2,
		    size_t start1, size_t start2, direction dir,
		    const int sm[64][64],
		    size_t antiDiagonal, size_t seq1pos );

  size_t fillBeg( size_t antiDiagonal ) const;
  size_t fillEnd( size_t antiDiagonal ) const;
  int cell( const matrix_t& matrix,
	    size_t antiDiagonal, size_t seq1pos ) const;
  int hori( const matrix_t& matrix,
	    size_t antiDiagonal, size_t seq1pos ) const;
  int vert( const matrix_t& matrix,
	    size_t antiDiagonal, size_t seq1pos ) const;
  int diag( const matrix_t& matrix,
	    size_t antiDiagonal, size_t seq1pos ) const;
  void reset( const GeneralizedAffineGapCosts& gap );
  void initScores( size_t antiDiagonal, size_t reservation );
  void updateBest( int score, size_t antiDiagonal, const int* x0 );
  size_t newFillBeg( size_t k1, size_t k2, size_t off1, size_t end1 ) const;
  size_t newFillEnd( size_t k1, size_t k2, size_t off1, size_t end1 ) const;
};

}  // end namespace cbrc
#endif  // XDROPALIGNER_HH
