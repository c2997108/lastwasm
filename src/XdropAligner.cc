// Copyright 2008 Martin C. Frith

#include "XdropAligner.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include <algorithm>  // max_element
#include <cassert>
//#include <iostream>  // for debugging

namespace{

int max3( int a, int b, int c ){
  return a > b ? std::max(a, c) : std::max(b, c);
}

int max4( int a, int b, int c, int d ){
  return a > b ? max3(a, c, d) : max3(b, c, d);
}

int maxIndex2( int a, int b ){
  return b > a ? 1 : 0;
}

int maxIndex3( int a, int b, int c ){
  return c > a ? maxIndex2( b, c ) + 1 : maxIndex2( a, b );
}

int maxIndex4( int a, int b, int c, int d ){
  return d > a ? maxIndex3( b, c, d ) + 1 : maxIndex3( a, b, c );
}

}

namespace cbrc{

// apply the X-drop to a score
int XdropAligner::drop( int score, int minScore ){
  return score < minScore ? -INF : score;
}

bool XdropAligner::isDelimiter( uchar c, const int sm[64][64] ){
  return sm[c][c] == -INF;  // assumes it's the same INF!
}

// get a pointer into a sequence, taking start and direction into account
const XdropAligner::uchar*
XdropAligner::seqPtr( const uchar* seq, size_t start,
		      direction dir, size_t pos ){
  assert( pos > 0 );
  if( dir == FORWARD ){
    return &seq[ start + pos - 1 ];
  }else{
    return &seq[ start - pos ];
  }
}

// get score for matching 2 residues at the given position
int XdropAligner::match( const uchar* seq1, const uchar* seq2,
			 size_t start1, size_t start2, direction dir,
			 const int sm[64][64],
			 size_t antiDiagonal, size_t seq1pos ){
  const size_t seq2pos = antiDiagonal - seq1pos;
  if( seq1pos > 0 && seq2pos > 0 ){
    return sm[ *seqPtr( seq1, start1, dir, seq1pos ) ]
             [ *seqPtr( seq2, start2, dir, seq2pos ) ];
  }else{
    return -INF;
  }
}

// get start of filled cells on an antidiagonal
size_t XdropAligner::fillBeg( size_t antiDiagonal ) const{
  return offsets[ antiDiagonal ];
}

// get end of filled cells on an antidiagonal
size_t XdropAligner::fillEnd( size_t antiDiagonal ) const{
  return offsets[ antiDiagonal ] + x[ antiDiagonal ].size();
}

// get DP matrix value at the given position, taking offsets into account
int XdropAligner::cell( const matrix_t& matrix,
			size_t antiDiagonal, size_t seq1pos ) const{
  assert( seq1pos >= fillBeg( antiDiagonal ) );
  assert( seq1pos < fillEnd( antiDiagonal ) );
  return matrix[ antiDiagonal ][ seq1pos - offsets[ antiDiagonal ] ];
}

// get DP matrix value "left of" the given position
int XdropAligner::hori( const matrix_t& matrix,
			size_t antiDiagonal, size_t seq1pos ) const{
  assert( antiDiagonal > 0 );
  if( seq1pos > fillBeg( antiDiagonal-1 ) ){
    return cell( matrix, antiDiagonal-1, seq1pos-1 );
  }
  else{
    return -INF;
  }
}

// get DP matrix value "above" the given position
int XdropAligner::vert( const matrix_t& matrix,
			size_t antiDiagonal, size_t seq1pos ) const{
  assert( antiDiagonal > 0 );
  if( seq1pos < fillEnd( antiDiagonal-1 ) ){
    return cell( matrix, antiDiagonal-1, seq1pos );
  }
  else{
    return -INF;
  }
}

// get DP matrix value "diagonal from" the given position
int XdropAligner::diag( const matrix_t& matrix,
			size_t antiDiagonal, size_t seq1pos ) const{
  if( antiDiagonal > 1 &&
      seq1pos > fillBeg( antiDiagonal-2 ) &&
      seq1pos <= fillEnd( antiDiagonal-2 ) ){
    return cell( matrix, antiDiagonal-2, seq1pos-1 );
  }else{
    return -INF;
  }
}

// do the first cell (first antidiagonal) of the DP matrix
void XdropAligner::reset( const GeneralizedAffineGapCosts& gap ){
  offsets.resize(1);
  initScores( 0, 1 );
  x[0][0] = 0;
  y[0][0] = -gap.first;
  z[0][0] = -gap.first;
  bestScore = 0;
  bestAntiDiagonal = 0;
  bestPos1 = 0;
}

// allocate memory for an antidiagonal
void XdropAligner::initScores( size_t antiDiagonal, size_t reservation ){
  if( antiDiagonal < x.size() ){
    x[ antiDiagonal ].resize( reservation );
    y[ antiDiagonal ].resize( reservation );
    z[ antiDiagonal ].resize( reservation );
  }
  else{
    assert( antiDiagonal == x.size() );
    x.push_back( std::vector< int >( reservation ) );
    y.push_back( std::vector< int >( reservation ) );
    z.push_back( std::vector< int >( reservation ) );
  }
}

void XdropAligner::updateBest( int score, size_t antiDiagonal, const int* x0 ){
  if( score > bestScore ){
    bestScore = score;
    bestAntiDiagonal = antiDiagonal;
    bestPos1 = x0 - &x[ antiDiagonal ][ 0 ] + offsets[ antiDiagonal ];
  }
}

// find start of cells that need to be filled on an antidiagonal
size_t XdropAligner::newFillBeg( size_t k1, size_t k2,
				 size_t off1, size_t end1 ) const{
  const int* x1 = &x[ k1 ][ 0 ];  // front() doesn't work!
  size_t seq1pos = off1;

  if( *x1++ == -INF ){
    const int* x2 = &x[ k2 ][ 0 ] + ( off1 - offsets[ k2 ] );
    for( ++seq1pos; seq1pos < end1; ++seq1pos ){
      if( *x1++ > -INF || *x2++ > -INF )  break;
    }
  }

  return seq1pos;
}

// find end of cells that need to be filled on an antidiagonal
size_t XdropAligner::newFillEnd( size_t k1, size_t k2,
				 size_t off1, size_t end1 ) const{
  const int* x1 = &x[ k1 ][ 0 ] + x[ k1 ].size();  // back() doesn't work!
  size_t seq1pos = end1;

  if( *--x1 == -INF ){
    const int* x2 = &x[ k2 ][ 0 ] + ( end1 - 1 - offsets[ k2 ] );
    for( --seq1pos; seq1pos > off1; --seq1pos ){
      if( *--x1 > -INF || *--x2 > -INF )  break;
    }
  }

  return seq1pos + 1;  // return one-past-the-end
}

// tried to make this fast, using low-level pointer operations
// uses the DP rearrangement from M Cameron, HE Williams, A Cannane 2004
int XdropAligner::fill( const uchar* seq1, const uchar* seq2,
			size_t start1, size_t start2, direction dir,
			const int sm[64][64], int maxDrop,
			const GeneralizedAffineGapCosts& gap ){

  const int seqIncrement = (dir == FORWARD) ? 1 : -1;

  reset( gap );

  for( size_t k = 1; /* noop */; ++k ){  // loop over antidiagonals
    const size_t k1 = k - 1;
    const size_t k2 = k - 2;  // might wrap around
    const size_t off1 = offsets[ k1 ];
    const size_t end1 = off1 + x[ k1 ].size();
    /* */ size_t off0 = newFillBeg( k1, k2, off1, end1 );
    /* */ size_t end0 = newFillEnd( k1, k2, off1, end1 );
    if( off0 < k && isDelimiter( *seqPtr( seq2, start2, dir, k-off0 ), sm ) ){
      ++off0;
    }
    if( end0 > 1 && isDelimiter( *seqPtr( seq1, start1, dir, end0-1 ), sm ) ){
      --end0;
    }
    if( off0 >= end0 )  break;
    offsets.push_back( off0 );
    initScores( k, end0 - off0 );
    const size_t loopBeg = off0 + ( off0 == off1 );
    const size_t loopEnd = end0 - ( end0 > end1 );
    const int minScore = bestScore - maxDrop;  // set scores < minScore to -INF

    int* x0 = &x[ k ][ 0 ];
    int* y0 = &y[ k ][ 0 ];
    int* z0 = &z[ k ][ 0 ];

    if( off0 == off1 ){  // do first cell on boundary
      const int xScore = z[ k1 ].front();
      *x0++ = drop( xScore, minScore );
      *y0++ = max3( xScore - gap.first,
		    diag( y, k, off1 ) - gap.extendPair,  // necessary???
		    x[ k1 ].front() - gap.firstPair );
      *z0++ = xScore - gap.extend;
    }

    if( loopBeg < loopEnd ){
      assert( k > 1 );
      const int* const x0end = x0 + (loopEnd - loopBeg);
      const uchar* s1 = seqPtr( seq1, start1, dir, loopBeg );
      const uchar* s2 = seqPtr( seq2, start2, dir, k - loopBeg );
      const size_t horiBeg = loopBeg - 1 - off1;
      const size_t diagBeg = loopBeg - 1 - offsets[ k2 ];
      const int* y1 = &y[ k1 ][ horiBeg ];
      const int* z1 = &z[ k1 ][ horiBeg ];
      const int* x2 = &x[ k2 ][ diagBeg ];

      // innermost loop: split into special cases for speed
      if( gap.firstPair >= 2 * gap.first ){  // standard affine gap costs
	do{
	  // if xScore <- yScore1, then newGap <= yScore1 - gap.extend:
	  // but using this logic doesn't seem to make it faster!
	  const int yScore1 = *y1++;
	  const int zScore1 = *++z1;
	  const int matchScore = *x2++ + sm[ *s1 ][ *s2 ];
	  const int xScore = max3( matchScore, yScore1, zScore1 );
	  const int newGap = xScore - gap.first;
	  *z0++ = std::max( newGap, zScore1 - gap.extend );
	  *y0++ = std::max( newGap, yScore1 - gap.extend );
	  updateBest( xScore, k, x0 );
	  *x0++ = drop( xScore, minScore );
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( x0 != x0end );
      }else{  // the general case (generalized affine gap costs)
	const int* x1 = &x[ k1 ][ horiBeg ];
	const int* y2 = &y[ k2 ][ diagBeg ];
	const int* z2 = &z[ k2 ][ diagBeg ];
	int newPair = *x1++ - gap.firstPair;
	do{
	  const int yScore1 = *y1++;
	  const int zScore1 = *++z1;
	  const int matchScore = *x2++ + sm[ *s1 ][ *s2 ];
	  const int xScore = max3( matchScore, yScore1, zScore1 );
	  const int newGap = xScore - gap.first;
	  *z0++ = max4( newGap, zScore1 - gap.extend,
			newPair, *z2++ - gap.extendPair );
	  newPair = *x1++ - gap.firstPair;
	  *y0++ = max4( newGap, yScore1 - gap.extend,
			newPair, *y2++ - gap.extendPair );
	  updateBest( xScore, k, x0 );
	  *x0++ = drop( xScore, minScore );
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( x0 != x0end );
      }
    }

    if( end0 > end1 ){  // do last cell on boundary
      const int xScore = y[ k1 ].back();
      *x0++ = drop( xScore, minScore );
      *y0++ = xScore - gap.extend;
      *z0++ = max3( xScore - gap.first,
		    diag( z, k, end1 ) - gap.extendPair,  // necessary???
		    x[ k1 ].back() - gap.firstPair );
    }
  }

  return bestScore;
}

void XdropAligner::traceback( std::vector< SegmentPair >& chunks,
			      const uchar* seq1, const uchar* seq2,
			      size_t start1, size_t start2, direction dir,
			      const int sm[64][64],
			      const GeneralizedAffineGapCosts& gap ) const{
  size_t k = bestAntiDiagonal;
  size_t i = bestPos1;
  size_t oldPos1 = i;
  int state = 0;  // enum? can be: 0 = match(x), 1 = gap(y), or 2 = gap(z)

  while( k > 0 ){
    if( state == 0 ){
      const int m =
	maxIndex3( diag( x, k, i ) +
		   match( seq1, seq2, start1, start2, dir, sm, k, i ),
		   hori( y, k, i ),
		   vert( z, k, i ) );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }

      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }

      if( m > 0 ){
	k -= 1;
	i -= (m < 2);
	state = m;
      }
    }
    else if( state == 1 ){
      const int m = maxIndex4( cell( x, k, i ) - gap.first,
			       hori( y, k, i ) - gap.extend,
			       vert( x, k, i ) - gap.firstPair,
			       diag( y, k, i ) - gap.extendPair );
      k -= (m + 1) / 2;
      i -= m % 2;
      state = m % 2;
      oldPos1 = i;
    }
    else{
      const int m = maxIndex4( cell( x, k, i ) - gap.first,
			       vert( z, k, i ) - gap.extend,
			       hori( x, k, i ) - gap.firstPair,
			       diag( z, k, i ) - gap.extendPair );
      k -= (m + 1) / 2;
      i -= (m > 1);
      state = (m % 2) * 2;
      oldPos1 = i;
    }
  }
}

}  // end namespace cbrc
