// Copyright 2009 Martin C. Frith

#include "XdropAligner.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include <algorithm>
//#include <iostream>  // for debugging

using namespace cbrc;

int XdropAligner::fill( const uchar* seq1, const uchar* seq2,
			size_t start1, size_t start2, direction dir,
			const int sm[MAT][MAT], int smMax, int maxDrop,
			const GeneralizedAffineGapCosts& gap,
			const int pssm2[][MAT] ){
  if( pssm2 ) return fillPssm( seq1, seq2, start1, start2, dir,
			       sm, smMax, maxDrop, gap, pssm2 );
  else return fill( seq1, seq2, start1, start2, dir, sm, smMax, maxDrop, gap );
}

// tried to make this fast, using low-level pointer operations
// uses the DP rearrangement from M Cameron, HE Williams, A Cannane 2004
int XdropAligner::fillPssm( const uchar* seq1, const uchar* seq2,
			    size_t start1, size_t start2, direction dir,
			    const int sm[MAT][MAT], int smMax, int maxDrop,
			    const GeneralizedAffineGapCosts& gap,
			    const int pssm2[][MAT] ){

  const int seqIncrement = (dir == FORWARD) ? 1 : -1;

  reset( gap );

  for( size_t k = 1; /* noop */; ++k ){  // loop over antidiagonals
    const size_t k1 = k - 1;
    const size_t k2 = k - 2;  // might wrap around
    const size_t off1 = offsets[ k1 ];
    const size_t end1 = off1 + x[ k1 ].size();
    /* */ size_t off0 = newFillBeg( k1, k2, off1, end1 );
    /* */ size_t end0 = newFillEnd( k1, k2, off1, end1 );
    const size_t loopBeg = off0 + ( off0 == off1 );
    const size_t loopEnd = end0 - ( end0 > end1 );
    const size_t loopSize = loopEnd - loopBeg;

    if( off0 < k && isDelimiter( *seqPtr( seq2, start2, dir, k-off0 ), sm ) ){
      ++off0;  // don't go past the end of the sequence
      // speedup: don't let the score drop by max-matches * max-match-score
      maxDrop = std::min( maxDrop, int(loopSize) * smMax - 1 );
    }

    if( end0 > 1 && isDelimiter( *seqPtr( seq1, start1, dir, end0-1 ), sm ) ){
      --end0;  // don't go past the end of the sequence
      // speedup: don't let the score drop by max-matches * max-match-score
      maxDrop = std::min( maxDrop, int(loopSize) * smMax - 1 );
    }

    if( off0 >= end0 )  break;

    offsets.push_back( off0 );
    initScores( k, end0 - off0 );
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

    if( loopSize > 0 ){
      //assert( k > 1 );  // without this, it's much faster (!?!?)
      const int* const x0end = x0 + loopSize;
      const uchar* s1 = seqPtr( seq1, start1, dir, loopBeg );
      const int (*p2)[MAT] = seqPtr( pssm2, start2, dir, k - loopBeg );
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
	  const int matchScore = *x2++ + ( *p2 )[ *s1 ];
	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	  const int yScore1 = *y1++;
	  const int zScore1 = *++z1;
	  const int xScore = max3( matchScore, yScore1, zScore1 );
	  if( xScore >= minScore ){
	    const int newGap = xScore - gap.first;
	    *z0 = std::max( newGap, zScore1 - gap.extend );
	    *y0 = std::max( newGap, yScore1 - gap.extend );
	    updateBest( xScore, k, x0 );
	    *x0 = xScore;
	  }
	  else *x0 = *y0 = *z0 = -INF;
	  ++x0;  ++y0;  ++z0;
	}while( x0 != x0end );
      }else{  // the general case (generalized affine gap costs)
	const int* x1 = &x[ k1 ][ horiBeg ];
	const int* y2 = &y[ k2 ][ diagBeg ];
	const int* z2 = &z[ k2 ][ diagBeg ];
	int newPair = *x1++ - gap.firstPair;
	do{
	  const int matchScore = *x2++ + ( *p2 )[ *s1 ];
	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	  const int yScore1 = *y1++;
	  const int zScore1 = *++z1;
	  const int xScore = max3( matchScore, yScore1, zScore1 );
	  const int newGap = xScore - gap.first;
	  *z0++ = max4( newGap, zScore1 - gap.extend,
			newPair, *z2++ - gap.extendPair );
	  newPair = *x1++ - gap.firstPair;
	  *y0++ = max4( newGap, yScore1 - gap.extend,
			newPair, *y2++ - gap.extendPair );
	  updateBest( xScore, k, x0 );
	  *x0++ = drop( xScore, minScore );
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
