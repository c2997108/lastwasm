// Copyright 2009 Toshiyuki Sato

#include "Xdrop3FrameAligner.hh"
#include "GeneticCode.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include <algorithm>  // max_element
//#include <iostream>  // for debugging

namespace cbrc{

// get a pointer into a sequence, taking start and direction into account
const Xdrop3FrameAligner::uchar*
Xdrop3FrameAligner::seqPtr2( const uchar* seq2, size_t start2,
			direction dir, size_t frameSize, size_t pos ){
  size_t s2pos;
  size_t startDna = aaToDna( start2, frameSize );

  if( dir == FORWARD ){
    size_t posDna = startDna + pos;
    s2pos = dnaToAa( posDna, frameSize ) - 1;
  }
  else{
    assert( startDna >= pos );
    size_t posDna = startDna - pos;
    s2pos = dnaToAa( posDna, frameSize );
  }

  return &seq2[ s2pos ];
}

// get score for matching 2 residues at the given position
int Xdrop3FrameAligner::match2( const uchar* seq1, const uchar* seq2,
				size_t start1, size_t start2, direction dir,
				const int sm[MAT][MAT], size_t frameSize,
				size_t antiDiagonal, size_t seq1pos ){
  const size_t seq2pos = antiDiagonal - seq1pos * 3;
  if( seq1pos > 0 && seq2pos > 0 ){
    const uchar* s1 = seqPtr( seq1, start1, dir, seq1pos );
    const uchar* s2 = seqPtr2( seq2, start2, dir, frameSize, seq2pos );
    return sm[ *s1 ][ *s2 ];
  }else{
    return -INF;
  }
}

// do the first three cells of the DP matrix
void Xdrop3FrameAligner::reset3( const GeneralizedAffineGapCosts& gap ){
  offsets.resize(3);
  initScores( 0, 1 );
  initScores( 1, 1 );
  initScores( 2, 1 );

  x[0][0] = 0;
  x[1][0] = -INF;
  x[2][0] = -INF;

  y[0][0] = -gap.first;
  y[1][0] = -INF;
  y[2][0] = -INF;

  z[0][0] = -gap.first;
  z[1][0] = -INF;
  z[2][0] = -INF;

  bestScore = 0;
  bestAntiDiagonal = 0;
  bestPos1 = 0;
}

// find start of finite cells on an antidiagonal
size_t Xdrop3FrameAligner::finiteBeg( size_t kn ){
  const int* xn = &x[ kn ][ 0 ];  // front() doesn't work!
  size_t offn = offsets[ kn ];
  size_t endn = offsets[ kn ] + x[ kn ].size();

  for( /* noop */; offn < endn; ++offn ){
    if( *xn++ > -INF ) return offn;
  }
  
  return INF;
}

// find end of finite cells on an antidiagonal
size_t Xdrop3FrameAligner::finiteEnd( size_t kn ){
  const int* xn = &x[ kn ][ 0 ] + x[ kn ].size();  // back() doesn't work!
  size_t offn = offsets[ kn ];
  size_t endn = offsets[ kn ] + x[ kn ].size();

  for( /* noop */; endn > offn; --endn ){
      if( *--xn > -INF ) return endn;
  }

  return 0;
}

// this routine is speed-critical, and it can probably be made faster
int Xdrop3FrameAligner::fillThreeFrame( const uchar* seq1, const uchar* seq2,
			size_t start1, size_t start2, direction dir,
			const int sm[MAT][MAT], int maxDrop,
			const GeneralizedAffineGapCosts& gap,
			int frameshiftCost, size_t frameSize ){

  const int seqIncrement = (dir == FORWARD) ? 1 : -1;
  reset3( gap );
  size_t consecutiveEmptyAntiDiagonals = 0;

  for( size_t k = 3; /* noop */; ++k ){  // loop over antidiagonals
    const size_t k3 = k - 3;
    const size_t k5 = k - 5;  // might wrap around
    const size_t k6 = k - 6;  // might wrap around
    const size_t k7 = k - 7;  // might wrap around

    size_t off3;
    size_t off5 = INF;
    size_t off6 = INF;
    size_t off7 = INF;
    size_t end3;
    size_t end5 = 0;
    size_t end6 = 0;
    size_t end7 = 0;

    size_t moff3;
    size_t moff5 = INF;
    size_t moff6 = INF;
    size_t moff7 = INF;
    size_t mend3;
    size_t mend5 = 0;
    size_t mend6 = 0;
    size_t mend7 = 0;

    off3 = offsets[ k3 ];
    end3 = offsets[ k3 ] + x[ k3 ].size();
    moff3 = finiteBeg( k3 );
    mend3 = finiteEnd( k3 );

    if( k >= 5 ){
      off5 = offsets[ k5 ];
      end5 = offsets[ k5 ] + x[ k5 ].size();
      moff5 = finiteBeg( k5 );
      mend5 = finiteEnd( k5 );

      if( k >= 6 ){
	off6 = offsets[ k6 ];
	end6 = offsets[ k6 ] + x[ k6 ].size();
	moff6 = finiteBeg( k6 );
	mend6 = finiteEnd( k6 );

	if( k >= 7 ){
	  off7 = offsets[ k7 ];
	  end7 = offsets[ k7 ] + x[ k7 ].size();
	  moff7 = finiteBeg( k7 );
	  mend7 = finiteEnd( k7 );
	}
      }
    }

    size_t off0 = min4( moff3,   moff5+1, moff6+1, moff7+1 );	// can be INF
    size_t end0 = max4( mend3+1, mend5+1, mend6+1, mend7+1 );

    if( off0*3 < k-3 && isDelimiter( *seqPtr2( seq2, start2, dir, frameSize, k-off0*3-3 ), sm ) ){
      ++off0;  // don't go past the end of the sequence
    }
    if( off0*3 < k   && isDelimiter( *seqPtr2( seq2, start2, dir, frameSize, k-off0*3   ), sm ) ){
      assert( off0*3 < k );
      ++off0;  // don't go past the end of the sequence
    }
    if( isDelimiter( *seqPtr( seq1, start1, dir, end0-1 ), sm ) ){
      assert( end0 > 1 );
      --end0;  // don't go past the end of the sequence
    }

    if( off0 >= end0 ){
      if( ++consecutiveEmptyAntiDiagonals < 3 ){
	offsets.push_back( 0 );
	initScores( k, 0 );
	continue;
      }
      else{
	break;
      }
    }
    else{
      consecutiveEmptyAntiDiagonals = 0;
    }

    size_t loopBeg = max5( off0, off3+1, off5+1, off6+1, off7+1 );
    size_t loopEnd = min5( end0, end3,   end5+1, end6+1, end7+1 );
    if( loopBeg > loopEnd ){
      loopBeg = off0;
      loopEnd = loopBeg;
    }
    size_t loopSize = loopEnd - loopBeg;

    offsets.push_back( off0 );
    initScores( k, end0 - off0 );
    const int minScore = bestScore - maxDrop;  // set scores < minScore to -INF

    int* x0 = &x[ k ][ 0 ];
    int* y0 = &y[ k ][ 0 ];
    int* z0 = &z[ k ][ 0 ];

    const uchar* s1 = seqPtr( seq1, start1, dir, off0 );
    const uchar* s2 = seqPtr2( seq2, start2, dir, frameSize, k-off0*3 );
    int newPair = hori3( x, k, off0 ) - gap.firstPair;

    // do first cells on boundary
    for( unsigned int i = off0 ; i < loopBeg ; i++ ){
      int yScore3 = hori3( y, k, i );
      int zScore3 = vert3( z, k, i );
      int matrixScore = sm[ *s1 ][ *s2 ];
      int secondScore = diag3( x, k + 1, i ) + matrixScore - frameshiftCost;
      int matchScore  = diag3( x, k    , i ) + matrixScore;
      int fourthScore = diag3( x, k - 1, i ) + matrixScore - frameshiftCost;
      const int xScore = max5( yScore3, zScore3,
			       secondScore, matchScore, fourthScore );
      const int newGap = xScore - gap.first;
      *z0++ = max4( newGap, zScore3 - gap.extend,
		    newPair, diag3( z, k, i ) - gap.extendPair );
      newPair = vert3( x, k, i ) - gap.firstPair;
      *y0++ = max4( newGap, yScore3 - gap.extend,
		    newPair, diag3( y, k, i ) - gap.extendPair );
      updateBest( xScore, k, x0 );
      *x0++ = drop( xScore, minScore );
      s1 += seqIncrement;
      s2 -= seqIncrement;
    }

    if( loopSize > 0 ){
      assert( loopBeg >= 1 + off3 );
      assert( loopBeg >= 1 + off5 );
      assert( loopBeg >= 1 + off6 );
      assert( loopBeg >= 1 + off7 );
      const int* y3 = &y[ k3 ][ loopBeg - 1 - off3 ];
      const int* z3 = &z[ k3 ][ loopBeg - 1 - off3 ];
      const int* x5 = &x[ k5 ][ loopBeg - 1 - off5 ];
      const int* x6 = &x[ k6 ][ loopBeg - 1 - off6 ];
      const int* x7 = &x[ k7 ][ loopBeg - 1 - off7 ];

      const int* const x0end = x0 + loopSize;

      // innermost loop: split into special cases for speed
      if( gap.firstPair >= 2 * gap.first ){  // standard affine gap costs
	do{
	  const int yScore3 = *y3++;
	  const int zScore3 = *++z3;
	  const int secondScore = *x5++ + sm[ *s1 ][ *s2 ] - frameshiftCost;
	  const int matchScore  = *x6++ + sm[ *s1 ][ *s2 ];
	  const int fourthScore = *x7++ + sm[ *s1 ][ *s2 ] - frameshiftCost;
	  const int xScore = max5( yScore3, zScore3,
				   secondScore, matchScore, fourthScore );
	  const int newGap = xScore - gap.first;
	  *z0++ = std::max( newGap, zScore3 - gap.extend );
	  *y0++ = std::max( newGap, yScore3 - gap.extend );
	  updateBest( xScore, k, x0 );
	  *x0++ = drop( xScore, minScore );
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( x0 != x0end );
      }else{  // the general case (generalized affine gap costs)
	const int* x3 = &x[ k3 ][ loopBeg - 1 - off3 ];
	const int* y6 = &y[ k6 ][ loopBeg - 1 - off6 ];
	const int* z6 = &z[ k6 ][ loopBeg - 1 - off6 ];
	newPair = *x3++ - gap.firstPair;
	do{
	  const int yScore3 = *y3++;
	  const int zScore3 = *++z3;
	  const int secondScore = *x5++ + sm[ *s1 ][ *s2 ] - frameshiftCost;
	  const int matchScore  = *x6++ + sm[ *s1 ][ *s2 ];
	  const int fourthScore = *x7++ + sm[ *s1 ][ *s2 ] - frameshiftCost;
	  const int xScore = max5( yScore3, zScore3,
				   secondScore, matchScore, fourthScore );
	  const int newGap = xScore - gap.first;
	  *z0++ = max4( newGap, zScore3 - gap.extend,
			newPair, *z6++ - gap.extendPair );
	  newPair = *x3++ - gap.firstPair;
	  *y0++ = max4( newGap, yScore3 - gap.extend,
			newPair, *y6++ - gap.extendPair );
	  updateBest( xScore, k, x0 );
	  *x0++ = drop( xScore, minScore );
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}while( x0 != x0end );
      }
    }

    newPair = hori3( x, k, loopEnd ) - gap.firstPair;

    // do last cells on boundary
    for( unsigned int i = loopEnd ; i < end0 ; i++ ){
      int yScore3 = hori3( y, k, i );
      int zScore3 = vert3( z, k, i );
      int matrixScore = sm[ *s1 ][ *s2 ];
      int secondScore = diag3( x, k + 1, i ) + matrixScore - frameshiftCost;
      int matchScore  = diag3( x, k    , i ) + matrixScore;
      int fourthScore = diag3( x, k - 1, i ) + matrixScore - frameshiftCost;
      const int xScore = max5( yScore3, zScore3,
			       secondScore, matchScore, fourthScore );
      const int newGap = xScore - gap.first;
      *z0++ = max4( newGap, zScore3 - gap.extend,
		    newPair, diag3( z, k, i ) - gap.extendPair );
      newPair = vert3( x, k, i ) - gap.firstPair;
      *y0++ = max4( newGap, yScore3 - gap.extend,
		    newPair, diag3( y, k, i ) - gap.extendPair );
      updateBest( xScore, k, x0 );
      *x0++ = drop( xScore, minScore );
      s1 += seqIncrement;
      s2 -= seqIncrement;
    }
  }
 
  return bestScore;
}

//
void Xdrop3FrameAligner::traceThreeFrame( std::vector< SegmentPair >& chunks,
			      const uchar* seq1, const uchar* seq2,
 			      size_t start1, size_t start2, direction dir,
			      const int sm[MAT][MAT],
			      const GeneralizedAffineGapCosts& gap,
			      int frameshiftCost, size_t frameSize ){
  size_t k = bestAntiDiagonal;
  size_t i = bestPos1;
  size_t oldPos1 = i;
  int state = 0;  // enum? can be: 0 = match(x), 1 = gap(y), or 2 = gap(z)

  while( k > 2 ){
    if( state == 0 ){
      const int matrixScore = match2( seq1, seq2, start1, start2,
				      dir, sm, frameSize, k, i );
      const int shiftedScore = matrixScore - frameshiftCost;
      const int m = maxIndex5( diag3( x, k,     i ) + matrixScore,
			       hori3( y, k,     i ),
			       vert3( z, k,     i ),
			       diag3( x, k + 1, i ) + shiftedScore,
			       diag3( x, k - 1, i ) + shiftedScore );
      if( m < 1 || m > 2 ){
	k -= 6;
	i -= 1;
      }

      if( (m > 0 && oldPos1 != i) || k < 3 ){
	chunks.push_back( SegmentPair( i, k - i * 3, oldPos1 - i ) );
      }

      if( m == 1 ){
	k -= 3;
	i -= 1;
      }
      else if( m == 2 ){
	k -= 3;
      }
      else if( m == 3 ){
	k += 1;
	oldPos1 = i;
      }
      else if( m == 4 ){
	k -= 1;
	oldPos1 = i;
      }

      if( m < 3 ){
	state = m;
      }
      else {
	state = 0;
      }
    }
    else if( state == 1 ){
      const int m = maxIndex4( cell ( x, k, i ) - gap.first,
			       hori3( y, k, i ) - gap.extend,
			       vert3( x, k, i ) - gap.firstPair,
			       diag3( y, k, i ) - gap.extendPair );
      if( m == 1 ){
	k -= 3;
      }
      else if( m == 2 ){
	k -= 3;
      }
      else if( m == 3 ){
	k -= 6;
      }

      i -= m % 2;
      state = m % 2;
      oldPos1 = i;
    }
    else if( state == 2 ){
      const int m = maxIndex4( cell ( x, k, i ) - gap.first,
			       vert3( z, k, i ) - gap.extend,
			       hori3( x, k, i ) - gap.firstPair,
			       diag3( z, k, i ) - gap.extendPair );
      if( m == 1 ){
	k -= 3;
      }
      else if( m == 2 ){
	k -= 3;
      }
      else if( m == 3 ){
	k -= 6;
      }

      i -= (m > 1);
      state = (m % 2) * 2;
      oldPos1 = i;
    }
  }
}

}  // end namespace cbrc
