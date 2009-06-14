// Copyright 2009 Martin C. Frith

#include "SegmentPair.hh"
#include <stdexcept>
//#include <iostream>  // for debugging

using namespace cbrc;

void SegmentPair::makeXdrop( indexT seed1, indexT seed2,
                             const uchar* seq1, const uchar* seq2,
                             const int scoreMatrix[MAT][MAT], int maxDrop,
                             const int pssm2[][MAT] ){
  if( pssm2 ) makeXdropPssm( seed1, seed2, seq1, pssm2, maxDrop );
  else makeXdrop( seed1, seed2, seq1, seq2, scoreMatrix, maxDrop );
}

void SegmentPair::makeXdropPssm( indexT seed1, indexT seed2,
				 const uchar* seq1,
				 const int pssm2[][MAT], int maxDrop ){
  score = 0;
  int scoreDrop = 0;  // drop in score since the maximum score seen so far
  const uchar* s1 = seq1 + seed1;  // starting position in the 1st sequence
  const int (*p2)[MAT] = pssm2 + seed2;  // starting position in the PSSM
  const uchar* bestEnd1 = s1;  // position of max score in the 1st sequence

  for( ;; ){  // extend forwards
    scoreDrop -= ( *p2++ )[ *s1++ ];
    if( scoreDrop < 0 ){
      score -= scoreDrop;  // overflow risk
      scoreDrop = 0;
      bestEnd1 = s1;
    }
    else if( scoreDrop > maxDrop ) break;
  }

  if( score < 0 ) throw std::runtime_error("score overflow detected!");

  scoreDrop = 0;
  s1 = seq1 + seed1;
  p2 = pssm2 + seed2;
  const uchar* bestBeg1 = s1;

  for( ;; ){  // extend backwards
    scoreDrop -= ( *--p2 )[ *--s1 ];
    if( scoreDrop < 0 ){
      score -= scoreDrop;  // overflow risk
      scoreDrop = 0;
      bestBeg1 = s1;
    }
    else if( scoreDrop > maxDrop ) break;
  }

  if( score < 0 ) throw std::runtime_error("score overflow detected!");

  start1 = bestBeg1 - seq1;
  start2 = start1 + seed2 - seed1;
  size = bestEnd1 - bestBeg1;
}
