// Copyright 2008 Martin C. Frith

#include "SegmentPair.hh"
#include <stdexcept>

namespace cbrc{

void SegmentPair::makeXdrop( indexT seed1, indexT seed2,
			     const uchar* seq1, const uchar* seq2,
			     const int scoreMatrix[MAT][MAT], int maxDrop ){
  score = 0;
  int scoreDrop = 0;  // drop in score since the maximum score seen so far
  const uchar* s1 = seq1 + seed1;  // starting position in the 1st sequence
  const uchar* s2 = seq2 + seed2;  // starting position in the 2nd sequence
  const uchar* bestEnd1 = s1;  // position of max score in the 1st sequence

  for( ;; ){  // extend forwards
    scoreDrop -= scoreMatrix[ *s1++ ][ *s2++ ];
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
  s2 = seq2 + seed2;
  const uchar* bestBeg1 = s1;

  for( ;; ){  // extend backwards
    scoreDrop -= scoreMatrix[ *--s1 ][ *--s2 ];
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

bool SegmentPair::isOptimal( const uchar* seq1, const uchar* seq2,
			     const int scoreMatrix[MAT][MAT],
			     int maxDrop ) const{
  int maxScore = 0;
  int runningScore = 0;
  const uchar* s1 = seq1 + start1;
  const uchar* s2 = seq2 + start2;
  const uchar* s1end = s1 + size;

  while( s1 < s1end ){
    runningScore += scoreMatrix[ *s1++ ][ *s2++ ];
    if( runningScore > maxScore ) maxScore = runningScore;
    else if( runningScore <= 0 ||                  // non-optimal prefix
	     s1 == s1end ||                        // non-optimal suffix
	     runningScore < maxScore - maxDrop ){  // excessive score drop
      return false;
    }
  }
  return true;
}

void SegmentPair::maxIdenticalRun( const uchar* seq1, const uchar* seq2,
				   const uchar* canonical,
				   const int scoreMatrix[MAT][MAT] ){
  const uchar* s1 = seq1 + start1;
  const uchar* s2 = seq2 + start2;
  const uchar* s1end = s1 + size;
  const uchar* runBeg1 = s1;
  int runningScore = 0;
  size = 0;

  while( s1 < s1end ){
    runningScore += scoreMatrix[ *s1 ][ *s2 ];
    if( canonical[ *s1++ ] == canonical[ *s2++ ] ){
      if( indexT(s1 - runBeg1) > size ){
	start1 = runBeg1 - seq1;
	start2 = (s2 - seq2) - (s1 - runBeg1);  // move this out of the loop?
	size = s1 - runBeg1;
	score = runningScore;
      }
    }
    else{
      runBeg1 = s1;
      runningScore = 0;
    }
  }
}

}  // end namespace cbrc
