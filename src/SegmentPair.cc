// Copyright 2008, 2009, 2011 Martin C. Frith

#include "SegmentPair.hh"

namespace cbrc{

void SegmentPair::maxIdenticalRun( const uchar* seq1, const uchar* seq2,
				   const uchar* canonical ){
  const uchar* s1 = seq1 + beg1();
  const uchar* s2 = seq2 + beg2();
  const uchar* e1 = seq1 + end1();
  const uchar* runBeg1 = s1;
  size = 0;

  while( s1 < e1 ){
    if( canonical[ *s1++ ] == canonical[ *s2++ ] ){
      if( indexT(s1 - runBeg1) > size ){
	start1 = runBeg1 - seq1;
	start2 = (s2 - seq2) - (s1 - runBeg1);  // move this out of the loop?
	size = s1 - runBeg1;
      }
    }
    else runBeg1 = s1;
  }
}

}
