// Copyright 2008, 2009 Martin C. Frith

#include "SuffixArray.hh"
#include "io.hh"
#include <cassert>
//#include <iostream>  // for debugging

namespace cbrc{

int SuffixArray::makeIndex( indexT beg, indexT end, indexT step,
			    std::size_t maxBytes ){
  assert( step > 0 );
  indexT oldSize = index.size();

  for( indexT i = beg; i < end; i += step ){
    if( text[i] < alphSize ){
      index.push_back(i);

      if( indexBytes() > maxBytes ){
	index.erase( index.begin() + oldSize, index.end() );
	return 0;
      }
    }
  }

  return 1;
}

void SuffixArray::sortIndex(){
  if( alphSize == 4 )
    radixSortAlph4( &index[0], &index[0] + index.size(), 0, &text[0] );
  else
    radixSort( &index[0], &index[0] + index.size(), 0, &text[0] );
}

void SuffixArray::clear(){
  index.clear();
  buckets.clear();
  bucketSteps.clear();
  bucketMask.clear();
}

void SuffixArray::fromFiles( const std::string& baseName,
			     indexT indexNum, indexT bucketDepth ){
  index.resize(indexNum);  // unwanted zero-fill
  vectorFromBinaryFile( index, baseName + ".suf" );

  makeBucketSteps(bucketDepth);
  makeBucketMask(bucketDepth);

  buckets.resize( bucketSteps[0] * alphSize + 1 );  // unwanted zero-fill
  vectorFromBinaryFile( buckets, baseName + ".bck" );
}

void SuffixArray::toFiles( const std::string& baseName ) const{
  vectorToBinaryFile( index, baseName + ".suf" );
  vectorToBinaryFile( buckets, baseName + ".bck" );
}

// use past results to speed up long matches?
// could & probably should return the match depth
void SuffixArray::match( const indexT*& beg, const indexT*& end,
			 const uchar* queryPtr,
			 indexT maxHits, indexT minDepth ) const{
  // match using buckets:
  indexT bucketDepth = bucketMask.size();
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();
  const indexT* bucketPtr = &buckets[0];

  for( indexT depth = 0; depth < bucketDepth; ++depth ){
    uchar symbol = *queryPtr;
    if( symbol < alphSize ){
      indexT step = bucketSteps[depth];
      bucketPtr += symbol * step;
      bucketBeg = *bucketPtr;
      bucketEnd = *(bucketPtr + step);
    }else{
      bucketBeg = bucketEnd;
      depth = minDepth;
    }
    if( bucketEnd - bucketBeg <= maxHits && depth+1 >= minDepth ){
      beg = &index[0] + bucketBeg;
      end = &index[0] + bucketEnd;
      return;
    }
    queryPtr += bucketMask[depth];
  }

  // match using binary search:
  beg = &index[0] + bucketBeg;
  end = &index[0] + bucketEnd;
  const uchar* textBase = &text[0] + bucketMaskTotal;

  for( indexT depth = bucketDepth; /* noop */; ++depth ){
    uchar symbol = *queryPtr;
    if( symbol < alphSize ){
      equalRange( beg, end, textBase, symbol );
    }else{
      beg = end;
      depth = minDepth;
    }
    if( indexT(end - beg) <= maxHits && depth+1 >= minDepth ) return;
    indexT off = mask[ depth % maskSize ];
    queryPtr += off;
    textBase += off;
  }
}

void SuffixArray::countMatches( std::vector<unsigned long long>& counts,
				const uchar* queryPtr ) const{
  // match using buckets:
  indexT bucketDepth = bucketMask.size();
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();
  const indexT* bucketPtr = &buckets[0];

  for( indexT depth = 0; depth < bucketDepth; ++depth ){
    uchar symbol = *queryPtr;
    if( symbol >= alphSize ) return;
    indexT step = bucketSteps[depth];
    bucketPtr += symbol * step;
    bucketBeg = *bucketPtr;
    bucketEnd = *(bucketPtr + step);
    if( bucketBeg == bucketEnd ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += bucketEnd - bucketBeg;
    queryPtr += bucketMask[depth];
  }

  // match using binary search:
  const indexT* beg = &index[0] + bucketBeg;
  const indexT* end = &index[0] + bucketEnd;
  const uchar* textBase = &text[0] + bucketMaskTotal;

  for( indexT depth = bucketDepth; /* noop */; ++depth ){
    uchar symbol = *queryPtr;
    if( symbol >= alphSize ) return;
    equalRange( beg, end, textBase, symbol );
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    indexT off = mask[ depth % maskSize ];
    queryPtr += off;
    textBase += off;
  }
}

void SuffixArray::equalRange( const indexT*& beg, const indexT*& end,
			      const uchar* textBase, uchar symbol ){
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    uchar s = textBase[ *mid ];
    if( s < symbol ){
      beg = mid + 1;
    }else if( s > symbol ){
      end = mid;
    }else{
      beg = lowerBound( beg, mid, textBase, symbol );
      end = upperBound( mid + 1, end, textBase, symbol );
      return;
    }
  }
}

const SuffixArray::indexT*
SuffixArray::lowerBound( const indexT* beg, const indexT* end,
			 const uchar* textBase, uchar symbol ){
  for( ;; ){
    std::size_t size = end - beg;
    if( size <= 4 ) break;  // 3,4 seem good for hg18 chr21 versus itself
    const indexT* mid = beg + size / 2;
    if( textBase[ *mid ] < symbol ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }

  while( textBase[ *beg ] < symbol ) ++beg;  // linear search

  return beg;
}

const SuffixArray::indexT*
SuffixArray::upperBound( const indexT* beg, const indexT* end,
			 const uchar* textBase, uchar symbol ){
  for( ;; ){
    std::size_t size = end - beg;
    if( size <= 4 ) break;  // 3,4 seem good for hg18 chr21 versus itself
    const indexT* mid = beg + size / 2;
    if( textBase[ *mid ] <= symbol ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }

  while( textBase[ *(end-1) ] > symbol ) --end;  // linear search

  return end;
}

void SuffixArray::makeBuckets( indexT bucketDepth ){
  makeBucketSteps(bucketDepth);
  makeBucketMask(bucketDepth);

  for( indexT i = 0; i < index.size(); ++i ){
    indexT bucketIndex = 0;
    const uchar* textPtr = &text[ index[i] ];

    for( indexT d = 0; d < bucketDepth; ++d ){
      unsigned symbol = *textPtr;
      indexT step = bucketSteps[d];
      bucketIndex += symbol * step;
      if( symbol >= alphSize ) break;
      textPtr += bucketMask[d];
    }

    buckets.resize( bucketIndex+1, i );
  }

  buckets.resize( bucketSteps[0] * alphSize + 1, index.size() );
}

void SuffixArray::makeBucketSteps( indexT bucketDepth ){
  bucketSteps.resize(bucketDepth);
  indexT depth = bucketDepth;
  indexT step = 1;
  while( depth > 0 ){
    --depth;
    bucketSteps[depth] = step;
    step = step * alphSize + 1;
  }
}

void SuffixArray::makeBucketMask( indexT bucketDepth ){
  bucketMaskTotal = 0;  // necessary!!!
  bucketMask.resize(bucketDepth);

  for( indexT i = 0; i < bucketDepth; ++i ){
    indexT offset = mask[ i % maskSize ];
    bucketMask[i] = offset;
    bucketMaskTotal += offset;
  }
}

// This radix sort is adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#define PUSH(b1, e1, d1, t1) sp->b = b1, sp->e = e1, sp->d = d1, (sp++)->t = t1
#define  POP(b1, e1, d1, t1) b1 = (--sp)->b, e1 = sp->e, d1 = sp->d, t1 = sp->t

static SuffixArray::Stack stack[1048576];  // big enough???

void SuffixArray::radixSort( indexT* beg, indexT* end,
			     indexT depth, const uchar* textBase ){
  static indexT bucketSize[256];  // initialized to zero at startup
  /*  */ indexT* bucketEnd[256];  // "static" makes little difference to speed
  Stack* sp = stack;
  PUSH( beg, end, depth, textBase );

  while( sp > stack ){
    POP( beg, end, depth, textBase );

    if( end - beg < 10 ){  // 10 seems good for hg18 chr21
      insertionSort( beg, end, depth, textBase );
      continue;
    }

    // get bucket sizes (i.e. letter counts):
    // The intermediate oracle array makes it faster (see "Engineering
    // Radix Sort for Strings" by J Karkkainen & T Rantala)
    for( indexT* i = beg; i < end; /* noop */ ){
      uchar oracle[256];
      uchar* oracleEnd =
	oracle + std::min( sizeof(oracle), std::size_t(end - i) );
      for( uchar* j = oracle; j < oracleEnd; ++j )
	*j = textBase[ *i++ ];
      for( uchar* j = oracle; j < oracleEnd; ++j )
	++bucketSize[ *j ];
    }

    // get the next textBase, and increment depth, for sorting within buckets:
    const uchar* textNext = textBase + mask[ depth % maskSize ];
    ++depth;

    // get bucket ends, and put buckets on the stack to sort within them later:
    // (could push biggest bucket first, to ensure logarithmic stack growth)
    indexT* pos = beg;
    for( unsigned i = 0; i < alphSize; ++i ){
      indexT* nextPos = pos + bucketSize[i];
      PUSH( pos, nextPos, depth, textNext );
      pos = nextPos;
      bucketEnd[i] = pos;
    }
    bucketEnd[alphSize] = end;  // don't sort within the delimiter bucket

    // permute items into the correct buckets:
    for( indexT* i = beg; i < end; /* noop */ ) {
      unsigned symbol;  // unsigned is faster than uchar!
      indexT holdOut = *i;
      while( --bucketEnd[ symbol = textBase[holdOut] ] > i ){
	std::swap( *bucketEnd[symbol], holdOut );
      }
      *i = holdOut;
      i += bucketSize[symbol];
      bucketSize[symbol] = 0;  // reset it so we can reuse it
    }
  }
}

// Specialized radix sort for alphSize = 4 (plus one delimiter symbol)
void SuffixArray::radixSortAlph4( indexT* beg, indexT* end,
				  indexT depth, const uchar* textBase ){
  Stack* sp = stack;
  PUSH( beg, end, depth, textBase );

  while( sp > stack ){
    POP( beg, end, depth, textBase );

    if( end - beg < 10 ){
      insertionSort( beg, end, depth, textBase );
      continue;
    }

    indexT* end0 = beg;  // end of '0's
    indexT* end1 = beg;  // end of '1's
    indexT* end2 = beg;  // end of '2's
    indexT* beg3 = end;  // beginning of '3's
    indexT* beg4 = end;  // beginning of '4's (delimiters)

    while( end2 < beg3 ){
      const indexT x = *end2;
      switch( textBase[x] ){
      case 0:
	*end2++ = *end1;
	*end1++ = *end0;
	*end0++ = x;
	break;
      case 1:
	*end2++ = *end1;
	*end1++ = x;
	break;
      case 2:
	end2++;
	break;
      case 3:
	*end2 = *--beg3;
	*beg3 = x;
	break;
      case 4:
	*end2 = *--beg3;
	*beg3 = *--beg4;
	*beg4 = x;
	break;
      }
    }

    // get the next textBase, and increment depth, for sorting within buckets:
    const uchar* textNext = textBase + mask[ depth % maskSize ];
    ++depth;

    // put buckets on the stack to sort within them later:
    // (could push biggest bucket first, to ensure logarithmic stack growth)
    PUSH( beg, end0, depth, textNext );   // the '0's
    PUSH( end0, end1, depth, textNext );  // the '1's
    PUSH( end1, beg3, depth, textNext );  // the '2's
    PUSH( beg3, beg4, depth, textNext );  // the '3's
    // don't sort within the delimiter bucket
  }
}

void SuffixArray::insertionSort( indexT* beg, indexT* end,
				 indexT depth, const uchar* textBase ){
  if( maskSize == 1 && mask[0] == 1 ){
    insertionSortSimple( beg, end, textBase );  // how much faster?
    return;
  }

  const indexT* maskPtr = &mask[ depth % maskSize ];
  const indexT* maskEnd = &mask[0] + maskSize;

  for( indexT* i = beg+1; i < end; ++i ){
    for( indexT* j = i; j > beg; --j ){
      const uchar* s = textBase + *(j-1);
      const uchar* t = textBase + *j;
      const indexT* m = maskPtr;

      while( *s == *t && *s < alphSize ){
	indexT off = *m;
	s += off;
	t += off;
	++m;
	if( m == maskEnd ) m = &mask[0];  // seems to be faster than modulus
      }

      if( *s <= *t ) break;
      std::iter_swap( j, j-1 );
    }
  }
}

void SuffixArray::insertionSortSimple( indexT* beg, indexT* end,
				       const uchar* textBase ){
  for( indexT* i = beg+1; i < end; ++i ){
    for( indexT* j = i; j > beg; --j ){
      const uchar* s = textBase + *(j-1);
      const uchar* t = textBase + *j;
      while( *s == *t && *s < alphSize ){
        ++s;
        ++t;
      }
      if( *s <= *t ) break;
      std::iter_swap( j, j-1 );
    }
  }
}

}  // end namespace cbrc
