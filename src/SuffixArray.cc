// Copyright 2008, 2009 Martin C. Frith

#include "SuffixArray.hh"
#include "PeriodicSpacedSeed.hh"
#include "io.hh"
#include <cassert>
//#include <iostream>  // for debugging

using namespace cbrc;

int SuffixArray::addIndices( const uchar* text,
			     indexT beg, indexT end, indexT step,
			     unsigned alphSize, std::size_t maxBytes ){
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

void SuffixArray::clear(){
  index.clear();
  buckets.clear();
  bucketSteps.clear();
  bucketMask.clear();
}

void SuffixArray::fromFiles( const std::string& baseName,
			     indexT indexNum, indexT bucketDepth,
			     const PeriodicSpacedSeed& seed,
			     unsigned alphSize ){
  index.resize(indexNum);  // unwanted zero-fill
  vectorFromBinaryFile( index, baseName + ".suf" );

  makeBucketSteps( alphSize, bucketDepth );
  makeBucketMask( seed, bucketDepth );

  buckets.resize( bucketSteps[0] );  // unwanted zero-fill
  vectorFromBinaryFile( buckets, baseName + ".bck" );
}

void SuffixArray::toFiles( const std::string& baseName ) const{
  vectorToBinaryFile( index, baseName + ".suf" );
  vectorToBinaryFile( buckets, baseName + ".bck" );
}

// use past results to speed up long matches?
// could & probably should return the match depth
void SuffixArray::match( const indexT*& beg, const indexT*& end,
			 const uchar* queryPtr, const uchar* text,
			 const PeriodicSpacedSeed& seed, unsigned alphSize,
			 indexT maxHits, indexT minDepth ) const{
  indexT depth = 0;

  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

  while( depth < bucketDepth ){
    if( bucketEnd - bucketBeg <= maxHits && depth >= minDepth ){
      beg = &index[0] + bucketBeg;
      end = &index[0] + bucketEnd;
      return;
    }
    uchar symbol = *queryPtr;
    if( symbol < alphSize ){
      queryPtr += bucketMask[depth];
      ++depth;
      indexT step = bucketSteps[depth];
      bucketPtr += symbol * step;
      bucketBeg = *bucketPtr;
      bucketEnd = *(bucketPtr + step);
    }else{  // we hit a delimiter in the query, so finish without any matches:
      bucketBeg = bucketEnd;
      minDepth = 0;
    }
  }

  // match using binary search:
  beg = &index[0] + bucketBeg;
  end = &index[0] + bucketEnd;
  const uchar* textBase = text + bucketMaskTotal;

  while( true ){
    if( indexT(end - beg) <= maxHits && depth >= minDepth ) return;
    uchar symbol = *queryPtr;
    if( symbol < alphSize ){
      equalRange( beg, end, textBase, symbol );
      indexT off = seed.offsets[ depth % seed.weight ];
      queryPtr += off;
      textBase += off;
      ++depth;
    }else{  // we hit a delimiter in the query, so finish without any matches:
      beg = end;
      minDepth = 0;
    }
  }
}

void SuffixArray::countMatches( std::vector<unsigned long long>& counts,
				const uchar* queryPtr, const uchar* text,
				const PeriodicSpacedSeed& seed,
				unsigned alphSize ) const{
  indexT depth = 0;

  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

  while( depth < bucketDepth ){
    if( bucketBeg == bucketEnd ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += bucketEnd - bucketBeg;
    uchar symbol = *queryPtr;
    if( symbol >= alphSize ) return;
    queryPtr += bucketMask[depth];
    ++depth;
    indexT step = bucketSteps[depth];
    bucketPtr += symbol * step;
    bucketBeg = *bucketPtr;
    bucketEnd = *(bucketPtr + step);
  }

  // match using binary search:
  const indexT* beg = &index[0] + bucketBeg;
  const indexT* end = &index[0] + bucketEnd;
  const uchar* textBase = text + bucketMaskTotal;

  while( true ){
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    uchar symbol = *queryPtr;
    if( symbol >= alphSize ) return;
    equalRange( beg, end, textBase, symbol );
    indexT off = seed.offsets[ depth % seed.weight ];
    queryPtr += off;
    textBase += off;
    ++depth;
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

void SuffixArray::makeBuckets( const uchar* text,
			       const PeriodicSpacedSeed& seed,
			       unsigned alphSize, indexT bucketDepth ){
  if( bucketDepth == -1u ) bucketDepth = defaultBucketDepth( alphSize );

  makeBucketSteps( alphSize, bucketDepth );
  makeBucketMask( seed, bucketDepth );

  for( indexT i = 0; i < index.size(); ++i ){
    const uchar* textPtr = text + index[i];
    indexT bucketIndex = 0;
    indexT depth = 0;

    while( depth < bucketDepth ){
      unsigned symbol = *textPtr;
      if( symbol >= alphSize ){
	bucketIndex += bucketSteps[depth] - 1;
	break;
      }
      textPtr += bucketMask[depth];
      ++depth;
      indexT step = bucketSteps[depth];
      bucketIndex += symbol * step;
    }

    buckets.resize( bucketIndex+1, i );
  }

  buckets.resize( bucketSteps[0], index.size() );
}

void SuffixArray::makeBucketSteps( unsigned alphSize, indexT bucketDepth ){
  indexT step = 0;
  indexT depth = bucketDepth + 1;
  bucketSteps.resize( depth );

  while( depth > 0 ){
    --depth;
    step = step * alphSize + 1;
    bucketSteps[depth] = step;
  }
}

void SuffixArray::makeBucketMask( const PeriodicSpacedSeed& seed,
				  indexT bucketDepth ){
  bucketMaskTotal = 0;  // necessary!!!
  bucketMask.resize(bucketDepth);

  for( indexT i = 0; i < bucketDepth; ++i ){
    indexT offset = seed.offsets[ i % seed.weight ];
    bucketMask[i] = offset;
    bucketMaskTotal += offset;
  }
}

SuffixArray::indexT SuffixArray::defaultBucketDepth( unsigned alphSize ){
  indexT maxBucketEntries = index.size() / 4;
  indexT bucketDepth = 0;
  indexT kmerEntries = 1;
  indexT bucketEntries = 1;

  while(true){
    if( kmerEntries > maxBucketEntries / alphSize ) return bucketDepth;
    kmerEntries *= alphSize;
    if( bucketEntries > maxBucketEntries - kmerEntries ) return bucketDepth;
    bucketEntries += kmerEntries;
    ++bucketDepth;
  }
}
