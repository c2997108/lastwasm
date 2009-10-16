// Copyright 2008, 2009 Martin C. Frith

#include "SuffixArray.hh"
#include "PeriodicSpacedSeed.hh"
#include "io.hh"
#include <cassert>
//#include <iostream>  // for debugging

using namespace cbrc;

int SuffixArray::addIndices( const uchar* text,
			     indexT beg, indexT end, indexT step,
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

void SuffixArray::clear(){
  index.clear();
  buckets.clear();
  bucketSteps.clear();
  bucketMask.clear();
}

void SuffixArray::fromFiles( const std::string& baseName,
			     indexT indexNum, indexT bucketDepth,
			     const PeriodicSpacedSeed& seed ){
  index.resize(indexNum);  // unwanted zero-fill
  vectorFromBinaryFile( index, baseName + ".suf" );

  makeBucketSteps( bucketDepth );
  makeBucketMask( seed, bucketDepth );
  if( bucketDepth == 0 ) return;

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
			 const uchar* queryPtr, const uchar* text,
			 const PeriodicSpacedSeed& seed,
			 indexT maxHits, indexT minDepth ) const{
  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

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
  const uchar* textBase = text + bucketMaskTotal;

  for( indexT depth = bucketDepth; /* noop */; ++depth ){
    uchar symbol = *queryPtr;
    if( symbol < alphSize ){
      equalRange( beg, end, textBase, symbol );
    }else{
      beg = end;
      depth = minDepth;
    }
    if( indexT(end - beg) <= maxHits && depth+1 >= minDepth ) return;
    indexT off = seed.offsets[ depth % seed.weight ];
    queryPtr += off;
    textBase += off;
  }
}

void SuffixArray::countMatches( std::vector<unsigned long long>& counts,
				const uchar* queryPtr, const uchar* text,
				const PeriodicSpacedSeed& seed ) const{
  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

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
  const uchar* textBase = text + bucketMaskTotal;

  for( indexT depth = bucketDepth; /* noop */; ++depth ){
    uchar symbol = *queryPtr;
    if( symbol >= alphSize ) return;
    equalRange( beg, end, textBase, symbol );
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    indexT off = seed.offsets[ depth % seed.weight ];
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

void SuffixArray::makeBuckets( const uchar* text,
			       const PeriodicSpacedSeed& seed,
			       indexT bucketDepth ){
  if( bucketDepth == -1u ) bucketDepth = defaultBucketDepth();

  makeBucketSteps( bucketDepth );
  makeBucketMask( seed, bucketDepth );
  if( bucketDepth == 0 ) return;

  for( indexT i = 0; i < index.size(); ++i ){
    indexT bucketIndex = 0;
    const uchar* textPtr = text + index[i];

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

SuffixArray::indexT SuffixArray::defaultBucketDepth(){
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
