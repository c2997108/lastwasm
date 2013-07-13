// Copyright 2008, 2009, 2010, 2013 Martin C. Frith

#include "SubsetSuffixArray.hh"
#include "io.hh"
#include <cassert>
//#include <iostream>  // for debugging

using namespace cbrc;

int SubsetSuffixArray::addIndices( const uchar* text,
				   indexT beg, indexT end, indexT step,
				   std::size_t maxBytes ){
  assert( step > 0 );
  indexT oldSize = index.v.size();
  const uchar* subsetMap = seed.firstMap();

  for( indexT i = beg; i < end; i += step ){
    if( subsetMap[ text[i] ] < CyclicSubsetSeed::DELIMITER ){
      index.v.push_back(i);

      if( indexBytes() > maxBytes ){
	index.v.erase( index.v.begin() + oldSize, index.v.end() );
	return 0;
      }
    }
  }

  return 1;
}

void SubsetSuffixArray::clear(){
  index.v.clear();
  buckets.v.clear();
  bucketSteps.clear();
}

void SubsetSuffixArray::fromFiles( const std::string& baseName,
				   indexT indexNum, indexT bucketDepth ){
  index.m.open( baseName + ".suf", indexNum );
  makeBucketSteps( bucketDepth );
  buckets.m.open( baseName + ".bck", bucketSteps[0] );
}

void SubsetSuffixArray::toFiles( const std::string& baseName,
				 indexT textLength ) const{
  std::string fileName = baseName + ".prj";
  std::ofstream f( fileName.c_str(), std::ios::app );

  f << "totallength=" << textLength << '\n';
  f << "specialcharacters=" << textLength - index.size() << '\n';
  f << "prefixlength=" << maxBucketPrefix() << '\n';

  if( !f ) throw std::runtime_error( "can't write file: " + fileName );

  memoryToBinaryFile( index.begin(), index.end(), baseName + ".suf" );
  memoryToBinaryFile( buckets.begin(), buckets.end(), baseName + ".bck" );
}

// use past results to speed up long matches?
// could & probably should return the match depth
void SubsetSuffixArray::match( const indexT*& beg, const indexT*& end,
                               const uchar* queryPtr, const uchar* text,
                               indexT maxHits, indexT minDepth ) const{
  indexT depth = 0;
  const uchar* subsetMap = seed.firstMap();

  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

  while( depth < bucketDepth ){
    if( bucketEnd - bucketBeg <= maxHits || depth >= minDepth ){
      beg = &index[0] + bucketBeg;
      end = &index[0] + bucketEnd;
      return;
    }
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset < CyclicSubsetSeed::DELIMITER ){
      ++depth;
      indexT step = bucketSteps[depth];
      bucketPtr += subset * step;
      bucketBeg = *bucketPtr;
      bucketEnd = *(bucketPtr + step);
      subsetMap = seed.nextMap( subsetMap );
    }else{  // we hit a delimiter in the query, so finish without any matches:
      bucketBeg = bucketEnd;
    }
  }

  // match using binary search:
  beg = &index[0] + bucketBeg;
  end = &index[0] + bucketEnd;

  while( true ){
    if( indexT(end - beg) <= maxHits || depth >= minDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset < CyclicSubsetSeed::DELIMITER ){
      equalRange( beg, end, text+depth, subsetMap, subset );
      ++depth;
      subsetMap = seed.nextMap( subsetMap );
    }else{  // we hit a delimiter in the query, so finish without any matches:
      beg = end;
    }
  }
}

void SubsetSuffixArray::countMatches( std::vector<unsigned long long>& counts,
				      const uchar* queryPtr,
				      const uchar* text ) const{
  indexT depth = 0;
  const uchar* subsetMap = seed.firstMap();

  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

  while( depth < bucketDepth ){
    if( bucketBeg == bucketEnd ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += bucketEnd - bucketBeg;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    ++depth;
    indexT step = bucketSteps[depth];
    bucketPtr += subset * step;
    bucketBeg = *bucketPtr;
    bucketEnd = *(bucketPtr + step);
    subsetMap = seed.nextMap( subsetMap );
  }

  // match using binary search:
  const indexT* beg = &index[0] + bucketBeg;
  const indexT* end = &index[0] + bucketEnd;

  while( true ){
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    equalRange( beg, end, text+depth, subsetMap, subset );
    ++depth;
    subsetMap = seed.nextMap( subsetMap );
  }
}

void SubsetSuffixArray::equalRange( const indexT*& beg, const indexT*& end,
				    const uchar* textBase,
				    const uchar* subsetMap, uchar subset ){
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    uchar s = subsetMap[ textBase[ *mid ] ];
    if( s < subset ){
      beg = mid + 1;
    }else if( s > subset ){
      end = mid;
    }else{
      if( subset > 0 )  // this "if" is unnecessary, but makes it a bit faster
	beg = lowerBound( beg, mid, textBase, subsetMap, subset );
      end = upperBound( mid + 1, end, textBase, subsetMap, subset );
      return;
    }
  }
}

const SubsetSuffixArray::indexT*
SubsetSuffixArray::lowerBound( const indexT* beg, const indexT* end,
			       const uchar* textBase,
			       const uchar* subsetMap, uchar subset ){
  for( ;; ){
    std::size_t size = end - beg;
    if( size <= 4 ) break;  // 3,4 seem good for hg18 chr21 versus itself
    const indexT* mid = beg + size / 2;
    if( subsetMap[ textBase[ *mid ] ] < subset ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }

  while( subsetMap[ textBase[ *beg ] ] < subset ) ++beg;  // linear search

  return beg;
}

const SubsetSuffixArray::indexT*
SubsetSuffixArray::upperBound( const indexT* beg, const indexT* end,
			       const uchar* textBase,
			       const uchar* subsetMap, uchar subset ){
  for( ;; ){
    std::size_t size = end - beg;
    if( size <= 4 ) break;  // 3,4 seem good for hg18 chr21 versus itself
    const indexT* mid = beg + size / 2;
    if( subsetMap[ textBase[ *mid ] ] <= subset ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }

  while( subsetMap[ textBase[ *(end-1) ] ] > subset ) --end;  // linear search

  return end;
}

void SubsetSuffixArray::makeBuckets( const uchar* text, indexT bucketDepth ){
  if( bucketDepth+1 == 0 ) bucketDepth = defaultBucketDepth();

  makeBucketSteps( bucketDepth );

  for( indexT i = 0; i < index.size(); ++i ){
    const uchar* textPtr = text + index[i];
    const uchar* subsetMap = seed.firstMap();
    indexT bucketIndex = 0;
    indexT depth = 0;

    while( depth < bucketDepth ){
      uchar subset = subsetMap[ *textPtr ];
      if( subset == CyclicSubsetSeed::DELIMITER ){
	bucketIndex += bucketSteps[depth] - 1;
	break;
      }
      ++textPtr;
      ++depth;
      indexT step = bucketSteps[depth];
      bucketIndex += subset * step;
      subsetMap = seed.nextMap( subsetMap );
    }

    buckets.v.resize( bucketIndex+1, i );
  }

  buckets.v.resize( bucketSteps[0], index.size() );
}

void SubsetSuffixArray::makeBucketSteps( indexT bucketDepth ){
  indexT step = 0;
  indexT depth = bucketDepth + 1;
  bucketSteps.resize( depth );

  while( depth > 0 ){
    --depth;
    step = step * seed.subsetCount(depth) + 1;
    bucketSteps[depth] = step;
  }
}

SubsetSuffixArray::indexT SubsetSuffixArray::defaultBucketDepth(){
  indexT maxBucketEntries = index.size() / 4;
  indexT bucketDepth = 0;
  indexT kmerEntries = 1;
  indexT bucketEntries = 1;

  while(true){
    indexT nextSubsetCount = seed.subsetCount(bucketDepth);
    if( kmerEntries > maxBucketEntries / nextSubsetCount ) return bucketDepth;
    kmerEntries *= nextSubsetCount;
    if( bucketEntries > maxBucketEntries - kmerEntries ) return bucketDepth;
    bucketEntries += kmerEntries;
    ++bucketDepth;
  }
}
