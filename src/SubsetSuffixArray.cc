// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

#include "SubsetSuffixArray.hh"
#include "SubsetMinimizerFinder.hh"
#include "io.hh"
#include <cassert>
#include <cstdio>  // remove
#include <sstream>

using namespace cbrc;

static void err(const std::string &s) {
  throw std::runtime_error(s);
}

static unsigned maxBucketDepth(const CyclicSubsetSeed &seed,
			       size_t maxBucketItems) {
  unsigned depth = 0;
  size_t kmerItems = 1;
  size_t bucketItems = 1;
  while (1) {
    size_t c = seed.subsetCount(depth);
    if (kmerItems > maxBucketItems / c) return depth;
    kmerItems *= c;
    if (bucketItems > maxBucketItems - kmerItems) return depth;
    bucketItems += kmerItems;
    ++depth;
  }
}

void SubsetSuffixArray::addPositions(const uchar* text, indexT beg, indexT end,
				     size_t step, size_t minimizerWindow) {
  if (beg >= end) return;
  assert(step > 0);
  const CyclicSubsetSeed &seed = seeds[0];
  const uchar *subsetMap = seed.firstMap();
  SubsetMinimizerFinder f;
  f.init(seed, text + beg, text + end);

  while (true) {
    if (minimizerWindow > 1) {
      if (f.isMinimizer(seed, text + beg, text + end, minimizerWindow))
	suffixArray.v.push_back(beg);
    } else {
      if (subsetMap[text[beg]] < CyclicSubsetSeed::DELIMITER)
	suffixArray.v.push_back(beg);
    }
    if (end - beg <= step) break;  // avoid overflow
    beg += step;
  }
}

void SubsetSuffixArray::fromFiles( const std::string& baseName,
				   bool isMaskLowercase,
				   const uchar letterCode[] ){
  size_t textLength = 0;  // 0 never occurs in a valid file
  size_t unindexedPositions = 0;  // 0 never occurs in a valid file
  unsigned bucketDepth = -1;
  unsigned version = 0;
  seeds.clear();
  seeds.resize(1);

  std::string fileName = baseName + ".prj";
  std::ifstream f( fileName.c_str() );
  if (!f) err("can't open file: " + fileName);

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "version" ) iss >> version;
    if( word == "totallength" ) iss >> textLength;
    if( word == "specialcharacters" ) iss >> unindexedPositions;
    if( word == "prefixlength" ) iss >> bucketDepth;
    if( word == "subsetseed" ){
      seeds.back().appendPosition(iss, isMaskLowercase, letterCode);
    }
  }

  if( textLength == 0 || unindexedPositions == 0 || bucketDepth+1 == 0 ||
      !seeds.back().span() || !f.eof() ){
    err("can't read file: " + fileName);
  }

  if (bucketDepth == 0 && version < 1087) {
    err("the lastdb files are too old: please re-run lastdb");
  }

  size_t indexedPositions = textLength - unindexedPositions;
  suffixArray.m.open( baseName + ".suf", indexedPositions );
  makeBucketSteps( bucketDepth );
  buckets.m.open(baseName + ".bck", bucketsSize());

  try{
    childTable.m.open( baseName + ".chi", indexedPositions );
  }catch( std::runtime_error& ){
    try{
      kiddyTable.m.open( baseName + ".chi2", indexedPositions );
    }catch( std::runtime_error& ){
      try{
	chibiTable.m.open( baseName + ".chi1", indexedPositions );
      }catch( std::runtime_error& ){}
    }
  }
}

void SubsetSuffixArray::toFiles( const std::string& baseName,
				 bool isAppendPrj, size_t textLength ) const{
  assert( textLength > suffixArray.size() );

  std::string fileName = baseName + ".prj";
  std::ofstream f( fileName.c_str(),
		   isAppendPrj ? std::ios::app : std::ios::out );

  f << "totallength=" << textLength << '\n';
  f << "specialcharacters=" << textLength - suffixArray.size() << '\n';
  f << "prefixlength=" << maxBucketPrefix() << '\n';

  for (size_t s = 0; s < seeds.size(); ++s) {
    for (size_t i = 0; i < seeds[s].span(); ++i) {
      f << "subsetseed=";
      seeds[s].writePosition(f, i);
      f << '\n';
    }
  }

  f.close();
  if (!f) err("can't write file: " + fileName);

  memoryToBinaryFile( suffixArray.begin(), suffixArray.end(),
		      baseName + ".suf" );
  memoryToBinaryFile( buckets.begin(), buckets.end(), baseName + ".bck" );

  fileName = baseName + ".chi";
  std::remove( fileName.c_str() );
  memoryToBinaryFile( childTable.begin(), childTable.end(), fileName );

  fileName = baseName + ".chi2";
  std::remove( fileName.c_str() );
  memoryToBinaryFile( kiddyTable.begin(), kiddyTable.end(), fileName );

  fileName = baseName + ".chi1";
  std::remove( fileName.c_str() );
  memoryToBinaryFile( chibiTable.begin(), chibiTable.end(), fileName );
}

void SubsetSuffixArray::makeBuckets(const uchar *text,
				    const size_t *cumulativeCounts,
				    unsigned bucketDepth) {
  std::vector<unsigned> bucketDepths(seeds.size(), bucketDepth);
  if (bucketDepth+1 == 0) {
    size_t minPositionsPerBucket = 4;
    size_t oldCount = 0;
    for (size_t s = 0; s < seeds.size(); ++s) {
      size_t newCount = cumulativeCounts[s];
      size_t maxBucketItems = (newCount - oldCount) / minPositionsPerBucket;
      bucketDepths[s] = maxBucketDepth(seeds[s], maxBucketItems);
      oldCount = newCount;
    }
  }

  makeBucketSteps(bucketDepths[0]);
  buckets.v.resize(bucketsSize());

  indexT *myBuckets = &buckets.v[0];
  indexT *bucketPtr = myBuckets;
  const CyclicSubsetSeed &seed = seeds[0];
  unsigned myBucketDepth = bucketDepths[0];

  for( indexT i = 0; i < suffixArray.size(); ++i ){
    const uchar* textPtr = text + suffixArray[i];
    const uchar* subsetMap = seed.firstMap();
    indexT bucketIndex = 0;
    unsigned depth = 0;

    while( depth < myBucketDepth ){
      uchar subset = subsetMap[ *textPtr ];
      if( subset == CyclicSubsetSeed::DELIMITER ){
	bucketIndex += bucketSteps[depth] - 1;  // depth > 0
	break;
      }
      ++textPtr;
      ++depth;
      bucketIndex += subset * bucketSteps[depth];
      subsetMap = seed.nextMap( subsetMap );
    }

    indexT *newBucketPtr = myBuckets + bucketIndex + 1;
    if (newBucketPtr > bucketPtr) {
      std::fill(bucketPtr, newBucketPtr, i);
    }
    bucketPtr = newBucketPtr;
  }

  std::fill(bucketPtr, myBuckets + bucketSteps[0] + 1, suffixArray.size());
}

void SubsetSuffixArray::makeBucketSteps(unsigned bucketDepth) {
  indexT step = 0;
  indexT depth = bucketDepth + 1;
  bucketSteps.resize( depth );
  const CyclicSubsetSeed &seed = seeds[0];

  while( depth > 0 ){
    --depth;
    step = step * seed.subsetCount(depth);
    if (depth != 0 || bucketDepth == 0) ++step;
    bucketSteps[depth] = step;
  }
}
