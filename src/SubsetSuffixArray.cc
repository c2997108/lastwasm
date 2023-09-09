// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

#include "SubsetSuffixArray.hh"
#include "io.hh"
#include <cassert>
#include <cstdio>  // remove
#include <sstream>

#include <thread>

using namespace cbrc;

static void err(const std::string &s) {
  throw std::runtime_error(s);
}

static void offSet(OffPart *p, size_t value) {
  for (int i = 0; i < offParts; ++i) {
    p[i] = value >> (i * sizeof(OffPart) * CHAR_BIT);
  }
}

static unsigned maxBucketDepth(const CyclicSubsetSeed &seed,
			       size_t maxBucketItems, unsigned wordLength) {
  unsigned depth = 0;
  size_t kmerItems = 1;
  while (depth < wordLength) {
    size_t c = seed.restrictedSubsetCount(depth);
    if (kmerItems > maxBucketItems / c) return depth;
    kmerItems *= c;
    ++depth;
  }
  size_t bucketItems = kmerItems;
  while (1) {
    size_t c = seed.unrestrictedSubsetCount(depth);
    if (kmerItems > maxBucketItems / c) return depth;
    kmerItems *= c;
    if (bucketItems > maxBucketItems - kmerItems) return depth;
    bucketItems += kmerItems;
    ++depth;
  }
}

void SubsetSuffixArray::setWordPositions(const DnaWordsFinder &finder,
					 size_t *cumulativeCounts,
					 const uchar *seqBeg,
					 const uchar *seqEnd) {
  size_t numOfSeeds = seeds.size();
  size_t wordLength = finder.wordLength;
  size_t sumOfCounts = 0;
  for (size_t i = 0; i < numOfSeeds; ++i) {
    std::swap(cumulativeCounts[i], sumOfCounts);
  }
  resizePositions(sumOfCounts);

  unsigned hash = 0;
  const uchar *seqPos = finder.init(seqBeg, seqEnd, &hash);
  while (seqPos < seqEnd) {
    unsigned c = finder.baseToCode[*seqPos];
    ++seqPos;
    if (c != dnaWordsFinderNull) {
      unsigned w = finder.next(&hash, c);
      if (w != dnaWordsFinderNull) {
	setPosition(cumulativeCounts[w]++, seqPos - seqBeg - wordLength);
      }
    } else {
      seqPos = finder.init(seqPos, seqEnd, &hash);
    }
  }
}

void SubsetSuffixArray::fromFiles( const std::string& baseName,
				   bool isMaskLowercase,
				   const uchar letterCode[],
				   const std::string &mainSequenceAlphabet ){
  size_t textLength = 0;  // 0 never occurs in a valid file
  size_t unindexedPositions = 0;  // 0 never occurs in a valid file
  unsigned version = 0;
  std::vector<unsigned> bucketDepths;
  seeds.clear();

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
    if( word == "prefixlength" ){
      if (!seeds.empty() && !seeds.back().span()) {
	err("can't read file: " + fileName);
      }
      unsigned d = 0;
      iss >> d;
      if (!iss) err("can't read file: " + fileName);
      bucketDepths.push_back(d);
      seeds.resize(seeds.size() + 1);
    }
    if( word == "subsetseed" ){
      if (seeds.empty()) err("can't read file: " + fileName);
      seeds.back().appendPosition(iss, isMaskLowercase, letterCode,
				  mainSequenceAlphabet);
    }
  }

  if (textLength == 0 || unindexedPositions == 0 || seeds.empty() ||
      !seeds.back().span() || !f.eof()) {
    err("can't read file: " + fileName);
  }

  if (bucketDepths[0] == 0 && version < 1087) {
    err("the lastdb files are too old: please re-run lastdb");
  }

  size_t indexedPositions = textLength - unindexedPositions;
  suffixArray.m.open(baseName + ".suf", indexedPositions * posParts);

  size_t wordLength = maxRestrictedSpan(&seeds[0], seeds.size());
  makeBucketSteps(&bucketDepths[0], wordLength);
  buckets.m.open(baseName + ".bck", offParts * bucketsSize());
  initBucketEnds();

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
  assert(textLength > size());

  std::string fileName = baseName + ".prj";
  std::ofstream f( fileName.c_str(),
		   isAppendPrj ? std::ios::app : std::ios::out );

  f << "totallength=" << textLength << '\n';
  f << "specialcharacters=" << textLength - size() << '\n';

  for (size_t s = 0; s < seeds.size(); ++s) {
    f << "prefixlength=" << maxBucketPrefix(s) << '\n';
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

static size_t bucketPos(const uchar *text, const CyclicSubsetSeed &seed,
			const size_t *steps, unsigned depth,
			const PosPart *sa, size_t saPos) {
  const uchar *textPtr = text + posGet(sa, saPos);

  size_t bucketIndex = 0;
  const uchar *subsetMap = seed.firstMap();
  unsigned d = 0;
  while (d < depth) {
    uchar subset = subsetMap[*textPtr];
    if (subset == CyclicSubsetSeed::DELIMITER) {
      return bucketIndex + steps[d];  // d > 0
    }
    ++textPtr;
    ++d;
    bucketIndex += subset * steps[d];
    subsetMap = seed.nextMap(subsetMap);
  }

  return bucketIndex + 1;
}

static void makeSomeBuckets(const uchar *text, const CyclicSubsetSeed &seed,
			    const size_t *steps, unsigned depth,
			    OffPart *buckets, size_t buckBeg, size_t buckIdx,
			    const PosPart *sa, size_t saBeg, size_t saEnd) {
  for (size_t i = saBeg; i < saEnd; ++i) {
    size_t b = buckBeg + bucketPos(text, seed, steps, depth, sa, i);
    for (; buckIdx < b; ++buckIdx) {
      offSet(buckets + offParts * buckIdx, i);
    }
  }
}

static void runThreads(const uchar *text, const CyclicSubsetSeed *seedPtr,
		       const size_t *steps, unsigned depth, OffPart *buckets,
		       size_t buckBeg, size_t buckIdx, const PosPart *sa,
		       size_t saBeg, size_t saEnd, size_t numOfThreads) {
  const CyclicSubsetSeed &seed = *seedPtr;
  if (numOfThreads > 1) {
    size_t len = (saEnd - saBeg + numOfThreads - 1) / numOfThreads;
    size_t mid = saBeg + len;
    size_t b = buckBeg + bucketPos(text, seed, steps, depth, sa, mid - 1);
    std::thread t(runThreads, text, seedPtr, steps, depth, buckets, buckBeg,
		  b, sa, mid, saEnd, numOfThreads - 1);
    makeSomeBuckets(text, seed, steps, depth, buckets, buckBeg,
		    buckIdx, sa, saBeg, mid);
    t.join();
  } else {
    makeSomeBuckets(text, seed, steps, depth, buckets, buckBeg,
		    buckIdx, sa, saBeg, saEnd);
  }
}

void SubsetSuffixArray::makeBuckets(const uchar *text,
				    unsigned wordLength,
				    const size_t *cumulativeCounts,
				    size_t minPositionsPerBucket,
				    unsigned bucketDepth,
				    size_t numOfThreads) {
  std::vector<unsigned> bucketDepths(seeds.size(), bucketDepth);
  if (bucketDepth+1 == 0) {
    assert(minPositionsPerBucket > 0);
    size_t oldCount = 0;
    for (size_t s = 0; s < seeds.size(); ++s) {
      size_t newCount = cumulativeCounts[s];
      size_t maxBucketItems = (newCount - oldCount) / minPositionsPerBucket;
      bucketDepths[s] = maxBucketDepth(seeds[s], maxBucketItems, wordLength);
      oldCount = newCount;
    }
  }

  makeBucketSteps(&bucketDepths[0], wordLength);
  buckets.v.resize(offParts * bucketsSize());
  initBucketEnds();

  const PosPart *sa = suffixArray.begin();
  OffPart *bucks = &buckets.v[0];
  size_t buckBeg = 0;
  size_t buckIdx = 0;
  size_t saBeg = 0;

  for (size_t s = 0; s < seeds.size(); ++s) {
    const CyclicSubsetSeed &seed = seeds[s];
    unsigned depth = bucketDepths[s];
    const size_t *steps = bucketStepEnds[s];
    size_t saEnd = cumulativeCounts[s];
    if (saEnd > saBeg) {
      runThreads(text, &seed, steps, depth, bucks, buckBeg,
		 buckIdx, sa, saBeg, saEnd, numOfThreads);
      buckIdx = buckBeg + bucketPos(text, seed, steps, depth, sa, saEnd - 1);
      saBeg = saEnd;
    }
    buckBeg += steps[0];
  }

  for (; buckIdx <= buckBeg; ++buckIdx) {
    offSet(bucks + offParts * buckIdx, saBeg);
  }
}

static void makeBucketStepsForOneSeed(size_t *steps, unsigned depth,
				      const CyclicSubsetSeed &seed,
				      size_t wordLength) {
  size_t step = 1;
  steps[depth] = step;

  while (depth > 0) {
    --depth;
    if (depth < wordLength) {
      step = step * seed.restrictedSubsetCount(depth);
    } else {
      step = step * seed.unrestrictedSubsetCount(depth) + (depth > 0);
      // Add one for delimiters, except when depth==0
    }
    steps[depth] = step;
  }
}

void SubsetSuffixArray::makeBucketSteps(const unsigned *bucketDepths,
					size_t wordLength) {
  size_t numOfSeeds = seeds.size();
  size_t numOfBucketSteps = numOfSeeds;
  for (size_t i = 0; i < numOfSeeds; ++i) {
    numOfBucketSteps += bucketDepths[i];
  }
  bucketSteps.resize(numOfBucketSteps);
  bucketStepEnds.resize(numOfSeeds + 1);
  size_t *steps = &bucketSteps[0];
  for (size_t i = 0; i < numOfSeeds; ++i) {
    bucketStepEnds[i] = steps;
    unsigned depth = bucketDepths[i];
    makeBucketStepsForOneSeed(steps, depth, seeds[i], wordLength);
    steps += depth + 1;
  }
  bucketStepEnds[numOfSeeds] = steps;
}
