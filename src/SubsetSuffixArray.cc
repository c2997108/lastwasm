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

void SubsetSuffixArray::setWordPositions(const DnaWordsFinder &finder,
					 size_t *cumulativeCounts,
					 const uchar *seqBeg,
					 const uchar *seqEnd,
					 size_t numOfThreads) {
  size_t numOfSeeds = seeds.size();
  size_t wordLength = finder.wordLength;
  size_t sumOfCounts = 0;
  for (size_t i = 0; i < numOfSeeds; ++i) {
    std::swap(cumulativeCounts[i], sumOfCounts);
  }
  resizePositions(sumOfCounts, seqEnd - seqBeg, numOfThreads);

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
				   int bitsPerInt,
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
  size_t wordLength = maxRestrictedSpan(&seeds[0], seeds.size());
  makeBucketStepsAndEnds(&bucketDepths[0], wordLength);
  size_t sufSize, bckSize, chiSize;

  if (bitsPerInt) {  // we have an old lastdb database :-(
    size_t chiBytes = bitsPerInt > 32 ? sizeof(size_t) : sizeof(unsigned);
    sufArray.bitsPerItem = 128 | bitsPerInt;
    bckArray.bitsPerItem = 128 | bitsPerInt;
    chiArray.bitsPerItem = 128 | chiBytes * CHAR_BIT;
    sufSize = bitsPerInt / CHAR_BIT * indexedPositions;
    bckSize = bitsPerInt / CHAR_BIT * bucketsSize();
    chiSize = chiBytes * indexedPositions;
  } else {  // we have a new lastdb database :-)
    sufArray.bitsPerItem = numOfBitsNeededFor(textLength - 1);
    bckArray.bitsPerItem = numOfBitsNeededFor(indexedPositions);
    chiArray.bitsPerItem = numOfBitsNeededFor(indexedPositions - 1);
    sufSize = numOfBytes(sufArray.bitsPerItem, indexedPositions);
    bckSize = numOfBytes(bckArray.bitsPerItem, bucketsSize());
    chiSize = numOfBytes(chiArray.bitsPerItem, indexedPositions);
  }

  suffixArray.m.open(baseName + ".suf", sufSize);
  sufArray.items = (const size_t *)suffixArray.m.begin();
  buckets.m.open(baseName + ".bck", bckSize);
  bckArray.items = (const size_t *)buckets.m.begin();

  // Prefer smaller child tables first (chi1 -> chi2 -> chi) to avoid
  // unnecessary failures in environments where one format may be absent.
  bool opened = false;
  try {
    chibiTable.m.open(baseName + ".chi1", indexedPositions);
    opened = true;
  } catch (std::runtime_error &) {}

  if (!opened) {
    try {
      kiddyTable.m.open(baseName + ".chi2", indexedPositions);
      opened = true;
    } catch (std::runtime_error &) {}
  }

  if (!opened) {
    try {
      childTable.m.open(baseName + ".chi", chiSize);
      chiArray.items = (const size_t *)childTable.m.begin();
      opened = true;
    } catch (std::runtime_error &) {}
  }
}

void SubsetSuffixArray::toFiles( const std::string& baseName,
				 bool isAppendPrj, size_t textLength ) const{
  size_t indexedPositions = bckArray[bucketEnds.back()];
  assert(textLength > indexedPositions);

  std::string fileName = baseName + ".prj";
  std::ofstream f( fileName.c_str(),
		   isAppendPrj ? std::ios::app : std::ios::out );

  f << "totallength=" << textLength << '\n';
  f << "specialcharacters=" << textLength - indexedPositions << '\n';

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

  size_t sufSize = numOfBytes(sufArray.bitsPerItem, indexedPositions);
  size_t bckSize = numOfBytes(bckArray.bitsPerItem, bucketsSize());
  memoryToBinaryFile(suffixArray.begin(), suffixArray.begin() + sufSize,
		     baseName + ".suf");
  memoryToBinaryFile(buckets.begin(), buckets.begin() + bckSize,
		     baseName + ".bck");

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
			const size_t *steps, int depth,
			ConstPackedArray sa, size_t saPos) {
  size_t position = sa[saPos];
  return bucketValue(seed, seed.firstMap(), steps, text + position, depth) + 1;
}

static void makeSomeBuckets(const uchar *text, const CyclicSubsetSeed &seed,
			    const size_t *steps, int depth,
			    PackedArray buckets, size_t buckBeg,
			    size_t buckIdx, size_t maxIdx, ConstPackedArray sa,
			    size_t saBeg, size_t saEnd) {
  for (size_t i = saBeg; buckIdx < maxIdx; ++i) {
    size_t b = buckBeg + bucketPos(text, seed, steps, depth, sa, i);
    b = std::min(b, maxIdx);  // for not-fully-sorted suffix array
    for (; buckIdx < b; ++buckIdx) {
      buckets.set(buckIdx, i);
    }
  }
}

static void runThreads(const uchar *text, const CyclicSubsetSeed *seedPtr,
		       const size_t *steps, int depth,
		       PackedArray buckets, size_t buckBeg,
		       size_t buckIdx, size_t maxIdx, ConstPackedArray sa,
		       size_t saBeg, size_t saEnd, size_t numOfThreads) {
  const CyclicSubsetSeed &seed = *seedPtr;
  if (numOfThreads > 1) {
    size_t len = (saEnd - saBeg) / numOfThreads;
    size_t mid = saEnd - len;
    size_t midIdx = buckBeg + bucketPos(text, seed, steps, depth, sa, mid - 1);
    midIdx = std::min(midIdx, maxIdx);  // for not-fully-sorted suffix array
    --numOfThreads;
    std::thread t(runThreads, text, seedPtr, steps, depth, buckets, buckBeg,
		  buckIdx, midIdx, sa, saBeg, mid, numOfThreads);
    size_t pad = numOfItemsBetweenWrites(buckets.bitsPerItem) * numOfThreads;
    makeSomeBuckets(text, seed, steps, depth, buckets, buckBeg+pad,
		    midIdx+pad, maxIdx+pad, sa, mid, saEnd);
    t.join();
    for (size_t i = midIdx; i < maxIdx; ++i) buckets.set(i, buckets[i+pad]);
    for (size_t i = maxIdx; i < maxIdx+pad; ++i) buckets.set(i, 0);
  } else {
    makeSomeBuckets(text, seed, steps, depth, buckets, buckBeg,
		    buckIdx, maxIdx, sa, saBeg, saEnd);
  }
}

void SubsetSuffixArray::makeBuckets(const uchar *text,
				    unsigned wordLength,
				    const size_t *cumulativeCounts,
				    size_t minPositionsPerBucket,
				    unsigned bucketDepth,
				    size_t numOfThreads) {
  bckArray.bitsPerItem = numOfBitsNeededFor(cumulativeCounts[seeds.size()-1]);
  std::vector<unsigned> bucketDepths(seeds.size(), bucketDepth);
  if (bucketDepth+1 == 0) {
    assert(minPositionsPerBucket > 0);
    unsigned long long pMem = sufArray.bitsPerItem;
    unsigned long long bMem = bckArray.bitsPerItem;
    size_t oldCount = 0;
    for (size_t s = 0; s < seeds.size(); ++s) {
      size_t newCount = cumulativeCounts[s];
      size_t m = (newCount - oldCount) * pMem / (minPositionsPerBucket * bMem);
      bucketDepths[s] = maxBucketDepth(seeds[s], 0, m, wordLength);
      oldCount = newCount;
    }
  }

  makeBucketStepsAndEnds(&bucketDepths[0], wordLength);
  size_t pad = totalItemsBetweenThreads(bckArray.bitsPerItem, numOfThreads);
  buckets.v.resize(numOfBytes(bckArray.bitsPerItem, bucketsSize() + pad));
  bckArray.items = (const size_t *)buckets.begin();

  PackedArray bucks = {(size_t *)&buckets.v[0], bckArray.bitsPerItem};
  size_t buckBeg = 0;
  size_t buckIdx = 0;
  size_t saBeg = 0;

  for (size_t s = 0; s < seeds.size(); ++s) {
    const CyclicSubsetSeed &seed = seeds[s];
    int depth = bucketDepths[s];
    const size_t *steps = bucketStepEnds[s];
    size_t saEnd = cumulativeCounts[s];
    if (saEnd > saBeg) {
      size_t maxIdx =
	buckBeg + bucketPos(text, seed, steps, depth, sufArray, saEnd - 1);
      runThreads(text, &seed, steps, depth, bucks, buckBeg,
		 buckIdx, maxIdx, sufArray, saBeg, saEnd, numOfThreads);
      buckIdx = maxIdx;
      saBeg = saEnd;
    }
    buckBeg += steps[0];
  }

  for (; buckIdx <= buckBeg; ++buckIdx) {
    bucks.set(buckIdx, saBeg);
  }
}

void SubsetSuffixArray::makeBucketStepsAndEnds(const unsigned *bucketDepths,
					       size_t wordLength) {
  size_t numOfSeeds = seeds.size();
  size_t numOfBucketSteps = numOfSeeds;
  for (size_t i = 0; i < numOfSeeds; ++i) {
    numOfBucketSteps += bucketDepths[i];
  }
  bucketSteps.resize(numOfBucketSteps);
  bucketStepEnds.resize(numOfSeeds + 1);
  bucketEnds.resize(numOfSeeds + 1);

  size_t *steps = &bucketSteps[0];
  size_t end = 0;

  for (size_t i = 0; i < numOfSeeds; ++i) {
    bucketStepEnds[i] = steps;
    bucketEnds[i] = end;
    unsigned depth = bucketDepths[i];
    makeBucketSteps(steps, seeds[i], 0, depth, wordLength);
    end += steps[0];
    steps += depth + 1;
  }

  bucketStepEnds[numOfSeeds] = steps;
  bucketEnds[numOfSeeds] = end;
}
