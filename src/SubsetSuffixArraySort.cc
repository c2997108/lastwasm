// Copyright 2008, 2009, 2010, 2011, 2013 Martin C. Frith

// Parts of this code are adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#include "SubsetSuffixArray.hh"

#include <string.h>  // memcpy

#include <algorithm>  // swap, min
//#include <iostream>  // for debugging

#ifdef HAS_CXX_THREADS
#include <thread>
#endif

using namespace cbrc;

namespace{
  typedef SubsetSuffixArray::Range Range;
}

static void pushRange(std::vector<Range> &v,
		      size_t beg, size_t end, size_t depth) {
  if (end - beg > 1) {
    Range r = {beg, end, depth};
    v.push_back(r);
  }
}

static void insertionSort(const uchar *text, const CyclicSubsetSeed &seed,
			  size_t *positions, size_t beg, size_t end,
			  const uchar *subsetMap) {
  for (size_t i = beg + 1; i < end; ++i) {
    size_t newPos = positions[i];
    size_t j = i;
    do {
      size_t oldPos = positions[j - 1];
      if (!seed.isLess(text + newPos, text + oldPos, subsetMap)) break;
      positions[j] = oldPos;
      --j;
    } while (j > beg);
    positions[j] = newPos;
  }
}

void SubsetSuffixArray::sort2(const uchar *text, const CyclicSubsetSeed &seed,
			      size_t *positions, size_t origin,
			      size_t beg, const uchar *subsetMap) {
  size_t mid = beg + 1;

  size_t b = positions[beg];
  size_t m = positions[mid];
  const uchar *s = text + b;
  const uchar *t = text + m;
  while (true) {
    uchar x = subsetMap[*s];
    uchar y = subsetMap[*t];
    if (x != y) {
      if (x > y) { positions[beg] = m; positions[mid] = b; }
      break;
    }
    if (x == CyclicSubsetSeed::DELIMITER) return;
    ++s;
    ++t;
    subsetMap = seed.nextMap(subsetMap);
  }

  if (isChildDirectionForward(origin + beg)) {
    setChildForward(origin + beg + 0, origin + mid);
  } else {
    setChildReverse(origin + beg + 2, origin + mid);
  }
}

// Specialized sort for 1 symbol + 1 delimiter.
// E.g. wildcard positions in spaced seeds.
void SubsetSuffixArray::radixSort1(std::vector<Range> &rangeStack,
				   const uchar *text, const uchar *subsetMap,
				   size_t beg, size_t end, size_t depth) {
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;

  size_t end0 = beg;  // end of '0's
  size_t begN = end;  // beginning of delimiters

  do {
    const size_t x = getBits(bits, a, end0);
    size_t y;
    switch (subsetMap[text[x]]) {
    case 0:
      end0++;
      break;
    default:  // the delimiter subset
      y = getBits(bits, a, --begN);
      setBits(bits, a, begN, x);
      setBits(bits, a, end0, y);
      break;
    }
  } while (end0 < begN);

  pushRange( rangeStack, beg, end0, depth );   // the '0's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
  }
}

// Specialized sort for 2 symbols + 1 delimiter.
// E.g. transition-constrained positions in subset seeds.
void SubsetSuffixArray::radixSort2(std::vector<Range> &rangeStack,
				   const uchar *text, const uchar *subsetMap,
				   size_t beg, size_t end, size_t depth) {
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;

  size_t end0 = beg;  // end of '0's
  size_t end1 = beg;  // end of '1's
  size_t begN = end;  // beginning of delimiters

  do {
    const size_t x = getBits(bits, a, end1);
    size_t y;
    switch (subsetMap[text[x]]) {
      case 0:
	y = getBits(bits, a, end0);
	setBits(bits, a, end0++, x);
	setBits(bits, a, end1++, y);
        break;
      case 1:
        end1++;
        break;
      default:  // the delimiter subset
	y = getBits(bits, a, --begN);
	setBits(bits, a, begN, x);
	setBits(bits, a, end1, y);
        break;
    }
  } while (end1 < begN);

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

// Specialized sort for 3 symbols + 1 delimiter.
// E.g. subset seeds for bisulfite-converted DNA.
void SubsetSuffixArray::radixSort3(std::vector<Range> &rangeStack,
				   const uchar *text, const uchar *subsetMap,
				   size_t beg, size_t end, size_t depth) {
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;

  size_t end0 = beg;  // end of '0's
  size_t end1 = beg;  // end of '1's
  size_t beg2 = end;  // beginning of '2's
  size_t begN = end;  // beginning of delimiters

  do {
    const size_t x = getBits(bits, a, end1);
    size_t y, z;
    switch (subsetMap[text[x]]) {
      case 0:
	y = getBits(bits, a, end0);
	setBits(bits, a, end0++, x);
	setBits(bits, a, end1++, y);
        break;
      case 1:
	end1++;
        break;
      case 2:
	y = getBits(bits, a, --beg2);
	setBits(bits, a, beg2, x);
	setBits(bits, a, end1, y);
        break;
      default:  // the delimiter subset
	y = getBits(bits, a, --beg2);
	z = getBits(bits, a, --begN);
	setBits(bits, a, end1, y);
	setBits(bits, a, beg2, z);
	setBits(bits, a, begN, x);
        break;
    }
  } while (end1 < beg2);

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's
  pushRange( rangeStack, beg2, begN, depth );  // the '2's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
    if( begN == end ) return;
    setChildForward( beg2, begN );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( beg2 == beg ) return;
    setChildReverse( begN, beg2 );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

// Specialized sort for 4 symbols + 1 delimiter.  E.g. DNA.
void SubsetSuffixArray::radixSort4(std::vector<Range> &rangeStack,
				   const uchar *text, const uchar *subsetMap,
				   size_t beg, size_t end, size_t depth) {
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;

  size_t end0 = beg;  // end of '0's
  size_t end1 = beg;  // end of '1's
  size_t end2 = beg;  // end of '2's
  size_t beg3 = end;  // beginning of '3's
  size_t begN = end;  // beginning of delimiters

  do {
    const size_t x = getBits(bits, a, end2);
    size_t y, z;
    switch (subsetMap[text[x]]) {
    case 0:
      y = getBits(bits, a, end1);
      z = getBits(bits, a, end0);
      setBits(bits, a, end2++, y);
      setBits(bits, a, end1++, z);
      setBits(bits, a, end0++, x);
      break;
    case 1:
      y = getBits(bits, a, end1);
      setBits(bits, a, end1++, x);
      setBits(bits, a, end2++, y);
      break;
    case 2:
      end2++;
      break;
    case 3:
      y = getBits(bits, a, --beg3);
      setBits(bits, a, beg3, x);
      setBits(bits, a, end2, y);
      break;
    default:  // the delimiter subset
      y = getBits(bits, a, --beg3);
      z = getBits(bits, a, --begN);
      setBits(bits, a, end2, y);
      setBits(bits, a, beg3, z);
      setBits(bits, a, begN, x);
      break;
    }
  } while (end2 < beg3);

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's
  pushRange( rangeStack, end1, end2, depth );  // the '2's
  pushRange( rangeStack, beg3, begN, depth );  // the '3's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
    if( end2 == end ) return;
    setChildForward( end1, end2 );
    if( begN == end ) return;
    setChildForward( beg3, begN );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( beg3 == beg ) return;
    setChildReverse( begN, beg3 );
    if( end1 == beg ) return;
    setChildReverse( end2, end1 );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

const int numOfBuckets = 2048;

void SubsetSuffixArray::radixSortN(std::vector<Range> &rangeStack,
				   const uchar *text, const uchar *subsetMap,
				   size_t beg, size_t end, size_t depth,
				   unsigned subsetCount, size_t *bucketSizes) {
  const bool isChildFwd = isChildDirectionForward(beg);
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;
  size_t whereTo[numOfBuckets];

  // get bucket sizes (i.e. letter counts):
  // The intermediate oracle array makes it faster (see "Engineering
  // Radix Sort for Strings" by J Karkkainen & T Rantala)
  for (size_t i = beg; i < end; /* noop */) {
    uchar oracle[256];
    size_t iEnd = i + std::min(sizeof(oracle), size_t(end - i));
    uchar *j = oracle;
    while (i < iEnd) {
      *j++ = subsetMap[text[getBits(bits, a, i++)]];
    }
    for (uchar *k = oracle; k < j; ++k) {
      ++bucketSizes[*k];
    }
  }

  size_t pos = beg;
  for (unsigned i = 0; i < subsetCount; ++i) {
    pos += bucketSizes[i];
    whereTo[i] = pos;
  }
  whereTo[CyclicSubsetSeed::DELIMITER] = end;

  // permute items into the correct buckets:
  for (size_t i = beg; i < pos; /* noop */) {
    size_t position = getBits(bits, a, i);
    unsigned subset;  // unsigned is faster than uchar!
    while (1) {
      subset = subsetMap[text[position]];
      size_t j = --whereTo[subset];
      if (j <= i) break;
      size_t p = getBits(bits, a, j);
      setBits(bits, a, j, position);
      position = p;
    }
    setBits(bits, a, i, position);
    size_t j = i + bucketSizes[subset];
    bucketSizes[subset] = 0;  // reset it so we can reuse it
    pushRange(rangeStack, i, j, depth);
    setChildLink(isChildFwd, 0, beg, end, i, j);
    i = j;
  }

  setChildLink(isChildFwd, 0, beg, end, pos, end);
}

static bool isMore(const size_t *steps, int minDelimDepth, size_t value) {
  int d = 0;
  while (value) {
    if (d >= minDelimDepth && value == steps[d] - 1) return false;
    ++d;
    value %= steps[d];
  }
  return true;
}

void SubsetSuffixArray::deepSort(std::vector<Range> &rangeStack,
				 const uchar *text, unsigned wordLength,
				 const CyclicSubsetSeed &seed,
				 const uchar *subsetMap,
				 size_t beg, size_t end, size_t depth,
				 size_t *bucketSizes) {
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;
  size_t whereTo[numOfBuckets];
  size_t steps[64];
  size_t depthEnd = maxBucketDepth(seed, depth, numOfBuckets, wordLength);
  size_t newDepth = depthEnd + seeds.size() - 1;
  int depthInc = depthEnd - depth;

  int minDelimDepth =
    (depth < wordLength) ? wordLength - depth : (depth < 1) ? 1 : 0;

  makeBucketSteps(steps, seed, depth, depthEnd, wordLength);

  for (size_t i = beg; i < end; ++i) {
    size_t position = getBits(bits, a, i);
    size_t v = bucketValue(seed, subsetMap, steps, text + position, depthInc);
    ++bucketSizes[v];
  }

  size_t pos = beg;
  for (size_t i = 0; i < steps[0]; ++i) {
    pos += bucketSizes[i];
    whereTo[i] = pos;
  }

  for (size_t i = beg; i < end; ) {
    size_t position = getBits(bits, a, i);
    size_t v;
    while (1) {
      v = bucketValue(seed, subsetMap, steps, text + position, depthInc);
      size_t j = --whereTo[v];
      if (j <= i) break;
      size_t p = getBits(bits, a, j);
      setBits(bits, a, j, position);
      position = p;
    }
    setBits(bits, a, i, position);
    size_t j = i + bucketSizes[v];
    bucketSizes[v] = 0;  // reset it so we can reuse it
    if (isMore(steps, minDelimDepth, v)) pushRange(rangeStack, i, j, newDepth);
    i = j;
  }
}

void SubsetSuffixArray::twoArraySort(std::vector<Range> &rangeStack,
				     const uchar *text, const uchar *subsetMap,
				     size_t origin, size_t beg, size_t end,
				     size_t depth, unsigned subsetCount,
				     size_t cacheSize,
				     size_t *positions, uchar *seqCache) {
  const bool isChildFwd = isChildDirectionForward(origin + beg);
  size_t whereTo[numOfBuckets];
  size_t *bucketSizes = positions + cacheSize;
  size_t *positions2 = bucketSizes + numOfBuckets;

  for (size_t i = beg; i < end; ++i) {
    seqCache[i] = subsetMap[text[positions[i]]];
  }  // this loop fission seems to make it much faster
  for (size_t i = beg; i < end; ++i) {
    ++bucketSizes[seqCache[i]];
  }

  size_t oldPos = beg;
  for (unsigned i = 0; i < subsetCount; ++i) {
    size_t newPos = oldPos + bucketSizes[i];
    bucketSizes[i] = 0;  // reset it so we can reuse it
    whereTo[i] = oldPos;
    pushRange(rangeStack, oldPos, newPos, depth);
    setChildLink(isChildFwd, origin, beg, end, oldPos, newPos);
    oldPos = newPos;
  }
  // don't sort within the delimiter bucket:
  bucketSizes[CyclicSubsetSeed::DELIMITER] = 0;
  whereTo[CyclicSubsetSeed::DELIMITER] = oldPos;
  setChildLink(isChildFwd, origin, beg, end, oldPos, end);

  for (size_t i = beg; i < end; ++i) {
    positions2[whereTo[seqCache[i]]++] = positions[i];
  }
  memcpy(positions + beg, positions2 + beg, (end - beg) * sizeof(size_t));
}

static size_t rangeSize(const Range &r) {
  return r.end - r.beg;
}

static size_t nextRangeSize(const std::vector<Range> &ranges) {
  return rangeSize(ranges.back());
}

static size_t rangeSizeSum(const std::vector<Range> &ranges) {
  size_t s = 0;
  for (size_t i = 0; i < ranges.size(); ++i) {
    s += rangeSize(ranges[i]);
  }
  return s;
}

static size_t numOfThreadsForOneRange(size_t numOfThreads,
				      size_t sizeOfThisRange,
				      size_t sizeOfAllRanges) {
  // We want:
  // min(t | sizeOfThisRange / t < sizeOfOtherRanges / (numOfThreads - (t+1)))
  // Or equivalently:
  // max(t | sizeOfThisRange / (t-1) >= sizeOfOtherRanges / (numOfThreads - t))
  double x = numOfThreads - 1;  // use double to avoid overflow
  return (x * sizeOfThisRange + sizeOfAllRanges) / sizeOfAllRanges;
}

static unsigned getSubsetCount(const CyclicSubsetSeed &seed,
			       size_t depth, unsigned wordLength) {
  return (depth < wordLength)
    ? seed.restrictedSubsetCount(depth)  // xxx inefficient
    : seed.unrestrictedSubsetCount(depth);
}

void SubsetSuffixArray::sortOutOfPlace(std::vector<Range> &stack,
				       size_t cacheSize,
				       size_t *intCache, uchar *seqCache,
				       const uchar *text, unsigned wordLength,
				       const CyclicSubsetSeed &seed,
				       size_t maxUnsortedInterval,
				       size_t origin) {
  size_t stackBase = stack.size();

  while (stack.size() >= stackBase) {
    size_t beg = stack.back().beg;
    size_t end = stack.back().end;
    size_t depth = stack.back().depth;
    stack.pop_back();
    size_t interval = end - beg;

    const size_t minLength = 1;
    if (interval <= maxUnsortedInterval && depth >= minLength) continue;

    const uchar *textBase = text + depth;
    const uchar *subsetMap = seed.subsetMap(depth);

    if (childTable.v.empty() && kiddyTable.v.empty() && chibiTable.v.empty()) {
      if (interval < 10) {  // ???
	insertionSort(textBase, seed, intCache, beg, end, subsetMap);
	continue;
      }
    } else {
      if (interval == 2) {
	sort2(textBase, seed, intCache, origin, beg, subsetMap);
	continue;
      }
    }

    unsigned subsetCount = getSubsetCount(seed, depth, wordLength);
    twoArraySort(stack, textBase, subsetMap, origin, beg, end,
		 depth + 1, subsetCount, cacheSize, intCache, seqCache);
  }
}

void SubsetSuffixArray::sortRanges(std::vector<Range> *stacks,
				   size_t cacheSize,
				   size_t *intCache, uchar *seqCache,
				   const uchar *text, unsigned wordLength,
				   unsigned seedNum,
				   size_t maxUnsortedInterval,
				   size_t numOfThreads) {
  const size_t numOfSeeds = seeds.size();
  size_t *a = (size_t *)&suffixArray.v[0];
  const int bits = sufArray.bitsPerItem;
  std::vector<Range> &myStack = stacks[0];

  while( !myStack.empty() ){
#ifdef HAS_CXX_THREADS
    size_t numOfChunks = std::min(numOfThreads, myStack.size());
    if (numOfChunks > 1) {
      size_t intCacheSize = cacheSize * 2 + numOfBuckets;
      size_t totalSize = rangeSizeSum(myStack);
      size_t numOfNewThreads = numOfChunks - 1;

      size_t thisSize = nextRangeSize(myStack);
      size_t t = numOfThreadsForOneRange(numOfThreads, thisSize, totalSize);
      size_t maxThreads = numOfThreads - numOfNewThreads;
      size_t thisThreads = std::min(t, maxThreads);
      numOfThreads -= thisThreads;
      size_t pad = numOfItemsBetweenWrites(bits) * numOfThreads;

      do {
	totalSize -= nextRangeSize(myStack);
	myStack.back().beg += pad;
	myStack.back().end += pad;
	stacks[numOfThreads].push_back(myStack.back());
	myStack.pop_back();
	thisSize += nextRangeSize(myStack);
      } while (myStack.size() > numOfThreads &&
	       thisSize <= totalSize / numOfThreads);
      // We want:
      // max(r | sizeSum(r) <= (totalSize - sizeSum(r-1)) / newNumOfThreads)

      size_t beg = stacks[numOfThreads].back().beg;
      size_t end = stacks[numOfThreads].front().end;

      for (size_t i = end; i-- > beg;) {
	setBits(bits, a, i, getBits(bits, a, i - pad));
      }

      std::thread myThread(&SubsetSuffixArray::sortRanges, this,
			   stacks + numOfThreads, cacheSize,
			   intCache + numOfThreads * intCacheSize,
			   seqCache + numOfThreads * cacheSize,
			   text, wordLength, seedNum,
			   maxUnsortedInterval, thisThreads);
      sortRanges(stacks, cacheSize, intCache, seqCache, text, wordLength,
		 seedNum, maxUnsortedInterval, numOfThreads);
      myThread.join();

      for (size_t i = beg; i < end; ++i) {
	setBits(bits, a, i - pad, getBits(bits, a, i));
      }

      // unnecessary: aims to make all output bits unchanged by multi-threading
      for (size_t i = end - pad; i < end; ++i) setBits(bits, a, i, 0);

      return;
    }
#endif

    size_t beg = myStack.back().beg;
    size_t end = myStack.back().end;
    size_t depth = myStack.back().depth;
    myStack.pop_back();

    if (depth < numOfSeeds) {
      seedNum = depth;
      depth = 0;
    } else {
      depth -= numOfSeeds - 1;
    }

    const CyclicSubsetSeed &seed = seeds[seedNum];
    size_t interval = end - beg;

    if (interval <= cacheSize) {
      unpackBits(bits, a, intCache, beg, end);
      pushRange(myStack, 0, interval, depth);
      sortOutOfPlace(myStack, cacheSize, intCache, seqCache, text,
		     wordLength, seed, maxUnsortedInterval, beg);
      packBits(bits, a, intCache, beg, end);
      continue;
    }

    const uchar *textBase = text + depth;
    const uchar *subsetMap = seed.subsetMap(depth);

    if (childTable.v.empty() && kiddyTable.v.empty() && chibiTable.v.empty()) {
      deepSort(myStack, textBase, wordLength, seed, subsetMap, beg, end, depth,
	       intCache + cacheSize);
      continue;
    }

    unsigned subsetCount = getSubsetCount(seed, depth, wordLength);
    depth += numOfSeeds;

    switch (subsetCount) {
    case 1:  radixSort1(myStack, textBase, subsetMap, beg, end, depth); break;
    case 2:  radixSort2(myStack, textBase, subsetMap, beg, end, depth); break;
    case 3:  radixSort3(myStack, textBase, subsetMap, beg, end, depth); break;
    case 4:  radixSort4(myStack, textBase, subsetMap, beg, end, depth); break;
    default: radixSortN(myStack, textBase, subsetMap, beg, end, depth,
			subsetCount, intCache + cacheSize);
    }
  }
}

void SubsetSuffixArray::sortIndex( const uchar* text,
				   unsigned wordLength,
				   const size_t* cumulativeCounts,
				   size_t maxUnsortedInterval,
				   int childTableType,
				   size_t numOfThreads ){
  size_t numOfSeeds = seeds.size();
  size_t total = cumulativeCounts[numOfSeeds - 1];
  size_t cacheSize = total / (32 * sizeof(size_t)) / numOfThreads;
  // tie-breaking depends on cacheSize, so it depends on numOfThreads

  if (childTableType == 1) chibiTable.v.assign(total, -1);
  if (childTableType == 2) kiddyTable.v.assign(total, -1);
  if (childTableType == 3) {
    childTable.v.assign(numOfBytes(chiArray.bitsPerItem, total), 0);
    chiArray.items = (const size_t *)childTable.begin();
  }

  std::vector< std::vector<Range> > stacks(numOfThreads);
  std::vector<size_t> intCache(numOfThreads * (cacheSize * 2 + numOfBuckets));
  std::vector<uchar> seqCache(numOfThreads * cacheSize);

  size_t beg = 0;
  for (size_t i = 0; i < numOfSeeds; ++i) {
    size_t end = cumulativeCounts[i];
    pushRange(stacks[0], beg, end, i);
    setChildReverse(end, beg);
    beg = end;
  }

  sortRanges(&stacks[0], cacheSize, &intCache[0], &seqCache[0], text,
	     wordLength, 0, maxUnsortedInterval, numOfThreads);
}
