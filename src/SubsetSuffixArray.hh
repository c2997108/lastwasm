// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

// This class holds a suffix array.  The suffix array is just a list
// of numbers indicating positions in a text, sorted according to the
// alphabetical order of the text suffixes starting at these
// positions.  A query sequence can then be matched incrementally to
// the suffix array using binary search.

// A "subset suffix array" means that, when comparing two suffixes, we
// consider subsets of letters to be equivalent.  For example, we
// might consider purines to be equivalent to each other, and
// pyrimidines to be equivalent to each other.  The subsets may vary
// from position to position as we compare two suffixes.

// There is always a special subset, called DELIMITER, which doesn't
// match anything.

// For faster matching, we use "buckets", which store the start and
// end in the suffix array of all size-k prefixes of the suffixes.
// They store this information for all values of k from 1 to, say, 12.

// This class can store multiple concatenated suffix arrays: each
// array holds suffixes starting with each pattern (of length
// "wordLength") in a DnaWordsFinder.  The endpoints of the
// concatenated arrays are specified by "cumulativeCounts".  Each
// array has its own letter-subsets.  "seedNum" specifies one of the
// arrays.

#ifndef SUBSET_SUFFIX_ARRAY_HH
#define SUBSET_SUFFIX_ARRAY_HH

#include "CyclicSubsetSeed.hh"
#include "dna_words_finder.hh"
#include "mcf_big_seq.hh"
#include "mcf_packed_array.hh"
#include "VectorOrMmap.hh"

#include <climits>

namespace cbrc {

using namespace mcf;

inline size_t numOfBytes(int bitsPerItem, size_t numOfItems) {
  return numOfWordsNeededFor(bitsPerItem, numOfItems) * sizeof(size_t);
}

inline size_t totalItemsBetweenThreads(int bitsPerItem, size_t numOfThreads) {
  return numOfItemsBetweenWrites(bitsPerItem) * (numOfThreads - 1);
}

inline size_t getItem(ConstPackedArray a, size_t i) {
  if ((a.bitsPerItem & 128) == 0) return a[i];
  // Stuff below here is just for reading old lastdb databases:
  if (a.bitsPerItem == 128 + 32) {
    const unsigned *b = (const unsigned *)a.items;
    return b[i];
  }
  const unsigned char *c = (const unsigned char *)a.items;
  return (size_t)c[i*5]
    |    (size_t)c[i*5 + 1] << 8
    |    (size_t)c[i*5 + 2] << 16
    |    (size_t)c[i*5 + 3] << 24
    |    (size_t)c[i*5 + 4] << 32;
}

inline size_t maxBucketDepth(const CyclicSubsetSeed &seed, size_t startDepth,
			     size_t maxBuckets, unsigned wordLength) {
  unsigned long long numOfBuckets = (startDepth > 0);  // delimiter if depth>0
  unsigned long long product = 1;
  for (size_t d = startDepth; ; ++d) {
    if (d < wordLength) {
      product *= seed.restrictedSubsetCount(d);
      numOfBuckets = product;
    } else {
      product *= seed.unrestrictedSubsetCount(d);
      numOfBuckets += product;
    }
    if (numOfBuckets > maxBuckets) return d;
  }
}

inline void makeBucketSteps(size_t *steps, const CyclicSubsetSeed &seed,
			    size_t depthBeg, size_t depthEnd,
			    size_t wordLength) {
  size_t step = 1;
  steps[depthEnd - depthBeg] = step;

  while (depthEnd-- > depthBeg) {
    step = (depthEnd < wordLength)
      ? step * seed.restrictedSubsetCount(depthEnd)
      : step * seed.unrestrictedSubsetCount(depthEnd) + (depthEnd > 0);
    // add one for delimiters, except when depth==0
    steps[depthEnd - depthBeg] = step;
  }
}

inline size_t bucketValue(const CyclicSubsetSeed &seed, const uchar *subsetMap,
			  const size_t *steps, const uchar *text, int depth) {
  size_t val = 0;
  for (int d = 0; d < depth; ) {
    unsigned s = subsetMap[text[d]];
    if (s == CyclicSubsetSeed::DELIMITER) return val + steps[d] - 1;
    ++d;
    val += s * steps[d];
    subsetMap = seed.nextMap(subsetMap);
  }
  return val;
}

class SubsetSuffixArray {
public:
  struct Range {size_t beg; size_t end; size_t depth;};

  std::vector<CyclicSubsetSeed> &getSeeds() { return seeds; }
  const std::vector<CyclicSubsetSeed> &getSeeds() const { return seeds; }

  void resizePositions(size_t numOfPositions, size_t seqLength,
		       size_t numOfThreads) {
    sufArray.bitsPerItem = numOfBitsNeededFor(seqLength - 1);
    bckArray.bitsPerItem = numOfBitsNeededFor(numOfPositions);
    chiArray.bitsPerItem = numOfBitsNeededFor(numOfPositions - 1);
    size_t n = totalItemsBetweenThreads(sufArray.bitsPerItem, numOfThreads);
    suffixArray.v.resize(numOfBytes(sufArray.bitsPerItem, numOfPositions + n));
    sufArray.items = (const size_t *)suffixArray.begin();
  }

  // Get the i-th item in the suffix array
  size_t getPosition(size_t i) const {
    return getItem(sufArray, i);
  }

  // Set the i-th item of the suffix array to x
  void setPosition(size_t i, size_t x) {
    setBits(sufArray.bitsPerItem, (size_t *)&suffixArray.v[0], i, x);
  }

  // Store positions in [seqBeg, seqEnd) where certain "words" start.
  // The cumulative word counts must be provided.  (cumulativeCounts
  // is internally modified and restored to its original values).
  void setWordPositions(const DnaWordsFinder &finder, size_t *cumulativeCounts,
			const uchar *seqBeg, const uchar *seqEnd,
			size_t numOfThreads);

  // Sort the suffix array (but don't make the buckets).
  void sortIndex(const uchar *text,
		 unsigned wordLength, const size_t *cumulativeCounts,
		 size_t maxUnsortedInterval, int childTableType,
		 size_t numOfThreads);

  // Make the buckets.  If bucketDepth+1 == 0, then the bucket depth
  // is: the maximum possible such that (memory use of buckets) <=
  // (memory use of stored positions) / minPositionsPerBucket.
  void makeBuckets(const uchar *text,
		   unsigned wordLength, const size_t *cumulativeCounts,
		   size_t minPositionsPerBucket, unsigned bucketDepth,
		   size_t numOfThreads);

  void fromFiles(const std::string &baseName, int bitsPerInt,
		 bool isMaskLowercase, const uchar letterCode[],
		 const std::string &mainSequenceAlphabet);

  void toFiles( const std::string& baseName,
		bool isAppendPrj, size_t textLength ) const;

  // Find the smallest match to the text, starting at the given
  // position in the query, such that there are at most maxHits
  // matches, and the match-depth is at least minDepth, or the
  // match-depth is maxDepth.  Return the range of matching indices
  // via beg and end.
  void match(size_t &beg, size_t &end,
	     const uchar *queryPtr, BigSeq text, unsigned seedNum,
	     size_t maxHits, size_t minDepth, size_t maxDepth) const;

  // Count matches of all sizes (up to maxDepth), starting at the
  // given position in the query.
  void countMatches( std::vector<unsigned long long>& counts,
		     const uchar *queryPtr, BigSeq text,
		     unsigned seedNum, size_t maxDepth ) const;

private:
  std::vector<CyclicSubsetSeed> seeds;
  std::vector<size_t> bucketEnds;
  std::vector<size_t> bucketSteps;  // step size for each k-mer
  std::vector<const size_t *> bucketStepEnds;

  // Use uchar instead of size_t just to support old lastdb databases:
  VectorOrMmap<uchar> suffixArray;  // sorted positions
  VectorOrMmap<uchar> buckets;

  VectorOrMmap<uchar> childTable;
  VectorOrMmap<unsigned short> kiddyTable;  // smaller child table
  VectorOrMmap<unsigned char> chibiTable;  // even smaller child table

  ConstPackedArray sufArray;
  ConstPackedArray bckArray;
  ConstPackedArray chiArray;

  enum ChildDirection { FORWARD, REVERSE, UNKNOWN };

  // This does the same thing as equalRange, but uses a child table:
  void childRange(size_t &beg, size_t &end, ChildDirection &childDirection,
		  BigPtr textBase, const uchar *subsetMap, uchar subset) const;

  // Return the maximum prefix size covered by the buckets.
  size_t maxBucketPrefix(unsigned seedNum) const
  { return bucketStepEnds[seedNum + 1] - bucketStepEnds[seedNum] - 1; }

  void makeBucketStepsAndEnds(const unsigned *bucketDepths, size_t wordLength);

  size_t bucketsSize() const { return bucketEnds.back() + 1; }

  void sort2(const uchar *text, const CyclicSubsetSeed &seed,
	     size_t *positions, size_t origin,
	     size_t beg, const uchar *subsetMap);

  void radixSort1(std::vector<Range> &rangeStack,
		  const uchar *text, const uchar *subsetMap,
		  size_t beg, size_t end, size_t depth);
  void radixSort2(std::vector<Range> &rangeStack,
		  const uchar *text, const uchar *subsetMap,
		  size_t beg, size_t end, size_t depth);
  void radixSort3(std::vector<Range> &rangeStack,
		  const uchar *text, const uchar *subsetMap,
		  size_t beg, size_t end, size_t depth);
  void radixSort4(std::vector<Range> &rangeStack,
		  const uchar *text, const uchar *subsetMap,
		  size_t beg, size_t end, size_t depth);
  void radixSortN(std::vector<Range> &rangeStack,
		  const uchar *text, const uchar *subsetMap,
		  size_t beg, size_t end, size_t depth,
		  unsigned subsetCount, size_t *bucketSizes);

  void deepSort(std::vector<Range> &rangeStack,
		const uchar *text, unsigned wordLength,
		const CyclicSubsetSeed &seed, const uchar *subsetMap,
		size_t beg, size_t end, size_t depth, size_t *bucketSizes);

  void twoArraySort(std::vector<Range> &rangeStack,
		    const CyclicSubsetSeed &seed,
		    const uchar *text, const uchar *subsetMap,
		    size_t origin, size_t beg, size_t end, size_t depth,
		    size_t cacheSize, size_t *positions, uchar *seqCache);

  void sortOutOfPlace(std::vector<Range> &stack, size_t cacheSize,
		      size_t *intCache, uchar *seqCache, const uchar *text,
		      const CyclicSubsetSeed &seed,
		      size_t maxUnsortedInterval, size_t origin);

  void sortRanges(std::vector<Range> *stacks, size_t cacheSize,
		  size_t *intCache, uchar *seqCache, const uchar *text,
		  unsigned wordLength, unsigned seedNum,
		  size_t maxUnsortedInterval, size_t numOfThreads,
		  size_t endOfRanges);

  size_t getChild(size_t i) const {
    if ((chiArray.bitsPerItem & 128) == 0) return chiArray[i];
    // Stuff below here is just for reading old lastdb databases:
    if (chiArray.bitsPerItem == 128 + 32) {
      const unsigned *b = (const unsigned *)chiArray.items;
      return b[i];
    }
    return chiArray.items[i];
  }

  size_t getChildForward(size_t from) const {
    return
      !childTable.empty() ? getChild(from) :
      !kiddyTable.empty() ? from + kiddyTable[from] :
      !chibiTable.empty() ? from + chibiTable[from] : from;
  }

  size_t getChildReverse(size_t from) const {
    return
      !childTable.empty() ? getChild(from - 1) :
      !kiddyTable.empty() ? from - kiddyTable[from - 1] :
      !chibiTable.empty() ? from - chibiTable[from - 1] : from;
  }

  void setChild(size_t index, size_t value) {
    setBits(chiArray.bitsPerItem, (size_t *)&childTable.v[0], index, value);
  }

  void setKiddy(size_t index, size_t value) {
    kiddyTable.v[index] = (value < USHRT_MAX) ? value : 0;
  }

  void setChibi(size_t index, size_t value) {
    chibiTable.v[index] = (value < UCHAR_MAX) ? value : 0;
  }

  void setChildForward(size_t from, size_t to) {
    if (to == from) return;
    /**/ if (!childTable.v.empty()) setChild(from, to);
    else if (!kiddyTable.v.empty()) setKiddy(from, to - from);
    else if (!chibiTable.v.empty()) setChibi(from, to - from);
  }

  void setChildReverse(size_t from, size_t to) {
    if (to == from) return;
    /**/ if (!childTable.v.empty()) setChild(from - 1, to);
    else if (!kiddyTable.v.empty()) setKiddy(from - 1, from - to);
    else if (!chibiTable.v.empty()) setChibi(from - 1, from - to);
  }

  void setChildLink(bool isFwd, size_t origin, size_t beg, size_t end,
		    size_t lo, size_t hi) {
    if (isFwd) {
      if (hi < end) setChildForward(origin + lo, origin + hi);
    } else {
      if (lo > beg) setChildReverse(origin + hi, origin + lo);
    }
  }

  bool isChildDirectionForward(size_t beg) const {
    return
      !childTable.v.empty() ? chiArray[beg] == 0 :
      !kiddyTable.v.empty() ? kiddyTable.v[beg] == USHRT_MAX :
      !chibiTable.v.empty() ? chibiTable.v[beg] == UCHAR_MAX : true;
  }
};

}

#endif
