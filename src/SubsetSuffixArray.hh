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
#include "VectorOrMmap.hh"

#include <climits>

namespace cbrc{

typedef LAST_INT_TYPE PosPart;
const int posParts = 1;

typedef LAST_INT_TYPE OffPart;
const int offParts = 1;

inline size_t posGet(const PosPart *p) {
  size_t x = 0;
  for (int i = 0; i < posParts; ++i) {
    size_t y = p[i];  // must convert to size_t before shifting!
    x += y << (i * sizeof(PosPart) * CHAR_BIT);
  }
  return x;
}

inline void posSet(PosPart *p, size_t value) {
  for (int i = 0; i < posParts; ++i) {
    p[i] = value >> (i * sizeof(PosPart) * CHAR_BIT);
  }
}

inline size_t posCount(const PosPart *beg, const PosPart *end) {
  size_t d = end - beg;
  return d / posParts;  // faster if the dividend is unsigned?
}

class SubsetSuffixArray{
public:
  typedef LAST_INT_TYPE indexT;

  struct Range {PosPart *beg; PosPart *end; indexT depth;};

  std::vector<CyclicSubsetSeed> &getSeeds() { return seeds; }
  const std::vector<CyclicSubsetSeed> &getSeeds() const { return seeds; }

  size_t size() const { return suffixArray.size() / posParts; }

  // Add every step-th text position in the range [beg,end).
  // Positions starting with delimiters aren't added.
  // The positions aren't sorted.
  // If minimizerWindow > 1 then the positions are added only if they
  // are "minimizers" for the given window and seed pattern.
  void addPositions( const uchar* text, indexT beg, indexT end,
		     size_t step, size_t minimizerWindow );

  // Store positions in [seqBeg, seqEnd) where certain "words" start.
  // The cumulative word counts must be provided.  (cumulativeCounts
  // is internally modified and restored to its original values).
  void setWordPositions(const DnaWordsFinder &finder, size_t *cumulativeCounts,
			const uchar *seqBeg, const uchar *seqEnd);

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
		   size_t minPositionsPerBucket, unsigned bucketDepth);

  void fromFiles(const std::string &baseName,
		 bool isMaskLowercase, const uchar letterCode[],
		 const std::string &mainSequenceAlphabet);

  void toFiles( const std::string& baseName,
		bool isAppendPrj, size_t textLength ) const;

  // Find the smallest match to the text, starting at the given
  // position in the query, such that there are at most maxHits
  // matches, and the match-depth is at least minDepth, or the
  // match-depth is maxDepth.  Return the range of matching indices
  // via begPtr and endPtr.
  void match( const PosPart *&begPtr, const PosPart *&endPtr,
              const uchar* queryPtr, const uchar* text, unsigned seedNum,
              size_t maxHits, size_t minDepth, size_t maxDepth ) const;

  // Count matches of all sizes (up to maxDepth), starting at the
  // given position in the query.
  void countMatches( std::vector<unsigned long long>& counts,
		     const uchar* queryPtr, const uchar* text,
		     unsigned seedNum, size_t maxDepth ) const;

private:
  std::vector<CyclicSubsetSeed> seeds;
  std::vector<const OffPart *> bucketEnds;
  std::vector<const indexT *> bucketStepEnds;

  VectorOrMmap<PosPart> suffixArray;  // sorted indices
  VectorOrMmap<OffPart> buckets;
  std::vector<indexT> bucketSteps;  // step size for each k-mer

  VectorOrMmap<indexT> childTable;
  VectorOrMmap<unsigned short> kiddyTable;  // smaller child table
  VectorOrMmap<unsigned char> chibiTable;  // even smaller child table

  enum ChildDirection { FORWARD, REVERSE, UNKNOWN };

  // These find the suffix array range of one letter, whose subset is
  // "subset", within the suffix array range [beg, end):
  void equalRange( indexT& beg, indexT& end, const uchar* textBase,
		   const uchar* subsetMap, uchar subset ) const;
  indexT lowerBound( indexT beg, indexT end, const uchar* textBase,
		     const uchar* subsetMap, uchar subset ) const;
  indexT upperBound( indexT beg, indexT end, const uchar* textBase,
		     const uchar* subsetMap, uchar subset ) const;

  // This does the same thing as equalRange, but uses a child table:
  void childRange( indexT& beg, indexT& end, ChildDirection& childDirection,
                   const uchar* textBase,
                   const uchar* subsetMap, uchar subset ) const;

  // These find the suffix array range of string [queryBeg, queryEnd)
  // within the suffix array range [beg, end):
  void equalRange2( indexT& beg, indexT& end,
		    const uchar* queryBeg, const uchar* queryEnd,
		    const uchar* textBase, const CyclicSubsetSeed& seed,
		    const uchar* subsetMap ) const;
  indexT lowerBound2( indexT beg, indexT end,
		      const uchar* queryBeg, const uchar* queryEnd,
		      const uchar* textBase, const CyclicSubsetSeed& seed,
		      const uchar* subsetMap ) const;
  indexT upperBound2( indexT beg, indexT end,
		      const uchar* queryBeg, const uchar* queryEnd,
		      const uchar* textBase, const CyclicSubsetSeed& seed,
		      const uchar* subsetMap ) const;

  // Return the maximum prefix size covered by the buckets.
  size_t maxBucketPrefix(unsigned seedNum) const
  { return bucketStepEnds[seedNum + 1] - bucketStepEnds[seedNum] - 1; }

  void makeBucketSteps(const unsigned *bucketDepths, size_t wordLength);

  size_t bucketsSize() const {
    size_t n = offParts;
    for (size_t i = 0; i < seeds.size(); ++i) {
      n += bucketStepEnds[i][0];
    }
    return n;
  }

  void initBucketEnds() {
    bucketEnds.resize(seeds.size());
    const OffPart *p = &buckets[0];
    for (size_t i = 0; i < seeds.size(); ++i) {
      bucketEnds[i] = p;
      p += bucketStepEnds[i][0];
    }
  }

  void sort2( const uchar* text, const CyclicSubsetSeed& seed,
	      PosPart *beg, const uchar* subsetMap );

  void radixSort1( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSort2( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSort3( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSort4( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSortN( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth,
		   unsigned subsetCount, indexT* bucketSize );

  void sortRanges( std::vector<Range>* stacks, indexT* bucketSizes,
		   const uchar* text,
		   unsigned wordLength, const CyclicSubsetSeed& seed,
		   size_t maxUnsortedInterval, size_t numOfThreads );

  // Same as the 1st equalRange, but uses more info and may be faster:
  void equalRange( indexT& beg, indexT& end, const uchar* textBase,
                   const uchar* subsetMap, uchar subset,
                   uchar begSubset, uchar endSubset,
                   indexT begOffset, indexT endOffset ) const{
    if( subset == begSubset ){
      end = upperBound( beg + begOffset, end - endOffset,
                        textBase, subsetMap, subset );
    }else if( subset == endSubset ){
      beg = lowerBound( beg + begOffset, end - endOffset,
                        textBase, subsetMap, subset );
    }else{
      beg += begOffset;
      end -= endOffset;
      equalRange( beg, end, textBase, subsetMap, subset );
    }
  }

  // Same as the 1st equalRange, but tries to be faster by checking endpoints:
  void fastEqualRange( indexT& beg, indexT& end, const uchar* textBase,
                       const uchar* subsetMap, uchar subset ) const{
    uchar b = subsetMap[ textBase[ suffixArray[ beg ] ] ];
    if( subset < b ){ end = beg; return; }
    uchar e = subsetMap[ textBase[ suffixArray[ end - 1 ] ] ];
    if( subset > e ){ beg = end; return; }
    if( b == e ) return;
    equalRange( beg, end, textBase, subsetMap, subset, b, e, 1, 1 );
  }

  indexT getChildForward( indexT from ) const{
    return
      !childTable.empty() ? childTable[ from ] :
      !kiddyTable.empty() ? from + kiddyTable[ from ] :
      !chibiTable.empty() ? from + chibiTable[ from ] : from;
  }

  indexT getChildReverse( indexT from ) const{
    return
      !childTable.empty() ? childTable[ from - 1 ] :
      !kiddyTable.empty() ? from - kiddyTable[ from - 1 ] :
      !chibiTable.empty() ? from - chibiTable[ from - 1 ] : from;
  }

  void setKiddy( indexT index, indexT value ){
    kiddyTable.v[ index ] = (value < USHRT_MAX) ? value : 0;
  }

  void setChibi( indexT index, indexT value ){
    chibiTable.v[ index ] = (value < UCHAR_MAX) ? value : 0;
  }

  void setChildForward(const PosPart *from, const PosPart *to) {
    if( to == from ) return;
    const PosPart *origin = &suffixArray.v[0];
    indexT i = posCount(origin, from);
    /**/ if (!childTable.v.empty()) childTable.v[i] = posCount(origin, to);
    else if (!kiddyTable.v.empty()) setKiddy(i, posCount(from, to));
    else if (!chibiTable.v.empty()) setChibi(i, posCount(from, to));
  }

  void setChildReverse(const PosPart *from, const PosPart *to) {
    if( to == from ) return;
    const PosPart *origin = &suffixArray.v[0];
    indexT i = posCount(origin, from) - 1;
    /**/ if (!childTable.v.empty()) childTable.v[i] = posCount(origin, to);
    else if (!kiddyTable.v.empty()) setKiddy(i, posCount(to, from));
    else if (!chibiTable.v.empty()) setChibi(i, posCount(to, from));
  }

  bool isChildDirectionForward(const PosPart *beg) const {
    const PosPart *origin = &suffixArray.v[0];
    indexT i = posCount(origin, beg);
    return
      !childTable.v.empty() ? childTable.v[ i ] == 0 :
      !kiddyTable.v.empty() ? kiddyTable.v[ i ] == USHRT_MAX :
      !chibiTable.v.empty() ? chibiTable.v[ i ] == UCHAR_MAX : true;
  }
};

}  // end namespace
#endif
