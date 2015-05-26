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

#ifndef SUBSET_SUFFIX_ARRAY_HH
#define SUBSET_SUFFIX_ARRAY_HH

#include "CyclicSubsetSeed.hh"
#include "VectorOrMmap.hh"

namespace cbrc{

class SubsetSuffixArray{
public:
  typedef unsigned indexT;

  CyclicSubsetSeed& getSeed() { return seed; }
  const CyclicSubsetSeed& getSeed() const { return seed; }

  // Add every step-th text position in the range [beg,end).
  // Positions starting with delimiters aren't added.
  // The positions aren't sorted.
  void addPositions( const uchar* text, indexT beg, indexT end, indexT step );

  // Sort the suffix array (but don't make the buckets).
  void sortIndex( const uchar* text, indexT maxUnsortedInterval );

  // Make the buckets.  If bucketDepth+1 == 0, then a default
  // bucketDepth is used.  The default is: the maximum possible
  // bucketDepth such that the number of bucket entries is at most 1/4
  // the number of suffix array entries.
  void makeBuckets( const uchar* text, indexT bucketDepth );

  // Clear the positions, so we can add new positions from scratch.
  void clearPositions();

  void fromFiles( const std::string& baseName,
		  bool isMaskLowercase, const uchar letterCode[] );

  void toFiles( const std::string& baseName,
		bool isAppendPrj, indexT textLength ) const;

  // Find the smallest match to the text, starting at the given
  // position in the query, such that there are at most maxHits
  // matches, and the match-depth is at least minDepth, or the
  // match-depth is maxDepth.  Return the range of matching indices
  // via begPtr and endPtr.
  void match( const indexT*& begPtr, const indexT*& endPtr,
              const uchar* queryPtr, const uchar* text,
              indexT maxHits, indexT minDepth, indexT maxDepth ) const;

  // Count matches of all sizes (up to maxDepth), starting at the
  // given position in the query.
  void countMatches( std::vector<unsigned long long>& counts,
		     const uchar* queryPtr, const uchar* text,
		     indexT maxDepth ) const;

private:
  CyclicSubsetSeed seed;
  VectorOrMmap<indexT> suffixArray;  // sorted indices
  VectorOrMmap<indexT> buckets;
  std::vector<indexT> bucketSteps;  // step size for each k-mer

  // These find the suffix array range of one letter, whose subset is
  // "subset", within the suffix array range [beg, end):
  void equalRange( indexT& beg, indexT& end, const uchar* textBase,
		   const uchar* subsetMap, uchar subset ) const;
  indexT lowerBound( indexT beg, indexT end, const uchar* textBase,
		     const uchar* subsetMap, uchar subset ) const;
  indexT upperBound( indexT beg, indexT end, const uchar* textBase,
		     const uchar* subsetMap, uchar subset ) const;

  // These find the suffix array range of string [queryBeg, queryEnd)
  // within the suffix array range [beg, end):
  void equalRange2( indexT& beg, indexT& end,
		    const uchar* queryBeg, const uchar* queryEnd,
		    const uchar* textBase, const uchar* subsetMap ) const;
  indexT lowerBound2( indexT beg, indexT end,
		      const uchar* queryBeg, const uchar* queryEnd,
		      const uchar* textBase, const uchar* subsetMap ) const;
  indexT upperBound2( indexT beg, indexT end,
		      const uchar* queryBeg, const uchar* queryEnd,
		      const uchar* textBase, const uchar* subsetMap ) const;

  // Return the maximum prefix size covered by the buckets.
  indexT maxBucketPrefix() const { return bucketSteps.size() - 1; }

  indexT defaultBucketDepth();

  void makeBucketSteps( indexT bucketDepth );

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
};

}  // end namespace
#endif
