// Copyright 2008, 2009 Martin C. Frith

// This class holds a suffix array.  The suffix array is just a list
// of numbers indicating positions in a text, sorted according to the
// alphabetical order of the text suffixes starting at these
// positions.  A query sequence can then be matched incrementally to
// the suffix array using binary search.

// For faster matching, we use "buckets", which store the start and
// end in the suffix array of all size-k prefixes of the suffixes.
// They store this information for all values of k from 1 to, say, 12.
// If the maximum k is B (e.g. 12), and alphabet size is A (e.g. 4 for
// DNA), then the buckets comprise A^B + A^(B-1) + ... + A + 1
// numbers.

// In addition, we allow skipping of positions while sorting and
// matching, using a periodic spaced seed.

#ifndef SUFFIXARRAY_HH
#define SUFFIXARRAY_HH
#include <string>
#include <vector>

namespace cbrc{

class PeriodicSpacedSeed;

class SuffixArray{
public:
  typedef unsigned indexT;
  typedef unsigned char uchar;

  // Add (unsorted) indices for text positions between beg and end,
  // and return 1.  Positions starting with delimiters aren't added.
  // If the index size would exceed maxBytes, don't add anything, and
  // return 0.
  int addIndices( const uchar* text, indexT beg, indexT end, indexT step,
		  unsigned alphSize, std::size_t maxBytes );

  indexT indexSize() const{ return index.size(); }
  std::size_t indexBytes() const{ return index.size() * sizeof(indexT); }

  // Sort the suffix array (but don't make the buckets).
  void sortIndex( const uchar* text, const PeriodicSpacedSeed& seed,
		  unsigned alphSize );

  // Make the buckets.  If bucketDepth == -1u, then a default
  // bucketDepth is used.  The default is: the maximum possible
  // bucketDepth such that the number of bucket entries is at most 1/4
  // the number of suffix array entries.
  void makeBuckets( const uchar* text, const PeriodicSpacedSeed& seed,
		    unsigned alphSize, indexT bucketDepth );

  // Return the maximum prefix size covered by the buckets.
  indexT maxBucketPrefix() const { return bucketSteps.size() - 1; }

  void clear();

  void fromFiles( const std::string& baseName, indexT indexNum,
		  indexT bucketDepth, const PeriodicSpacedSeed& seed,
		  unsigned alphSize );

  void toFiles( const std::string& baseName ) const;

  // Find the smallest match to the text, starting at the given
  // position in the query, such that there are at most maxHits
  // matches, and the match-depth is at least minDepth.  Return the
  // range of matching indices via beg and end.
  void match( const indexT*& beg, const indexT*& end,
	      const uchar* queryPtr, const uchar* text,
	      const PeriodicSpacedSeed& seed, unsigned alphSize,
	      indexT maxHits, indexT minDepth ) const;

  // Count matches of all sizes, starting at the given position in the
  // query.  Don't try this for large self-comparisons!
  void countMatches( std::vector<unsigned long long>& counts,
		     const uchar* queryPtr, const uchar* text,
		     const PeriodicSpacedSeed& seed, unsigned alphSize ) const;

private:
  std::vector<indexT> index;  // sorted indices
  std::vector<indexT> buckets;
  std::vector<indexT> bucketSteps;  // step size for each k-mer
  std::vector<indexT> bucketMask;  // extended mask, to avoid modulus
  indexT bucketMaskTotal;

  static void equalRange( const indexT*& beg, const indexT*& end,
			  const uchar* textBase, uchar symbol );
  static const indexT* lowerBound( const indexT* beg, const indexT* end,
				   const uchar* textBase, uchar symbol );
  static const indexT* upperBound( const indexT* beg, const indexT* end,
				   const uchar* textBase, uchar symbol );

  indexT defaultBucketDepth( unsigned alphSize );
  void makeBucketSteps( unsigned alphSize, indexT bucketDepth );
  void makeBucketMask( const PeriodicSpacedSeed& seed, indexT bucketDepth );
};

}  // end namespace
#endif
