// Copyright 2008, 2009 Martin C. Frith

// This struct holds a suffix array.  The suffix array is just a list
// of numbers indicating positions in a text, sorted according to the
// alphabetical order of the text suffixes starting at these
// positions.  A query sequence can then be matched incrementally to
// the suffix array using binary search.

// For faster matching, we use "buckets", which store the matching
// regions for all k-mers for k=1 to, say, 12.  If the maximum k is B
// (e.g. 12), and alphabet size is A (e.g. 4 for DNA), then the
// buckets comprise A^B + A^(B-1) + ... + A + 1 numbers.  These
// numbers give the start and end offsets in the suffix array of every
// k-mer.

// In addition, we allow skipping of positions while sorting and
// matching, using a periodic spaced seed.

#ifndef SUFFIXARRAY_HH
#define SUFFIXARRAY_HH
#include <string>
#include <vector>

namespace cbrc{

struct SuffixArray{
  typedef unsigned indexT;
  typedef unsigned char uchar;

  SuffixArray( const std::vector<uchar>& t, const std::vector<indexT>& m,
	       unsigned a )
    : text(t), mask(m), alphSize(a), maskSize( m.size() ){}

  // Add (unsorted) indices for text positions between beg and end,
  // and return 1.  Positions starting with delimiters aren't added.
  // If the index size would exceed maxBytes, don't add anything, and
  // return 0.
  int makeIndex( indexT beg, indexT end, indexT step, std::size_t maxBytes );

  std::size_t indexBytes() const{ return index.size() * sizeof(indexT); }
  void sortIndex();
  void makeBuckets( indexT bucketDepth );
  void clear();

  void fromFiles( const std::string& baseName,
		  indexT indexNum, indexT bucketDepth );
  void toFiles( const std::string& baseName ) const;

  // Find the smallest match to the text, starting at the given
  // position in the query, such that there are at most maxHits
  // matches, and the match-depth is at least minDepth.  Return the
  // range of matching indices via beg and end.
  void match( const indexT*& beg, const indexT*& end,
	      const uchar* queryPtr,
	      indexT maxHits, indexT minDepth ) const;

  // Count matches of all sizes, starting at the given position in the
  // query.  Don't try this for large self-comparisons!
  void countMatches( std::vector<unsigned long long>& counts,
		     const uchar* queryPtr ) const;

  const std::vector<uchar>& text;
  const std::vector<indexT>& mask;
  unsigned alphSize;  // excluding delimiter symbol
  unsigned maskSize;  // same as mask.size(), but faster!!!
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

  void makeBucketSteps( indexT bucketDepth );
  void makeBucketMask( indexT bucketDepth );

  struct Stack{ indexT* b; indexT* e; indexT d; const uchar* t; };

  void radixSort( indexT* beg, indexT* end,
		  indexT depth, const uchar* textBase );
  void radixSortAlph4( indexT* beg, indexT* end,
		       indexT depth, const uchar* textBase );
  void insertionSort( indexT* beg, indexT* end,
		      indexT depth, const uchar* textBase );
  void insertionSortSimple( indexT* beg, indexT* end,
			    const uchar* textBase );
};

}  // end namespace cbrc
#endif  // SUFFIXARRAY_HH
