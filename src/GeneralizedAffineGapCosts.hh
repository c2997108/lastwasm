// Copyright 2008, 2009 Martin C. Frith

// This struct holds parameters for so-called generalized affine gap
// costs (for pair-wise sequence alignment).  In this scheme, a "gap"
// may consist of unaligned regions in both sequences.  If these
// unaligned regions have sizes j and k, where j <= k, the cost is:

// a + b*(k-j) + c*j

// If c >= a + 2b, it reduces to standard affine gaps.  For more
// information, see: SF Altschul 1998 Proteins 32(1):88-96.

#ifndef GENERALIZEDAFFINEGAPCOSTS_HH
#define GENERALIZEDAFFINEGAPCOSTS_HH

namespace cbrc{

struct GeneralizedAffineGapCosts{
  int exist;
  int extend;
  int extendPair;
  int first;
  int firstPair;

  void assign( int a, int b, int c )
  { exist = a; extend = b; extendPair = c; first = a + b; firstPair = a + c; }

  // Will standard affine gaps always suffice for maximal alignment scores?
  bool isAffine() const { return (firstPair >= 2 * first); }

  // Return the score of a gap with the given sizes in a pair of
  // sequences, considering that it might be either one "generalized"
  // gap or two neighbouring "affine" gaps.
  int cost( int gapSize1, int gapSize2 ) const;
};

}  // end namespace cbrc
#endif  // GENERALIZEDAFFINEGAPCOSTS_HH
