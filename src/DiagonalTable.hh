// Copyright 2008 Martin C. Frith

// This struct records coverage of "diagonals" by gapless alignments,
// when comparing two sequences.  The diagonal is the coordinate in
// one sequence minus the coordinate in the other sequence.  We assume
// that one sequence is scanned sequentially, and the other is
// accessed randomly.  We record the furthest sequential position
// covered so far in each diagonal.  This lets us avoid triggering
// gapless alignments in places that are already covered.

// Since the number of diagonals may be huge, we map them to a smaller
// number of bins.  To keep the bin populations small, when checking
// if a position is covered, we discard information about earlier
// sequential positions.

// This stuff uses a non-negligible fraction of the running time, and
// I have a hunch it can be done more efficiently...

#ifndef DIAGONALTABLE_HH
#define DIAGONALTABLE_HH
#include <vector>

namespace cbrc{

struct DiagonalTable{
  typedef unsigned indexT;
  typedef std::pair<indexT, indexT> pairT;

  enum { BINS = 65536 };  // use a power-of-two for faster modulus (maybe)

  // is this position on this diagonal already covered by an alignment?
  bool isCovered( indexT sequentialPos, indexT randomPos );

  // add an alignment endpoint to the table:
  void addEndpoint( indexT sequentialPos, indexT randomPos );

  std::vector<pairT> hits[BINS];
};

}  // end namespace cbrc
#endif  // DIAGONALTABLE_HH
