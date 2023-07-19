// Author Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

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

#ifndef DIAGONALTABLE_HH
#define DIAGONALTABLE_HH

#include <stddef.h>
#include <utility>  // pair
#include <vector>

namespace cbrc {

struct DiagonalTable {

  typedef std::pair<size_t, size_t> pairT;

  enum { BINS = 256 };  // use a power-of-two for faster modulus (maybe)
                        // 256 is much faster than 65536 in my tests

  // is this position on this diagonal already covered by an alignment?
  bool isCovered(size_t diagonal, size_t sequentialPos) {
    std::vector<pairT> &v = hits[diagonal % BINS];

    for (std::vector<pairT>::iterator i = v.begin(); i < v.end(); ) {
      if (i->first >= sequentialPos) {
	if (i->second == diagonal) return true;
	++i;
      } else {
	i = v.erase(i);  // hopefully we rarely get here
      }
    }

    return false;
  }

  // add an alignment endpoint to the table:
  void addEndpoint(size_t diagonal, size_t sequentialPos) {
    hits[diagonal % BINS].push_back(pairT(sequentialPos, diagonal));
  }

  std::vector<pairT> hits[BINS];
};

}

#endif
