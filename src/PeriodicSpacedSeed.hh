// Copyright 2008 Martin C. Frith

// This struct holds a periodic spaced seed, represented by a string
// of 1s and 0s, such as "110".  This means: skip every third position
// (when comparing two sequences, starting at some pair of positions).

#ifndef PERIODICSPACEDSEED_HH
#define PERIODICSPACEDSEED_HH
#include <vector>
#include <string>
#include <iosfwd>

namespace cbrc{

struct PeriodicSpacedSeed{
  void fromString( const std::string& s );
  void init();

  std::string pattern;            // e.g. "110"
  std::vector<unsigned> offsets;  // e.g. {1, 2}
  unsigned maxOffset;             // the maximum offset
};

std::ostream& operator<<( std::ostream& s, const PeriodicSpacedSeed& m );
std::istream& operator>>( std::istream& s, PeriodicSpacedSeed& m );

}  // end namespace cbrc
#endif  // PERIODICSPACEDSEED_HH
