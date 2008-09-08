// Copyright 2008 Martin C. Frith

#include "PeriodicSpacedSeed.hh"
#include <stdexcept>
#include <istream>
#include <ostream>

namespace cbrc{

void PeriodicSpacedSeed::fromString( const std::string& s ){
  pattern = s;
  init();
}

void PeriodicSpacedSeed::init(){
  if( pattern.empty() || pattern[0] != '1' ){
    throw std::runtime_error("bad mask: " + pattern);
  }

  offsets.clear();
  unsigned off = 1;

  for( std::size_t i = 1; i < pattern.size(); ++i ){
    if( pattern[i] == '1' ){
      offsets.push_back(off);
      off = 0;
    }
    ++off;
  }

  offsets.push_back(off);

  maxOffset = *max_element( offsets.begin(), offsets.end() );
}

std::ostream& operator<<( std::ostream& s, const PeriodicSpacedSeed& m ){
  return s << m.pattern;
}

std::istream& operator>>( std::istream& s, PeriodicSpacedSeed& m ){
  s >> m.pattern;
  if( s ) m.init();
  return s;
}

}  // end namespace cbrc
