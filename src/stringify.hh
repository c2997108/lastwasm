// Copyright 2008 Martin C. Frith

// Convert things (e.g. numbers) to and from strings.  These are
// replacements for boost::lexical_cast: by avoiding dependency on
// boost, we can distribute code more easily.

#ifndef STRINGIFY_HH
#define STRINGIFY_HH
#include <sstream>
#include <string>
#include <stdexcept>
#include <cassert>

namespace cbrc{

template<typename T>
std::string stringify( const T& x ){
  std::ostringstream oss;
  oss << x;
  assert(oss);
  return oss.str();
}

template<typename T>
void unstringify( T& x, const std::string& s ){
  std::istringstream iss(s);
  if( !(iss >> x) || !(iss >> std::ws).eof() ){
    throw std::runtime_error( "can't interpret: " + s );
  }
}

}  // end namespace cbrc
#endif  // STRINGIFY_HH
