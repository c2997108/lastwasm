// Copyright 2008 Martin C. Frith

// Generally useful input/output functions, mostly for binary reading
// and writing of vectors.

#ifndef IO_H
#define IO_H

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <cassert>

namespace cbrc{

// open an input file, but if the name is "-", just return cin
std::istream& openIn( const std::string& fileName, std::ifstream& ifs );

// open an output file, but if the name is "-", just return cout
std::ostream& openOut( const std::string& fileName, std::ofstream& ofs );

template <typename T>  // T should be a vector-iterator or a pointer
void memoryFromStream( T beg, T end, std::istream& s ){
  assert( beg < end );
  s.read( (char*)&(*beg), (end - beg) * sizeof(*beg) );
}

template<typename T>
void vectorFromStream( std::vector<T>& v, std::istream& s ){
  memoryFromStream( v.begin(), v.end(), s );
}

template<typename T>  // T should be a vector-iterator or a pointer
void memoryToStream( T beg, T end, std::ostream& s ){
  assert( beg < end );
  s.write( (const char*)&(*beg), (end - beg) * sizeof(*beg) );
}

template<typename T>
void vectorToStream( const std::vector<T>& v, std::ostream& s ){
  memoryToStream( v.begin(), v.end(), s );
}

template<typename T>  // T should be a vector-iterator or a pointer
void memoryFromBinaryFile( T beg, T end, const std::string& fileName ){
  if( beg == end ) return;
  std::ifstream file( fileName.c_str(), std::ios::binary );
  memoryFromStream( beg, end, file );
  if( !file ) throw std::runtime_error("can't read file: " + fileName);
}

template<typename T>
void vectorFromBinaryFile( std::vector<T>& v,
			   const std::string& fileName ){
  memoryFromBinaryFile( v.begin(), v.end(), fileName );
}

template<typename T>  // T should be a vector-iterator or a pointer
void memoryToBinaryFile( T beg, T end, const std::string& fileName ){
  if( beg == end ) return;
  std::ofstream file( fileName.c_str(), std::ios::binary );
  memoryToStream( beg, end, file );
  if( !file ) throw std::runtime_error("can't write file: " + fileName);
}

template<typename T>
void vectorToBinaryFile( const std::vector<T>& v,
			 const std::string& fileName ){
  memoryToBinaryFile( v.begin(), v.end(), fileName );
}

}  // end namespace cbrc
#endif  // IO_H
