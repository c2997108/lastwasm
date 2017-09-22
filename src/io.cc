// Copyright 2008, 2009, 2010 Martin C. Frith

#include "io.hh"
#include <iostream>
#include <iterator>  // istreambuf_iterator

namespace cbrc{

std::istream& openIn( const std::string& fileName, mcf::izstream& ifs ){
  if( fileName == "-" ) return std::cin;
  ifs.open( fileName.c_str() );
  if( !ifs ) throw std::runtime_error("can't open file: " + fileName);
  return ifs;
}

std::string slurp( const std::string& fileName ){
  mcf::izstream inFileStream;
  std::istream& in = openIn( fileName, inFileStream );
  std::istreambuf_iterator<char> beg(in), end; 
  return std::string( beg, end );
}

}  // end namespace cbrc
