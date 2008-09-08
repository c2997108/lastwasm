// Copyright 2008 Martin C. Frith

#include "io.hh"
#include <iostream>

namespace cbrc{

std::istream& openIn( const std::string& fileName, std::ifstream& ifs ){
  if( fileName == "-" ) return std::cin;
  ifs.open( fileName.c_str() );
  if( !ifs ) throw std::runtime_error("can't open file: " + fileName);
  return ifs;
}

std::ostream& openOut( const std::string& fileName, std::ofstream& ofs ){
  if( fileName == "-" ) return std::cout;
  ofs.open( fileName.c_str() );
  if( !ofs ) throw std::runtime_error("can't open file: " + fileName);
  return ofs;
}

}  // end namespace cbrc
