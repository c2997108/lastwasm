// Copyright 2010 Martin C. Frith

#include "fileMap.hh"

#include "stringify.hh"

#include <cassert>
#include <stdexcept>

// File mapping requires non-standard-C++ library functions.  I think
// this code will not work on all platforms, e.g. windows.  Hopefully,
// it's easy to rewrite this code for those platforms.

#include <fcntl.h>  // open
#include <unistd.h>  // close
#include <sys/mman.h>  // mmap, munmap

// This function tries to force the file-mapping to actually get
// loaded into memory, by reading it sequentially.  Without this,
// random access can be horribly slow (at least on two Linux 2.6
// systems).
static void primeMemory( void* begin, std::size_t bytes ){
  // use "static" to stop the compiler optimizing the whole function away:
  static unsigned z = 0;
  std::size_t stepSize = 1024;
  const char* x = static_cast<char*>(begin);
  const char* y = x + (bytes / stepSize) * stepSize;
  while( x < y ){
    z += *x;
    x += stepSize;
  }
}

namespace cbrc{

void* openFileMap( const std::string& fileName, std::size_t bytes ){
  if( bytes == 0 ) return 0;

  int f = open( fileName.c_str(), O_RDONLY );
  void* m = mmap( 0, bytes, PROT_READ, MAP_SHARED, f, 0 );
  int e = close(f);

  if( m == MAP_FAILED || e < 0 )
    throw std::runtime_error( "can't map file: " + fileName );

  primeMemory( m, bytes );

  return m;
}

void closeFileMap( void* begin, std::size_t bytes ){
  if( bytes == 0 ) return;
  int e = munmap( begin, bytes );
  if( e < 0 ) throw std::runtime_error( "failed to \"munmap\" " +
                                        stringify(bytes) + " bytes");
}

}  // end namespace
