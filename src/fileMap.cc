// Copyright 2010, 2013 Martin C. Frith

#include "fileMap.hh"

#include "stringify.hh"

#include <cstring>  // strerror
#include <cerrno>
#include <stdexcept>

// File mapping requires non-standard-C++ library functions.  I think
// this code will not work on all platforms, e.g. windows.  Hopefully,
// it's easy to rewrite this code for those platforms.

#ifndef __EMSCRIPTEN__
#include <fcntl.h>  // open
#include <unistd.h>  // close
#include <sys/mman.h>  // mmap, munmap
#else
#include <cstdio>
#include <cstdlib>
#include <iostream>
#endif

static void err( const std::string& s ) {
  throw std::runtime_error( s + ": " + std::strerror(errno) );
}

// This function tries to force the file-mapping to actually get
// loaded into memory, by reading it sequentially.  Without this,
// random access can be horribly slow (at least on two Linux 2.6
// systems).
static void primeMemory( void* begin, size_t bytes ){
  unsigned z = 0;
  size_t stepSize = 1024;
  const char* x = static_cast<char*>(begin);
  const char* y = x + (bytes / stepSize) * stepSize;
  while( x < y ){
    z += *x;
    x += stepSize;
  }
  volatile unsigned dontOptimizeMeAway = z;
  dontOptimizeMeAway = dontOptimizeMeAway;  // ??? prevents compiler warning
}

namespace cbrc{

void* openFileMap( const std::string& fileName, size_t bytes ){
  if( bytes == 0 ) return 0;
#ifndef __EMSCRIPTEN__
  int f = open( fileName.c_str(), O_RDONLY );
  if( f < 0 ) err( "can't open file " + fileName );

  void* m = mmap( 0, bytes, PROT_READ, MAP_SHARED, f, 0 );
  if( m == MAP_FAILED ) err( "can't map file " + fileName );

  int e = close(f);
  if( e < 0 ) err( "can't close file " + fileName );

  primeMemory( m, bytes );

  return m;
#else
  std::cerr << "lastal.js: fileMap open '" << fileName << "' bytes=" << bytes << "\n";
  FILE* f = std::fopen(fileName.c_str(), "rb");
  if (!f) err("can't open file " + fileName);
  void* m = std::malloc(bytes);
  if (!m) { std::fclose(f); throw std::runtime_error("out of memory mapping " + fileName); }
  size_t n = std::fread(m, 1, bytes, f);
  std::fclose(f);
  if (n != bytes) {
    std::free(m);
    throw std::runtime_error(std::string("short read: ") + fileName +
                             " got=" + stringify(n) + " expected=" + stringify(bytes));
  }
  std::cerr << "lastal.js: fileMap mapped '" << fileName << "' ok\n";
  return m;
#endif
}

void closeFileMap( void* begin, size_t bytes ){
  if( bytes == 0 ) return;
#ifndef __EMSCRIPTEN__
  int e = munmap( begin, bytes );
  if( e < 0 ) err( "failed to \"munmap\" " + stringify(bytes) + " bytes" );
#else
  std::free(begin);
#endif
}

}  // end namespace
