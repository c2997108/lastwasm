// Copyright 2008, 2009 Martin C. Frith

#include "LastdbArguments.hh"
#include "stringify.hh"
#include <unistd.h>  // getopt
#include <stdexcept>

static void badopt( char opt, const char* arg ){
  throw std::runtime_error( std::string("bad option value: -") +
                            opt + ' ' + arg );
}

namespace cbrc{

LastdbArguments::LastdbArguments() :
  isProtein(false),
  isCaseSensitive(false),
  maskPattern("1"),
  indexStep(1),
  volumeSize(0x50000000),  // 1.25 Gbytes
  userAlphabet(""),
  bucketDepth(-1u),  // means: use the default (adapts to the data)
  verbosity(0){}

void LastdbArguments::fromArgs( int argc, char** argv ){
  std::string usage = "\
usage: lastdb [options] output-name fasta-sequence-file(s)\n\
\n\
Main Options (default settings):\n\
-p: interpret the sequences as proteins\n\
-c: read the sequences case-sensitively\n\
-m: periodic spaced-seed pattern (" + maskPattern + ")\n\
\n\
Advanced Options (default settings):\n\
-w: index step (" + stringify(indexStep) + ")\n\
-s: volume size (" + stringify(volumeSize) + ")\n\
-a: user-defined alphabet\n\
-b: bucket depth\n\
-v: be verbose: write messages about what lastdb is doing\n\
";

  int c;
  while( (c = getopt(argc, argv, "pcm:w:s:a:b:v")) != -1 ) {
    switch(c){
    case 'p':
      isProtein = true;
      break;
    case 'c':
      isCaseSensitive = true;
      break;
    case 'm':
      maskPattern = optarg;
      break;
    case 'w':
      unstringify( indexStep, optarg );
      if( indexStep < 1 ) badopt( c, optarg );
      break;
    case 's':
      unstringify( volumeSize, optarg );
      break;
    case 'a':
      userAlphabet = optarg;
      break;
    case 'b':
      unstringify( bucketDepth, optarg );
      break;
    case 'v':
      ++verbosity;
      break;
    case '?':
      throw std::runtime_error("bad option");
    }
  }

  if( optind == argc ) throw std::runtime_error(usage);
  lastdbName = argv[optind++];
  inputStart = optind;
}

}  // end namespace cbrc
