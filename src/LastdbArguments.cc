// Copyright 2008 Martin C. Frith

#include "LastdbArguments.hh"
#include "stringify.hh"
#include <stdexcept>

namespace cbrc{

void LastdbArguments::badopt( char opt, const char* arg ){
  throw std::runtime_error( std::string("bad option value: -") +
                            opt + ' ' + arg );
}

LastdbArguments::LastdbArguments( int argc, char** argv ) :
  // default values:
  isProtein(false),
  isCaseSensitive(false),
  userAlphabet(""),
  maskPattern("1"),
  indexStep(1),
  bucketDepth(-1u),  // means: use the default (adapts to the data)
  volumeSize(0x50000000),  // 1.25 Gbytes
  verbosity(0)
{
  std::string usage = "\
usage: lastdb [options] output-name fasta-sequence-file(s)\n\
\n\
Options (default settings):\n\
-p: interpret the sequences as proteins\n\
-c: read the sequences case-sensitively\n\
-a: user-defined alphabet\n\
-m: periodic spaced-seed pattern (" + maskPattern + ")\n\
-w: index step (" + stringify(indexStep) + ")\n\
-b: bucket depth\n\
-s: volume size (" + stringify(volumeSize) + ")\n\
-v: be verbose: write messages about what lastdb is doing\n\
";

  int c;
  while( (c = getopt(argc, argv, "pca:m:w:b:s:v")) != -1 ) {
    switch(c){
    case 'p':
      isProtein = true;
      break;
    case 'c':
      isCaseSensitive = true;
      break;
    case 'a':
      userAlphabet = optarg;
      break;
    case 'm':
      maskPattern = optarg;
      break;
    case 'w':
      unstringify( indexStep, optarg );
      if( indexStep < 1 ) badopt( c, optarg );
      break;
    case 'b':
      unstringify( bucketDepth, optarg );
      break;
    case 's':
      unstringify( volumeSize, optarg );
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
