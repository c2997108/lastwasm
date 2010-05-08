// Copyright 2008, 2009, 2010 Martin C. Frith

#include "LastdbArguments.hh"
#include "stringify.hh"
#include <unistd.h>  // getopt
#include <iostream>
#include <stdexcept>
#include <cstdlib>  // EXIT_SUCCESS

#define ERR(x) throw std::runtime_error(x)

static void badopt( char opt, const char* arg ){
  ERR( std::string("bad option value: -") + opt + ' ' + arg );
}

using namespace cbrc;

LastdbArguments::LastdbArguments() :
  isProtein(false),
  isCaseSensitive(false),
  spacedSeed(""),
  volumeSize(0x50000000),  // 1.25 Gbytes
  indexStep(1),
  subsetSeedFile(""),
  userAlphabet(""),
  bucketDepth(indexT(-1)),  // means: use the default (adapts to the data)
  verbosity(0){}

void LastdbArguments::fromArgs( int argc, char** argv ){
  std::string usage = "\
Usage: lastdb [options] output-name fasta-sequence-file(s)\n\
Prepare sequences for subsequent alignment with lastal.\n\
\n\
Main Options:\n\
-h: show all options and their default settings\n\
-p: interpret the sequences as proteins\n\
-c: soft-mask lowercase letters";

    std::string help = usage + "\n\
\n\
Advanced Options (default settings):\n\
-s: volume size (1280 MiB)\n\
-m: spaced seed pattern\n\
-u: subset seed file (yass.seed)\n\
-w: index step (" + stringify(indexStep) + ")\n\
-a: user-defined alphabet\n\
-b: bucket depth\n\
-v: be verbose: write messages about what lastdb is doing\n\
\n\
Report bugs to: last (ATmark) cbrc (dot) jp\n\
LAST home page: http://last.cbrc.jp/\n\
";

  int c;
  while( (c = getopt(argc, argv, "hpcm:s:w:u:a:b:v")) != -1 ) {
    switch(c){
    case 'h':
      std::cout << help;
      throw EXIT_SUCCESS;
    case 'p':
      isProtein = true;
      break;
    case 'c':
      isCaseSensitive = true;
      break;
    case 'm':
      spacedSeed = optarg;
      break;
    case 's':
      unstringifySize( volumeSize, optarg );
      break;
    case 'w':
      unstringify( indexStep, optarg );
      if( indexStep < 1 ) badopt( c, optarg );
      break;
    case 'u':
      subsetSeedFile = optarg;
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
      ERR( "bad option" );
    }
  }

  if( optind + 1 >= argc )
    ERR( "please give me an output name and sequence file(s)\n\n" + usage );
  lastdbName = argv[optind++];
  inputStart = optind;
}
