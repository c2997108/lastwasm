// Copyright 2008 Martin C. Frith

// This struct holds the command line arguments for lastdb.

#ifndef LASTDBARGUMENTS_HH
#define LASTDBARGUMENTS_HH
#include <string>

namespace cbrc{

struct LastdbArguments{
  typedef unsigned indexT;

  LastdbArguments( int argc, char** argv );

  static void badopt( char opt, const char* arg );

  // options:
  bool isProtein;
  bool isCaseSensitive;
  std::string userAlphabet;
  std::string maskPattern;
  indexT indexStep;
  indexT bucketDepth;
  std::size_t volumeSize;  // type?
  int verbosity;

  // positional arguments:
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace cbrc
#endif  // LASTDBARGUMENTS_HH
