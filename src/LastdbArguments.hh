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
  std::string maskPattern;
  indexT indexStep;
  std::size_t volumeSize;  // type?
  std::string userAlphabet;
  indexT bucketDepth;
  int verbosity;

  // positional arguments:
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace cbrc
#endif  // LASTDBARGUMENTS_HH
