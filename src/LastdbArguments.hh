// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

// This struct holds the command line arguments for lastdb.

#ifndef LASTDBARGUMENTS_HH
#define LASTDBARGUMENTS_HH

#include "SequenceFormat.hh"

#include <string>

namespace cbrc{

struct LastdbArguments{
  typedef unsigned indexT;

  // set the parameters to their default values:
  LastdbArguments();

  // set parameters from a list of arguments:
  void fromArgs( int argc, char** argv );

  // options:
  bool isProtein;
  bool isCaseSensitive;
  std::string spacedSeed;
  std::size_t volumeSize;  // type?
  indexT indexStep;
  std::string subsetSeedFile;
  std::string userAlphabet;
  indexT bucketDepth;
  bool isCountsOnly;
  int verbosity;
  sequenceFormat::Enum inputFormat;

  // positional arguments:
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace
#endif
