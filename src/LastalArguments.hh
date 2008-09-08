// Copyright 2008 Martin C. Frith

// This struct holds the command line arguments for lastal.

#ifndef LASTALARGUMENTS_HH
#define LASTALARGUMENTS_HH
#include <string>

namespace cbrc{

struct LastalArguments{
  typedef unsigned indexT;

  LastalArguments( int argc, char** argv );

  // set default option values that depend on input files:
  void setDefaults( bool isDna, bool isProtein, int maxMatchScore );

  void writeCommented( std::ostream& stream ) const;

  static void badopt( char opt, const char* arg );

  // options:
  std::string outFile;
  int outputFormat;
  int outputType;
  int strand;
  int maskLowercase;
  int minScoreGapped;
  int minScoreGapless;
  int matchScore;
  int mismatchCost;
  int gapExistCost;
  int gapExtendCost;
  int gapPairCost;
  std::string matrixFile;
  int maxDropGapped;
  int maxDropGapless;
  indexT minHitDepth;
  indexT oneHitMultiplicity;
  indexT queryStep;
  indexT batchSize;  // approx size of query sequences to scan in 1 batch
  indexT maxRepeatDistance;  // supress repeats <= this distance apart
  int verbosity;

  // positional arguments:
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace cbrc
#endif  // LASTALARGUMENTS_HH
