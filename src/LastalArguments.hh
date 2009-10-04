// Copyright 2008, 2009 Martin C. Frith

// This struct holds the command line arguments for lastal.

#ifndef LASTALARGUMENTS_HH
#define LASTALARGUMENTS_HH
#include <string>
#include <iosfwd>

namespace cbrc{

struct LastalArguments{
  typedef unsigned indexT;

  // set the parameters to their default values:
  LastalArguments();

  // set parameters from a list of arguments:
  void fromArgs( int argc, char** argv, bool optionsOnly = false );

  // set parameters from a command line (by splitting it into arguments):
  void fromLine( const std::string& line, bool optionsOnly = true );

  // set parameters from lines beginning with "#lastal":
  void fromStream( std::istream& is, bool optionsOnly = true );

  // set parameters from lines beginning with "#lastal":
  void fromString( const std::string& s, bool optionsOnly = true );

  // set default option values that depend on input files:
  void setDefaultsFromAlphabet( bool isDna, bool isProtein );
  void setDefaultsFromMatrix( double lambda );

  // write the parameter settings, starting each line with "#":
  void writeCommented( std::ostream& stream ) const;

  // are we doing translated alignment (DNA versus protein)?
  bool isTranslated() const{ return frameshiftCost > 0; }

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
  int frameshiftCost;
  std::string matrixFile;
  int maxDropGapped;
  int maxDropGapless;
  int inputFormat;
  indexT minHitDepth;
  indexT oneHitMultiplicity;
  indexT queryStep;
  indexT batchSize;  // approx size of query sequences to scan in 1 batch
  indexT maxRepeatDistance;  // supress repeats <= this distance apart
  double temperature;  // probability = exp( score / temperature ) / Z
  double gamma;        // parameter for gamma-centroid alignment
  std::string geneticCodeFile;
  int verbosity;

  // positional arguments:
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace cbrc
#endif  // LASTALARGUMENTS_HH
