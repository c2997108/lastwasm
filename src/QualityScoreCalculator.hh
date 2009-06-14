// Copyright 2009 Martin C. Frith

// This class calculates a Position-Specific-Scoring-Matrix from
// sequence quality scores.  The "init" routine creates some lookup
// tables, so that "makePssm" is fast.

#ifndef QUALITY_SCORE_CALCULATOR
#define QUALITY_SCORE_CALCULATOR

#include <cstddef>  // size_t

namespace cbrc{

class QualityScoreCalculator{
  typedef unsigned char uchar;
  enum { MAT = 64 };

public:
  void init( const int scoreMatrix[MAT][MAT],
	     unsigned alphabetSize,
	     double temperature,
	     bool isCaseSensitiveMatrix,
	     bool isMatchMismatchMatrix,
	     int matchScore,
	     int mismatchScore,  // score not cost!
	     const uchar* toUppercaseCode,
	     bool isPhred,  // phred scores or solexa scores?
	     int asciiOffset );  // how are the quality scores encoded?

  void makePssm( int pssm[][MAT],
		 const uchar* qualityScores,
		 const uchar* sequence,
		 std::size_t seqSize,
		 bool isSingleQualities );  // one quality score per position?

private:
  const int (*scoreMatrix)[MAT];
  unsigned alphabetSize;
  double temperature;
  bool isCaseSensitiveMatrix;
  bool isMatchMismatchMatrix;
  int matchScore;
  int mismatchScore;
  const uchar* toUppercaseCode;

  double qualityToProbCorrect[256];
  double qualityToProbOther[256];
  int qualityToMatchScore[256];
  int qualityToMismatchScore[256];
  double likelihoodRatioMatrix[MAT][MAT];
  bool isUniformRow[MAT];
  uchar toLowercaseCode[256];

  int probCorrectToMatchScore( double probCorrect );
};

}

#endif
