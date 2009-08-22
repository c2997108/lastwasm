// Copyright 2009 Martin C. Frith

#include "QualityScoreCalculator.hh"

//#include <iostream>  // for debugging
#include <numeric>  // inner_product
#include <algorithm>  // fill_n
#include <cmath>
//#include <cassert>  // for debugging

using namespace cbrc;

static int nearestInt( double x ){
  if( x > 0 ) return int( 0.5 + x );
  else return -int( 0.5 - x );
}

static double phredScoreToProbCorrect( int score ){
  // XXX real data sometimes has score=0, and perhaps this really
  // means that the error probability is 3/4, not 1.
  if( score < 0 ) return 0;  // there shouldn't be any scores < 0
  return 1 - std::pow(10, -0.1 * score);
}

static double solexaScoreToProbCorrect( int score ){
  return 1 / (1 + std::pow(10, -0.1 * score));
}

int QualityScoreCalculator::probCorrectToMatchScore( double probCorrect ){
  if( probCorrect <= 0 ) return mismatchScore;  // avoid underflow & log(0)
  double matchLR = std::exp( matchScore / temperature );
  double mismatchLR = std::exp( mismatchScore / temperature );
  double averageLR = probCorrect * matchLR + (1-probCorrect) * mismatchLR;
  return nearestInt( temperature * std::log( averageLR ) );
}

void QualityScoreCalculator::init( const int scoreMatrix[MAT][MAT],
				   unsigned alphabetSize,
				   double temperature,
				   bool isCaseSensitiveMatrix,
				   bool isMatchMismatchMatrix,
				   int matchScore,
				   int mismatchScore,
				   const uchar* toUppercaseCode,
				   bool isPhred,
				   int asciiOffset ){
  this->scoreMatrix = scoreMatrix;
  this->alphabetSize = alphabetSize;
  this->temperature = temperature;
  this->isCaseSensitiveMatrix = isCaseSensitiveMatrix;
  this->isMatchMismatchMatrix = isMatchMismatchMatrix;
  this->matchScore = matchScore;
  this->mismatchScore = mismatchScore;
  this->toUppercaseCode = toUppercaseCode;

  for( unsigned i = 0; i < 256; ++i ){
    int qualityScore = int(i) - asciiOffset;
    double probCorrect = ( isPhred ? phredScoreToProbCorrect(qualityScore)
			   :         solexaScoreToProbCorrect(qualityScore) );
    double probOther = (1 - probCorrect) / (alphabetSize - 1);
    qualityToProbCorrect[i] = probCorrect;
    qualityToProbOther[i] = probOther;
    qualityToMatchScore[i] = probCorrectToMatchScore(probCorrect);
    qualityToMismatchScore[i] = probCorrectToMatchScore(probOther);
  }

  for( unsigned i = 0; i < MAT; ++i ){
    for( unsigned j = 0; j < MAT; ++j ){
      likelihoodRatioMatrix[i][j]
	= std::exp( scoreMatrix[i][j] / temperature );  // can underflow to 0
    }
  }

  for( unsigned i = 0; i < MAT; ++i ){
    isUniformRow[i] = true;
    for( unsigned j = 1; j < alphabetSize; ++j ){
      if( scoreMatrix[i][j] != scoreMatrix[i][0] ) isUniformRow[i] = false;
    }
  }

  for( unsigned i = 0; i < 256; ++i ){
    unsigned j = toUppercaseCode[i];
    if( isCaseSensitiveMatrix && j != i ) continue;
    toLowercaseCode[j] = i;
  }
}

void QualityScoreCalculator::makePssm( int pssm[][MAT],
				       const uchar* qualityScores,
				       const uchar* sequence,
				       std::size_t seqSize,
				       bool isSingleQualities ) const{
  for( std::size_t i = 0; i < seqSize; ++i ){
    unsigned letter = sequence[i];
    if( !isCaseSensitiveMatrix ) letter = toUppercaseCode[letter];

    if( letter >= alphabetSize ){
      for( unsigned j = 0; j < MAT; ++j ){
        pssm[i][j] = scoreMatrix[j][letter];
      }
      continue;
    }

    // This special case for match-mismatch matrices is unnecessary,
    // but faster.  Unfortunately, it can give slightly different
    // results for PRB data (if the PRB probabilities don't add up to
    // exactly 1).
    if( isMatchMismatchMatrix ){
      std::fill_n( pssm[i], unsigned(MAT), mismatchScore );
      if( isSingleQualities ){
	uchar q = qualityScores[i];
	for( unsigned j = 0; j < alphabetSize; ++j ){
	  pssm[i][j] = qualityToMismatchScore[q];
	  pssm[i][ toLowercaseCode[j] ] = qualityToMismatchScore[q];
	}
	pssm[i][letter] = qualityToMatchScore[q];
	pssm[i][ toLowercaseCode[letter] ] = qualityToMatchScore[q];
      }
      else{
	for( unsigned j = 0; j < alphabetSize; ++j ){
	  uchar q = qualityScores[ i * alphabetSize + j ];
	  pssm[i][j] = qualityToMatchScore[q];
	  pssm[i][ toLowercaseCode[j] ] = qualityToMatchScore[q];
	}
      }
      pssm[i][alphabetSize] = scoreMatrix[alphabetSize][0];
      continue;
    }

    double letterProbs[256];

    if( isSingleQualities ){
      uchar q = qualityScores[i];
      std::fill_n( letterProbs, alphabetSize, qualityToProbOther[q] );
      letterProbs[letter] = qualityToProbCorrect[q];
    }
    else{
      for( unsigned j = 0; j < alphabetSize; ++j ){
	uchar q = qualityScores[ i * alphabetSize + j ];
	letterProbs[j] = qualityToProbCorrect[q];
      }
    }

    for( unsigned j = 0; j < MAT; ++j ){
      if( isUniformRow[j] ){
	pssm[i][j] = scoreMatrix[j][0];
	continue;
      }
      // I think this is the right way round for a non-symmetric score matrix:
      double likelihoodRatio = std::inner_product( letterProbs,
						   letterProbs + alphabetSize,
						   likelihoodRatioMatrix[j],
						   0.0 );
      // XXX likelihoodRatio can be zero (if we have very negative
      // match/mismatch scores, and maybe a Phred score of zero).
      pssm[i][j] = nearestInt( temperature * std::log(likelihoodRatio) );
    }
  }
}
