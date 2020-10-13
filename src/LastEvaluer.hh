// Copyright 2015 Martin C. Frith

// This class calculates E-values for pair-wise local alignment.
// It has 2 states: "good" and "bad".  It starts in the "bad" state.

// "database" = sequence1, "query" = sequence2

// For DNA-versus-protein alignment:
// protein = "database" = sequence1, DNA = "query" = sequence2

// "deletion" means deletion in sequence 2 relative to sequence 1
// "insertion" means insertion in sequence 2 relative to sequence 1

// A length-k deletion costs delOpen + k * delEpen
// A length-k insertion costs insOpen + k * insEpen

#ifndef LAST_EVALUER_HH
#define LAST_EVALUER_HH

#include "ScoreMatrixRow.hh"

#include "alp/sls_alignment_evaluer.hpp"
#include "alp/sls_falp_alignment_evaluer.hpp"

namespace cbrc {

class GeneticCode;

class LastEvaluer {
public:
  // This routine tries to initialize the object for a given set of
  // alignment parameters.  It may fail, i.e. set the object to the
  // "bad" state and throw an Sls::error.
  // These arguments are only used to lookup pre-calculated cases:
  // matrixName, matchScore, mismatchCost, geneticCodeName.
  // DNA-versus-protein alignment is indicated by: frameshiftCost >= 0.
  // As a special case, frameshiftCost==0 means no frameshifts.
  // For DNA-versus-protein alignment, letterFreqs2 should either be
  // NULL or point to 64 codon frequencies (aaa, aac, etc).
  void init(const char *matrixName,
	    int matchScore,
	    int mismatchCost,
	    const char *alphabet,
	    const ScoreMatrixRow *scoreMatrix,  // score[sequence1][sequence2]
	    const double *letterFreqs1,
	    const double *letterFreqs2,
	    bool isGapped,
	    int delOpen,
	    int delEpen,
	    int insOpen,
	    int insEpen,
	    int frameshiftCost,
	    const GeneticCode &geneticCode,
	    const char *geneticCodeName,
	    int verbosity);

  void setSearchSpace(double databaseLength,  // number of database letters
		      double databaseMaxSeqLength,  // length of longest seq
		      double numOfStrands) {  // 1 or 2
    if (databaseMaxSeqLength > 0) {
      databaseSeqLen = databaseMaxSeqLength;
      databaseSeqNum = databaseLength / databaseMaxSeqLength * numOfStrands;
    } else {
      this->databaseSeqLen = 1;  // ALP doesn't like 0
      this->databaseSeqNum = 0;
    }
  }

  bool isGood() const
  { return evaluer.isGood() || frameshiftEvaluer.isGood(); }

  // Don't call this in the "bad" state:
  double evaluePerArea(double score) const {
    return evaluer.isGood() ?
      evaluer.evaluePerArea(score) : frameshiftEvaluer.evaluePerArea(score);
  }

  // Don't call this in the "bad" state or before calling setSearchSpace:
  double area(double score, double queryLength) const {
    return databaseSeqNum *
      (evaluer.isGood()
       ?           evaluer.area(score, queryLength, databaseSeqLen)
       : frameshiftEvaluer.area(score, queryLength, databaseSeqLen));
  }

  // Don't call this in the "bad" state:
  double bitScore(double score) const {
    return evaluer.isGood() ?
      evaluer.bitScore(score) : frameshiftEvaluer.bitScore(score);
  }

  // In the "good" state: returns the minimum score with E-value <=
  // "evalue", which is always > 0.  Otherwise: returns -1.
  int minScore(double evalue, double area) const;

  // In the "good" state, after calling setSearchSpace: returns the
  // minimum positive score with E-value <= 1 per this many query
  // letters.  In the "bad" state: returns -1.
  int minScore(double queryLettersPerRandomAlignment) const;

  // Writes some parameters preceded by "#".  Does nothing in the "bad" state.
  void writeCommented(std::ostream& out) const;

  // Writes all parameters in full precision.  Does nothing in the "bad" state.
  void writeParameters(std::ostream& out) const;

private:
  Sls::AlignmentEvaluer evaluer;
  Sls::FrameshiftAlignmentEvaluer frameshiftEvaluer;
  double databaseSeqLen;
  double databaseSeqNum;
};

}

#endif
