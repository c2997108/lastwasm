// Copyright 2015 Martin C. Frith

#include "LastEvaluer.hh"

#include "GeneticCode.hh"
#include "mcf_alignment_path_adder.hh"

#include "alp/sls_falp_alignment_evaluer.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>

#include <random>

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

namespace cbrc {

#include "LastEvaluerData.hh"

static bool isEqual(const char *x, const char *y) {
  return std::strcmp(x, y) == 0;
}

static bool isHit(const EvalueParametersByName &p,
		  const char *n, int a, int b) {
  return isEqual(p.matrixName, n) && p.gapOpen == a && p.gapEpen == b;
}

static bool isHit(const EvalueParametersByScore &p,
		  int r, int q, int a, int b) {
  return p.matchScore == r && p.mismatchCost == q &&
    p.gapOpen == a && p.gapEpen == b;
}

static bool isHit(const FrameshiftEvalueParameters &p,
		  const char *g, const char *n, int a, int b, int f) {
  return isEqual(p.geneticCodeName, g) && isEqual(p.matrixName, n) &&
    p.gapOpen == a && p.gapEpen == b && p.frameshiftCost == f;
}

static bool isProtein(const char *alphabet) {
  return isEqual(alphabet, "ACDEFGHIKLMNPQRSTVWY");
}

static bool isDna(const char *alphabet) {
  return isEqual(alphabet, "ACGT");
}

static void makeMatrix(size_t size, long *data, long **rows) {
  for (size_t i = 0; i < size; ++i)
    rows[i] = data + i * size;
}

// We have to transpose the matrix, because: for DNA-versus-protein
// alignment the input matrix is score[protein][DNA], but the ALP
// library wants score[DNA][protein].

static void copyMatrix(size_t size, const ScoreMatrixRow *in, long **out) {
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j)
      out[i][j] = in[j][i];  // transpose
}

static void copyMatrix(size_t aaNum, size_t txNum, const ScoreMatrixRow *in,
		       long **out, const size_t *usedLetters) {
  for (size_t i = 0; i < txNum; ++i)
    for (size_t j = 0; j < aaNum; ++j)
      out[i][j] = in[j][usedLetters[i]];  // transpose
}

static void makeCodonTable(long *codonTable, const GeneticCode &geneticCode) {
  unsigned char c[3];
  for (c[0] = 0; c[0] < 4; ++c[0])
    for (c[1] = 0; c[1] < 4; ++c[1])
      for (c[2] = 0; c[2] < 4; ++c[2])
	*codonTable++ = geneticCode.translation(c);
}

static size_t findUsedLetters(size_t *usedLetters, size_t alphabetSize,
			      const long *codonTable) {
  const long *end = codonTable + 64;
  size_t j = 0;
  for (size_t i = 0; i < scoreMatrixRowSize; ++i)
    if (i < alphabetSize || std::find(codonTable, end, i) < end)
      usedLetters[j++] = i;
  return j;  // returns the number of used letters
}

static void shrinkCodonTable(long *codonTable, const size_t *usedLettersBeg,
			     const size_t *usedLettersEnd) {
  for (size_t i = 0; i < 64; ++i)
    codonTable[i] = std::find(usedLettersBeg,
			      usedLettersEnd, codonTable[i]) - usedLettersBeg;
}

static void makeNtFreqs(double *ntFreqs, const double *txFreqs) {
  ntFreqs[0] = ntFreqs[1] = ntFreqs[2] = ntFreqs[3] = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
	ntFreqs[i] += *txFreqs;
	ntFreqs[j] += *txFreqs;
	ntFreqs[k] += *txFreqs;
	++txFreqs;
      }
    }
  }
}

static void makeTxFreqs(double *txFreqs, const double *ntFreqs,
			const long *codonTable) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
	txFreqs[*codonTable++] += ntFreqs[i] * ntFreqs[j] * ntFreqs[k];
}

void LastEvaluer::init(const char *matrixName,
		       int matchScore,
		       int mismatchCost,
		       const char *alphabet,
		       const ScoreMatrixRow *scoreMatrix,
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
		       int verbosity) {
  const double lambdaTolerance = 0.01;
  const double kTolerance = 0.05;
  const double maxMegabytes = 500;
  const long randomSeed = 1;

  const double maxSeconds = 60.0;  // seems to work nicely on my computer
  size_t alphabetSize = std::strlen(alphabet);

  if (frameshiftCost >= 0) {  // DNA-versus-protein alignment:
    if (isGapped && insOpen == delOpen && insEpen == delEpen) {
      if (isProtein(alphabet)) {
	for (size_t i = 0; i < COUNTOF(frameshiftEvalueParameters); ++i) {
	  const FrameshiftEvalueParameters &p = frameshiftEvalueParameters[i];
	  if (isHit(p, geneticCodeName, matrixName, delOpen, delEpen,
		    frameshiftCost))
	    return evaluer.initParameters(p.parameters);
	}
      }
    }

    long codonTable[64];
    makeCodonTable(codonTable, geneticCode);
    size_t usedLetters[scoreMatrixRowSize];
    size_t matrixSize = findUsedLetters(usedLetters, alphabetSize, codonTable);
    shrinkCodonTable(codonTable, usedLetters, usedLetters + matrixSize);

    std::vector<long> matrixData(matrixSize * matrixSize);
    std::vector<long*> matrix(matrixSize);
    makeMatrix(matrixSize, &matrixData[0], &matrix[0]);
    copyMatrix(alphabetSize, matrixSize, scoreMatrix, &matrix[0], usedLetters);

    std::vector<double> ntFreqs(4, 0.25);
    if (letterFreqs2) makeNtFreqs(&ntFreqs[0], letterFreqs2);
    std::vector<double> aaFreqs(matrixSize);
    copy(letterFreqs1, letterFreqs1 + alphabetSize, aaFreqs.begin());

    if (isGapped && frameshiftCost > 0) {  // with frameshifts:
      Sls::FrameshiftAlignmentEvaluer frameshiftEvaluer;
      frameshiftEvaluer.initFrameshift(4, matrixSize, codonTable,
				       &matrix[0], &ntFreqs[0], &aaFreqs[0],
				       delOpen, delEpen, insOpen, insEpen,
				       frameshiftCost, true,
				       lambdaTolerance, kTolerance,
				       0, maxMegabytes, randomSeed);
      const Sls::FALP_set_of_parameters &p = frameshiftEvaluer.parameters();
      Sls::AlignmentEvaluerParameters q = {p.lambda, p.K,
					   p.a_I, p.b_I,  // !!! no flip
					   p.a_J, p.b_J,
					   p.alpha_I, p.beta_I,
					   p.alpha_J, p.beta_J,
					   p.sigma, p.tau};
      evaluer.initParameters(q);
    } else {  // without frameshifts:
      std::vector<double> txFreqs(matrixSize);
      makeTxFreqs(&txFreqs[0], &ntFreqs[0], codonTable);
      if (isGapped) {
	evaluer.set_gapped_computation_parameters_simplified(maxSeconds);
	evaluer.initGapped(matrixSize, &matrix[0], &txFreqs[0], &aaFreqs[0],
			   delOpen, delEpen, insOpen, insEpen,
			   true, lambdaTolerance, kTolerance,
			   0, maxMegabytes, randomSeed);
      } else {
	evaluer.initGapless(matrixSize, &matrix[0], &txFreqs[0], &aaFreqs[0]);
      }
      const Sls::ALP_set_of_parameters &p = evaluer.parameters();
      Sls::AlignmentEvaluerParameters q = {p.lambda, p.K,
					   p.a_J * 3, p.b_J * 3,
					   p.a_I, p.b_I,  // !!! flip (I, J)
					   p.alpha_J * 9, p.beta_J * 9,
					   p.alpha_I, p.beta_I,
					   p.sigma * 3, p.tau * 3};
      evaluer.initParameters(q);
    }
  } else {  // ordinary alignment:
    if (isGapped && insOpen == delOpen && insEpen == delEpen) {
      if (isProtein(alphabet)) {
	for (size_t i = 0; i < COUNTOF(proteinParameters); ++i) {
	  const EvalueParametersByName &p = proteinParameters[i];
	  if (isHit(p, matrixName, delOpen, delEpen))
	    return evaluer.initParameters(p.parameters);
	}
      }
      if (isDna(alphabet)) {
	if (*matrixName) {
	  for (size_t i = 0; i < COUNTOF(dnaParametersByName); ++i) {
	    const EvalueParametersByName &p = dnaParametersByName[i];
	    if (isHit(p, matrixName, delOpen, delEpen))
	      return evaluer.initParameters(p.parameters);
	  }
	} else {
	  for (size_t i = 0; i < COUNTOF(dnaParametersByScore); ++i) {
	    const EvalueParametersByScore &p = dnaParametersByScore[i];
	    if (isHit(p, matchScore, mismatchCost, delOpen, delEpen))
	      return evaluer.initParameters(p.parameters);
	  }
	}
      }
    }

    std::vector<long> matrixData(alphabetSize * alphabetSize);
    std::vector<long*> matrix(alphabetSize);
    makeMatrix(alphabetSize, &matrixData[0], &matrix[0]);
    copyMatrix(alphabetSize, scoreMatrix, &matrix[0]);

    if (isGapped) {
      evaluer.set_gapped_computation_parameters_simplified(maxSeconds);
      for (int i = 0; ; ++i) {
	double t = Sls::default_importance_sampling_temperature + 0.01 * i;
	if (verbosity > 0) std::cerr << "try temperature=" << t << " ";
	try {
	  evaluer.initGapped(alphabetSize, &matrix[0],
			     letterFreqs2, letterFreqs1,
			     delOpen, delEpen, insOpen, insEpen,
			     true, lambdaTolerance, kTolerance,
			     0, maxMegabytes, randomSeed, t);
	  if (verbosity > 0) std::cerr << "OK\n";
	  break;
	} catch (const Sls::error& e) {
	  if (verbosity > 0) std::cerr << "NG\n";
	  if (verbosity > 1) {
	    std::cerr << "ALP: " << e.error_code << ": " << e.st;
	  }
	  if (i == 20) throw;
	}
      }
    } else {
      evaluer.initGapless(alphabetSize, &matrix[0],
			  letterFreqs2, letterFreqs1);
    }
  }
}

void LastEvaluer::initFullScores(const const_dbl_ptr *substitutionProbs,
				 const double *letterFreqs1, int alphabetSize1,
				 const double *letterFreqs2, int alphabetSize2,
				 const GapCosts &gapCosts, double scale,
				 int numOfAlignments, int seqLength,
				 int verbosity, bool isFrameshift) {
  int border = isFrameshift ? 0 : seqLength / 4;
  int seqLength1 = seqLength + border;
  int seqLength2 = seqLength1;
  int seqLength3 = seqLength2;
  if (isFrameshift) {
    seqLength2 = seqLength1 * 3;
    seqLength3 = seqLength2 + 2;
  }

  const GapCosts::ProbPiece &del = gapCosts.delProbPieces[0];
  const GapCosts::ProbPiece &ins = gapCosts.insProbPieces[0];
  AlignmentPathAdder scorer;
  FrameshiftXdropAligner frameshiftScorer;

  std::mt19937_64 randGen;
  std::discrete_distribution<> dist1(letterFreqs1,
				     letterFreqs1 + alphabetSize1);
  std::discrete_distribution<> dist2(letterFreqs2,
				     letterFreqs2 + alphabetSize2);

  std::vector<uchar> seqs(seqLength1 + seqLength2);
  uchar *seq1 = &seqs[0];
  uchar *seq2 = seq1 + seqLength1;

  double probRatioSum = 0;

  for (int i = 0; i < numOfAlignments; ++i) {
    for (int j = 0; j < seqLength1; ++j) seq1[j] = dist1(randGen);
    for (int j = 0; j < seqLength2; ++j) seq2[j] = dist2(randGen);
    double p = isFrameshift
      ? frameshiftScorer.maxSumOfProbRatios(seq1, seqLength1, seq2, seqLength3,
					    substitutionProbs, gapCosts)
      : scorer.maxSum(seq1, seqLength1, seq2, seqLength3, substitutionProbs,
		      del.openProb, del.growProb, ins.openProb, ins.growProb,
		      border);
    probRatioSum += 1 / p;
    if (verbosity > 1) std::cerr << "simScore: " << (log(p) / scale) << "\n";
  }

  // max likelihood k  =  1 / (m * n * avg[exp(-lambda * score)])
  int area = (seqLength1 - border) * (seqLength3 - border);
  double k = numOfAlignments / (area * probRatioSum);

  if (verbosity > 1) {
    std::cerr << "lambda k m n: " << scale << " " << k << " "
	      << (seqLength1 - border) << " " << (seqLength3 - border) << "\n";
  }

  Sls::AlignmentEvaluerParameters p = {scale, k, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  evaluer.initParameters(p);
}

double LastEvaluer::minScore(double evalue, double area) const {
  const Sls::ALP_set_of_parameters &p = evaluer.parameters();
  // do log of evalue separately, to reduce the risk of overflow:
  double s = (log(p.K * area) - log(evalue)) / p.lambda;
  return std::max(0.0, s);
}

double LastEvaluer::minScore(double queryLettersPerRandomAlignment) const {
  double huge = 1e9;
  const Sls::ALP_set_of_parameters &p = evaluer.parameters();
  double x = queryLettersPerRandomAlignment * databaseSeqNum;
  double beg = 0;
  double len = log(x * databaseSeqLen * p.K) / p.lambda;

  while (1) {
    len /= 2;
    double mid = beg + len;
    if (mid <= beg) return mid;
    if (evaluer.evalue(mid, huge, databaseSeqLen) >= huge / x)
      beg = mid;
  }
}

void LastEvaluer::writeCommented(std::ostream& out) const {
  if (evaluer.isGood()) {
    const Sls::ALP_set_of_parameters &p = evaluer.parameters();
    out << "# lambda=" << p.lambda << " K=" << p.K << "\n";
  }
}

void LastEvaluer::writeParameters(std::ostream &out) const {
  std::streamsize prec = out.precision(17);
  if (evaluer.isGood()) {
    const Sls::ALP_set_of_parameters &p = evaluer.parameters();
    // !!! flip (I, J) order, to match AlignmentEvaluer::initParameters
    out << p.lambda << ", " << p.K << ",\n"
	<< p.a_J << ", " << p.b_J << ",\n"
	<< p.a_I << ", " << p.b_I << ",\n"
	<< p.alpha_J << ", " << p.beta_J << ",\n"
	<< p.alpha_I << ", " << p.beta_I << ",\n"
	<< p.sigma << ", " << p.tau << "\n";
  }
  out.precision(prec);
}

}
