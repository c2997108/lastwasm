// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "last.hh"

#include "LastalArguments.hh"
#include "QualityPssmMaker.hh"
#include "OneQualityScoreMatrix.hh"
#include "TwoQualityScoreMatrix.hh"
#include "LastEvaluer.hh"
#include "GeneticCode.hh"
#include "AlignmentPot.hh"
#include "Alignment.hh"
#include "SegmentPairPot.hh"
#include "ScoreMatrix.hh"
#include "TantanMasker.hh"
#include "DiagonalTable.hh"
#include "gaplessXdrop.hh"
#include "gaplessPssmXdrop.hh"
#include "gaplessTwoQualityXdrop.hh"
#include "mcf_substitution_matrix_stats.hh"
#include "zio.hh"
#include "stringify.hh"
#include "threadUtil.hh"
#include "split/mcf_last_splitter.hh"

#include <math.h>
#include <signal.h>
#include <stdlib.h>  // EXIT_SUCCESS, EXIT_FAILURE

#include <iomanip>  // setw
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <mutex>

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << args.programName << ": " << x << '\n'
#define LOG2(x) if( args.verbosity > 1 ) std::cerr << args.programName << ": " << x << '\n'

static void warn( const char* programName, const char* s ){
  std::cerr << programName << ": " << s << '\n';
}

using namespace cbrc;

typedef unsigned long long countT;

struct LastAligner {  // data that changes between queries
  Aligners engines;
  LastSplitter splitter;
  std::vector<int> qualityPssm;
  std::vector<AlignmentText> textAlns;
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  countT numOfNormalLetters;
  countT numOfSequences;
};

struct SubstitutionMatrices {
  mcf::SubstitutionMatrixStats stats;
  QualityPssmMaker maker;

  int scores[scoreMatrixRowSize][scoreMatrixRowSize];
  double ratiosMatrix[scoreMatrixRowSize][scoreMatrixRowSize];
  double *ratios[scoreMatrixRowSize];
  OneQualityScoreMatrix oneQual;
  OneQualityExpMatrix oneQualExp;
  TwoQualityScoreMatrix twoQual;

  int scoresMasked[scoreMatrixRowSize][scoreMatrixRowSize];
  double ratiosMaskedMatrix[scoreMatrixRowSize][scoreMatrixRowSize];
  double *ratiosMasked[scoreMatrixRowSize];
  OneQualityScoreMatrix oneQualMasked;
  OneQualityExpMatrix oneQualExpMasked;
  TwoQualityScoreMatrix twoQualMasked;
};

struct SeqData {
  size_t seqNum;
  size_t padLen;
  size_t seqBeg;
  size_t seqEnd;
  size_t frameSize;
  uchar *seq;
  const uchar *seqPadBeg;
  const uchar *seqPadEnd;
  const uchar *qual;
  int *qualityPssm;
  const ScoreMatrixRow *pssm;
};

namespace {
  std::mutex inputMutex;
  std::mutex outputMutex;

  char **querySequenceFileNames;
  mcf::izstream querySequenceFile;
  mcf::izstream querySequenceFile2;  // for files with paired sequences

  LastalArguments args;
  Alphabet alph;
  Alphabet queryAlph;  // for translated alignment
  TantanMasker tantanMasker;
  GeneticCode geneticCode;
  SplitAlignerParams splitParams;
  const unsigned maxNumOfIndexes = 16;
  SubsetSuffixArray suffixArrays[maxNumOfIndexes];
  DnaWordsFinder wordsFinder;
  ScoreMatrix scoreMatrix;
  SubstitutionMatrices fwdMatrices;
  SubstitutionMatrices revMatrices;
  mcf::GapCosts gapCosts;
  std::vector<LastAligner> aligners;
  LastEvaluer evaluer;
  LastEvaluer gaplessEvaluer;
  MultiSequence qrySeqsGlobal;  // sequence that hasn't been indexed by lastdb
  MultiSequence refSeqs;  // sequence that has been indexed by lastdb
  sequenceFormat::Enum referenceFormat = sequenceFormat::fasta;
  int minScoreGapless;
  unsigned numOfVolumes = -1;
  unsigned numOfIndexes = 1;  // assume this value, if unspecified
}

void complementMatrix(const ScoreMatrixRow *from, ScoreMatrixRow *to) {
  for (unsigned i = 0; i < scoreMatrixRowSize; ++i)
    for (unsigned j = 0; j < scoreMatrixRowSize; ++j)
      to[i][j] = from[alph.complement[i]][alph.complement[j]];
}

// Meaningless for PSSMs, unless they have the same scale as the score matrix
static void
calculateSubstitutionScoreMatrixStatistics(const std::string &matrixName) {
  int *scoreMat[scoreMatrixRowSize];
  // the case-sensitivity of the matrix makes no difference here
  std::copy(scoreMatrix.caseSensitive,
	    scoreMatrix.caseSensitive + alph.size, scoreMat);

  mcf::SubstitutionMatrixStats &stats = fwdMatrices.stats;
  if (scoreMatrix.hasLetterFrequencies()) {
    unsigned alphSize2 = scoreMatrix.isCodonCols() ? 64 : alph.size;
    double *p1 = stats.sizedLetterProbs1(alph.size);
    double *p2 = stats.sizedLetterProbs2(alphSize2);
    scoreMatrix.calcLetterProbs(p1, alph.size, p2, alphSize2, alph.encode);
    if (args.temperature < 0) {
      ERR("not implemented");
    } else {
      stats.calcBias(scoreMat, alph.size, alphSize2, args.temperature);
      LOG("score matrix bias=" << stats.bias());
    }
  } else {
    if (scoreMatrix.isCodonCols()) {
      static const char msg[] = "can't calculate probabilities";
      if (args.isSumOfPaths()) ERR(msg);  // xxx what if temperature is given?
      LOG(msg);
      return;
    } else if (args.temperature < 0) {
      const char *canonicalMatrixName = ScoreMatrix::canonicalName(matrixName);
      LOG("calculating matrix probabilities...");
      stats.calcUnbiased(canonicalMatrixName, scoreMat, alph.size);
      if (stats.isBad()) {
	static const char msg[] = "can't calculate probabilities: "
	  "maybe the mismatch costs are too weak";
	if (isUseQuality(args.inputFormat) || args.isSumOfPaths()) ERR(msg);
	LOG(msg);
	return;
      }
    } else {
      stats.calcFromScale(scoreMat, alph.size, args.temperature);
      LOG("score matrix bias=" << stats.bias());
    }
  }

  revMatrices.stats = stats;
  revMatrices.stats.flipDnaStrands(alph.complement);

  const double *p1 = stats.letterProbs1();
  const double *p2 = stats.letterProbs2();

  LOG( "matrix lambda=" << stats.lambda() );
  if (scoreMatrix.isCodonCols()) return;
  LOG( "matrix letter frequencies (upper=reference, lower=query):" );
  if( args.verbosity > 0 ){
    std::cerr << std::left;
    std::streamsize p = std::cerr.precision(2);
    unsigned e = alph.size;
    for( unsigned i = 0; i < e; ++i )
      std::cerr << std::setw(3) << alph.letters[i] << (i + 1 < e ? " " : "\n");
    for( unsigned i = 0; i < e; ++i )
      std::cerr << std::setw(3) << 100 * p1[i] << (i + 1 < e ? " " : "\n");
    for( unsigned i = 0; i < e; ++i )
      std::cerr << std::setw(3) << 100 * p2[i] << (i + 1 < e ? " " : "\n");
    std::cerr.precision(p);
    std::cerr << std::right;
  }
}

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( const std::string& matrixName,
		      const std::string& matrixFile ){
  if( !matrixName.empty() && !args.isGreedy ){
    scoreMatrix.fromString( matrixFile );
    if (scoreMatrix.isCodonRows())
      err("unsupported score matrix");
    if (scoreMatrix.isCodonCols() && !args.isTranslated())
      err("unsuitable score matrix");
  } else {
    scoreMatrix.setMatchMismatch( args.matchScore, args.mismatchCost,
				  alph.letters );
  }

  scoreMatrix.init(alph.encode);
  const mcf::SubstitutionMatrixStats &stats = fwdMatrices.stats;
  if (args.outputType > 0) {
    calculateSubstitutionScoreMatrixStatistics(matrixName);
    if (!stats.isBad()) {
      scoreMatrix.addAmbiguousScores(alph.letters == alph.dna,
				     args.ambiguousLetterOpt % 2,
				     args.ambiguousLetterOpt / 2,
				     alph.encode, stats.lambda(),
				     stats.letterProbs1(),
				     stats.letterProbs2());
      scoreMatrix.init(alph.encode);
    }
  }

  // If the input is a PSSM, the score matrix is not used, and its
  // maximum score should not be used.  Here, we try to set it to a
  // high enough value that it has no effect.  This is a kludge - it
  // would be nice to use the maximum PSSM score.
  if( args.inputFormat == sequenceFormat::pssm ) scoreMatrix.maxScore = 10000;
  // This would work, except the maxDrops aren't finalized yet:
  // maxScore = std::max(args.maxDropGapped, args.maxDropFinal) + 1;

  for (unsigned i = 0; i < scoreMatrixRowSize; ++i) {  // xxx
    for (unsigned j = 0; j < scoreMatrixRowSize; ++j) {
      fwdMatrices.scores[i][j] = scoreMatrix.caseInsensitive[i][j];
      fwdMatrices.scoresMasked[i][j] = scoreMatrix.caseSensitive[i][j];
    }
  }

  if( args.isQueryStrandMatrix && args.strand != 1 ){
    complementMatrix(fwdMatrices.scores, revMatrices.scores);
    complementMatrix(fwdMatrices.scoresMasked, revMatrices.scoresMasked);
  }

  double s = stats.lambda();
  for (unsigned i = 0; i < scoreMatrixRowSize; ++i) {
    fwdMatrices.ratios[i] = fwdMatrices.ratiosMatrix[i];
    revMatrices.ratios[i] = revMatrices.ratiosMatrix[i];
    fwdMatrices.ratiosMasked[i] = fwdMatrices.ratiosMaskedMatrix[i];
    revMatrices.ratiosMasked[i] = revMatrices.ratiosMaskedMatrix[i];
    for (unsigned j = 0; j < scoreMatrixRowSize; ++j) {
      fwdMatrices.ratios[i][j] = exp(s * fwdMatrices.scores[i][j]);
      revMatrices.ratios[i][j] = exp(s * revMatrices.scores[i][j]);
      fwdMatrices.ratiosMasked[i][j] = exp(s * fwdMatrices.scoresMasked[i][j]);
      revMatrices.ratiosMasked[i][j] = exp(s * revMatrices.scoresMasked[i][j]);
    }
  }
}

void makeQualityScorers(SubstitutionMatrices &m, bool isCheck) {
  bool isMatchMismatch = (args.matrixFile.empty() && args.matchScore > 0);
  bool isPhred1 = isPhred( referenceFormat );
  int offset1 = qualityOffset( referenceFormat );
  bool isPhred2 = isPhred( args.inputFormat );
  int offset2 = qualityOffset( args.inputFormat );
  const uchar *toUnmasked = alph.numbersToUppercase;

  if( !isUseFastq( referenceFormat ) ){
    if( isUseFastq( args.inputFormat ) ){
      if (isCheck) LOG("calculating per-quality scores...");
      if (args.maskLowercase > 0) {
	m.oneQualMasked.init(m.scores, alph.size, m.stats,
			     isPhred2, offset2, toUnmasked, true);
	if (args.isSumOfPaths())
	  m.oneQualExpMasked.init(m.oneQualMasked, args.temperature);
      }
      if (args.maskLowercase < 3) {
	m.oneQual.init(m.scores, alph.size, m.stats,
		       isPhred2, offset2, toUnmasked, false);
	if (args.isSumOfPaths())
	  m.oneQualExp.init(m.oneQual, args.temperature);
      }
      const OneQualityScoreMatrix &q = (args.maskLowercase < 3) ?
	m.oneQual : m.oneQualMasked;
      if( isCheck && args.verbosity > 0 )
	writeOneQualityScoreMatrix( q, alph.letters.c_str(),
				    offset2, std::cerr );
    }
    if( isUseQuality(args.inputFormat) ){
      m.maker.init(m.scores, alph.size, m.stats.lambda(),
		   isMatchMismatch, args.matchScore, -args.mismatchCost,
		   isPhred2, offset2, toUnmasked);
    }
  } else {
    if( isUseFastq( args.inputFormat ) ){
      if( args.maskLowercase > 0 )
	m.twoQualMasked.init(m.scores, m.stats,
			     isPhred1, offset1, isPhred2, offset2,
			     toUnmasked, true, isMatchMismatch);
      if( args.maskLowercase < 3 )
	m.twoQual.init(m.scores, m.stats,
		       isPhred1, offset1, isPhred2, offset2,
		       toUnmasked, false, isMatchMismatch);
      if (args.isSumOfPaths())
        ERR("fastq-versus-fastq alignment probabilities not implemented");
    } else if (isCheck) {
      warn(args.programName,
	   "quality data not used for non-fastq query versus fastq reference");
    }
  }
}

void makeQualityScorers(){
  if( args.isGreedy ) return;

  if( args.isTranslated() )
    if( isUseQuality( args.inputFormat ) || isUseQuality( referenceFormat ) )
      return warn( args.programName,
		   "quality data not used for DNA-versus-protein alignment" );

  makeQualityScorers(fwdMatrices, true);
  if (args.isQueryStrandMatrix && args.strand != 1) {
    makeQualityScorers(revMatrices, false);
  }
}

// Calculate statistical parameters for the alignment scoring scheme
static void calculateScoreStatistics(const std::string& matrixName,
				     countT refLetters, countT refMaxSeqLen) {
  const mcf::SubstitutionMatrixStats &stats = fwdMatrices.stats;
  if (stats.isBad()) return;
  const char *canonicalMatrixName = ScoreMatrix::canonicalName( matrixName );
  if (args.temperature > 0 && !matrixName.empty()) canonicalMatrixName = " ";
  bool isGapped = (args.outputType > 1);
  const ScoreMatrixRow *scoreMat = fwdMatrices.scoresMasked;
  const double *p1 = stats.letterProbs1();
  const double *p2 = stats.letterProbs2();
  if (args.isTranslated() && !scoreMatrix.isCodonCols()) p2 = 0;
  int fsCost = gapCosts.isNewFrameshifts() ? 0 : gapCosts.frameshiftCost;
  LOG( "getting E-value parameters..." );
  try{
    gaplessEvaluer.init(0, 0, 0, alph.letters.c_str(), scoreMat, p1, p2, false,
			0, 0, 0, 0, fsCost, geneticCode, 0, 0);
    if (args.scoreType != 0 && isGapped) {
      if (args.temperature < 0 &&
	  gapCosts.delProbPieces[0].openProb +
	  gapCosts.insProbPieces[0].openProb > 0) return;
      unsigned alphSize2 = scoreMatrix.isCodonCols() ? 64 : alph.size;
      evaluer.initFullScores(fwdMatrices.ratios, p1, alph.size,
			     stats.letterProbs2(), alphSize2,
			     gapCosts, stats.lambda(),
			     args.gumbelSimAlignmentCount,
			     args.gumbelSimSequenceLength, args.verbosity,
			     gapCosts.isNewFrameshifts());
    } else {
      const mcf::GapCosts::Piece &del = gapCosts.delPieces[0];
      const mcf::GapCosts::Piece &ins = gapCosts.insPieces[0];
      evaluer.init(canonicalMatrixName, args.matchScore, args.mismatchCost,
		   alph.letters.c_str(), scoreMat, p1, p2, isGapped,
		   del.openCost, del.growCost, ins.openCost, ins.growCost,
		   fsCost, geneticCode,
		   args.geneticCodeFile.c_str(), args.verbosity);
    }
    countT m = std::min(refMaxSeqLen, refLetters);
    evaluer.setSearchSpace(refLetters, m, args.numOfStrands());
    if( args.verbosity > 0 ) evaluer.writeParameters( std::cerr );
  }catch( const Sls::error& e ){
    LOG( "can't get E-value parameters for this scoring scheme" );
  }
}

struct LastdbData {
  countT numOfSeqs;
  countT numOfLetters;
  countT maxSeqLen;
  size_t minSeedLimit;
  size_t minimizerWindow;
  bool isKeepLowercase;
  int tantanSetting;
  int maxRepeatUnit;
  int isCaseSensitive;
  int strand;
  int bitsPerBase;
  int bitsPerInt;
};

// Read the .prj file for the whole database
static LastdbData readOuterPrj(const std::string &fileName) {
  countT z = -1;
  LastdbData d = {z, z, z, 0, 1, true, 0, 0, -1, 1, CHAR_BIT, 0};
  int version = 0;
  std::string alphabetLetters;

  std::ifstream f(fileName.c_str());
  if (!f) ERR("can't open file: " + fileName);

  std::string trigger = "#lastal";
  std::string line, word;
  while (getline(f, line)) {
    if (line.compare(0, trigger.size(), trigger) == 0) {
      args.fromLine(line);
      continue;
    }
    std::istringstream iss(line);
    getline(iss, word, '=');
    if (word == "version") iss >> version;
    if (word == "alphabet") iss >> alphabetLetters;
    if (word == "strand") iss >> d.strand;
    if (word == "numofsequences") iss >> d.numOfSeqs;
    if (word == "numofletters") iss >> d.numOfLetters;
    if (word == "maxsequenceletters") iss >> d.maxSeqLen;
    if (word == "maxunsortedinterval") iss >> d.minSeedLimit;
    if (word == "keeplowercase") iss >> d.isKeepLowercase;
    if (word == "tantansetting") iss >> d.tantanSetting;
    if (word == "maxrepeatunit") iss >> d.maxRepeatUnit;
    if (word == "masklowercase") iss >> d.isCaseSensitive;
    if (word == "sequenceformat") iss >> referenceFormat;
    if (word == "minimizerwindow") iss >> d.minimizerWindow;
    if (word == "volumes") iss >> numOfVolumes;
    if (word == "numofindexes") iss >> numOfIndexes;
    if (word == "integersize") iss >> d.bitsPerInt;
    if (word == "symbolsize") iss >> d.bitsPerBase;
  }

  if (f.eof() && !f.bad()) f.clear();
  if (alphabetLetters.empty() || d.numOfSeqs+1 == 0 || d.numOfLetters+1 == 0 ||
      (d.maxSeqLen == 0 && d.numOfLetters != 0) ||
      d.isCaseSensitive < 0 || numOfIndexes > maxNumOfIndexes ||
      referenceFormat == sequenceFormat::prb ||
      referenceFormat == sequenceFormat::pssm) {
    f.setstate(std::ios::failbit);
  }
  if (!f) ERR("can't read file: " + fileName);
  if (version < 294 && version > 0) {
    ERR("the lastdb files are old: please re-run lastdb");
  }

  if (d.bitsPerInt < 1 && version < 999) d.bitsPerInt = 32;
  alph.init(alphabetLetters, d.bitsPerBase == 4);
  return d;
}

// Read a per-volume .prj file, with info about a database volume
void readInnerPrj(const std::string &fileName,
		  size_t &seqCount, size_t &seqLen) {
  std::ifstream f( fileName.c_str() );
  if( !f ) ERR( "can't open file: " + fileName );

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "numofsequences" ) iss >> seqCount;
    if( word == "numofletters" ) iss >> seqLen;
    if( word == "numofindexes" ) iss >> numOfIndexes;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( seqCount+1 == 0 || seqLen+1 == 0 || numOfIndexes > maxNumOfIndexes ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) ERR( "can't read file: " + fileName );
}

static size_t seedSearchEnd(size_t seqEnd) {
  size_t d = wordsFinder.wordLength ? wordsFinder.wordLength : 1;
  size_t x = args.minHitDepth - std::min(d, args.minHitDepth);
  return seqEnd - std::min(x, seqEnd);
}

// Write match counts for each query sequence
void writeCounts(const std::vector< std::vector<countT> > &matchCounts,
		 const MultiSequence &qrySeqs, size_t firstSequence) {
  for (size_t i = 0; i < matchCounts.size(); ++i) {
    std::cout << qrySeqs.seqName(firstSequence + i) << '\n';
    for (size_t j = args.minHitDepth; j < matchCounts[i].size(); ++j) {
      std::cout << j << '\t' << matchCounts[i][j] << '\n';
    }
    std::cout << '\n';  // blank line afterwards
  }
}

// Count all matches, of all sizes, of a query sequence against a suffix array
void countMatches(std::vector<countT> &counts, const SeqData &qryData) {
  if (wordsFinder.wordLength) {  // YAGNI
    err("can't count initial matches with word-restricted seeds, sorry");
  }
  size_t loopEnd = seedSearchEnd(qryData.seqEnd);

  for (size_t i = qryData.seqBeg; i < loopEnd; i += args.queryStep) {
    for (unsigned x = 0; x < numOfIndexes; ++x) {
      suffixArrays[x].countMatches(counts, qryData.seq + i,
				   refSeqs.seqPtr(), 0, args.maxHitDepth);
    }
  }
}

static int *qualityPssmSpace(LastAligner &aligner, size_t padLen) {
  if (args.outputType == 0 || args.isGreedy || args.isTranslated() ||
      !isUseQuality(args.inputFormat) || isUseQuality(referenceFormat)) {
    return 0;
  }
  aligner.qualityPssm.resize(padLen * scoreMatrixRowSize);
  return &aligner.qualityPssm[0];
}

static const ScoreMatrixRow *getQueryPssm(const int *qualityPssm,
					  const MultiSequence &qrySeqs,
					  size_t padBeg) {
  if (args.isGreedy) return 0;

  if (args.inputFormat == sequenceFormat::pssm) {
    return qrySeqs.pssmReader() + padBeg;
  }

  return reinterpret_cast<const ScoreMatrixRow *>(qualityPssm);
}

namespace Phase{ enum Enum{ gapless, pregapped, gapped, postgapped }; }

static bool isMaskLowercase(Phase::Enum e) {
  return (e < 1 && args.maskLowercase > 0)
    || (e < 3 && args.maskLowercase > 1 && args.scoreType != 0)
    || args.maskLowercase > 2;
}

struct Dispatcher{
  BigSeq a;  // the reference sequence
  const uchar* b;  // the query sequence
  const uchar* i;  // the reference quality data
  const uchar* j;  // the query quality data
  const ScoreMatrixRow* p;  // the query PSSM
  const ScoreMatrixRow* m;  // the score matrix
  const const_dbl_ptr* r;   // the substitution probability ratios
  const TwoQualityScoreMatrix& t;
  int d;  // the maximum score drop
  int z;

  Dispatcher(Phase::Enum e, const SeqData &qryData,
	     const SubstitutionMatrices &matrices) :
      a( refSeqs.seqPtr() ),
      b( qryData.seq ),
      i( refSeqs.qualityReader() ),
      j( qryData.qual ),
      p( qryData.pssm ),
      m( isMaskLowercase(e) ? matrices.scoresMasked : matrices.scores ),
      r( isMaskLowercase(e) ? matrices.ratiosMasked : matrices.ratios ),
      t( isMaskLowercase(e) ? matrices.twoQualMasked : matrices.twoQual ),
      d( (e == Phase::gapless) ? args.maxDropGapless :
         (e == Phase::pregapped ) ? args.maxDropGapped : args.maxDropFinal ),
      z( t ? 2 : p ? 1 : 0 ){}

  int gaplessOverlap(size_t x, size_t y, size_t &rev, size_t &fwd) const {
    if (z==0) return gaplessXdropOverlap(a+x, b+y, m, d, rev, fwd);
    if (z==1) return gaplessPssmXdropOverlap(a+x, p+y, d, rev, fwd);
    return gaplessTwoQualityXdropOverlap(a+x, i+x, b+y, j+y, t, d, rev, fwd);
  }

  void gaplessExtensionScores(size_t rPos, size_t qPos,
			      int &fwdScore, int &revScore) const {
    if (z == 0) {
      gaplessXdropScores(a+rPos, b+qPos, m, d, fwdScore, revScore);
    } else if (z == 1) {
      gaplessPssmXdropScores(a+rPos, p+qPos, d, fwdScore, revScore);
    } else {
      gaplessTwoQualityXdropScores(a+rPos, i+rPos, b+qPos, j+qPos, t, d,
				   fwdScore, revScore);
    }
  }

  bool gaplessEnds(int fwdScore, int revScore,
		   size_t &rPos, size_t &qPos, size_t &length) const {
    return (z == 0) ? gaplessXdropEnds(a, b, m, d, fwdScore, revScore,
				       rPos, qPos, length)
      :    (z == 1) ? gaplessPssmXdropEnds(a, p, d, fwdScore, revScore,
					   rPos, qPos, length)
      :               gaplessTwoQualityXdropEnds(a, i, b, j, t, d, fwdScore,
						 revScore, rPos, qPos, length);
  }

  int gaplessScore(size_t x, size_t y, size_t length) const {
    if (z==0) return gaplessAlignmentScore(a+x, b+y, m, length);
    if (z==1) return gaplessPssmAlignmentScore(a+x, p+y, length);
    return gaplessTwoQualityAlignmentScore(a+x, i+x, b+y, j+y, t, length);
  }
};

static bool isCollatedAlignments() {
  return args.outputFormat == 'b' || args.outputFormat == 'B' ||
    args.cullingLimitForFinalAlignments + 1 || numOfVolumes > 1;
}

static void writeAlignment(LastAligner &aligner, const MultiSequence &qrySeqs,
			   const SeqData &qryData, const Alignment &aln,
			   const AlignmentExtras &extras = AlignmentExtras()) {
  int translationType = scoreMatrix.isCodonCols() ? 2 : args.isTranslated();
  AlignmentText a = aln.write(refSeqs, qrySeqs, qryData.seqNum, qryData.seq,
			      alph, queryAlph,
			      translationType, geneticCode.getCodonToAmino(),
			      evaluer, args.outputFormat, extras);
  if (isCollatedAlignments() || aligners.size() > 1 || args.isSplit) {
    aligner.textAlns.push_back(a);
  } else {
    std::cout << a.text;
    delete[] a.text;
  }
}

static void writeSegmentPair(LastAligner &aligner,
			     const MultiSequence &qrySeqs,
			     const SeqData &qryData, const SegmentPair &s) {
  Alignment a;
  a.fromSegmentPair(s);
  writeAlignment(aligner, qrySeqs, qryData, a);
}

struct GaplessAlignmentCounts {
  countT matchCount;
  countT gaplessExtensionCount;
  countT gaplessAlignmentCount;
  size_t maxSignificantAlignments;
};

// Get seed hits and gapless alignments at one query-sequence position
void alignGapless1(LastAligner &aligner, SegmentPairPot &gaplessAlns,
		   const MultiSequence &qrySeqs, const SeqData &qryData,
		   const Dispatcher &dis, DiagonalTable &dt,
		   GaplessAlignmentCounts &counts, const SubsetSuffixArray &sa,
		   const uchar *qryPtr, unsigned seedNum) {
  const bool isOverlap = (args.globality && args.outputType == 1);

  size_t beg;
  size_t end;
  sa.match(beg, end, qryPtr, dis.a, seedNum,
	   args.oneHitMultiplicity, args.minHitDepth, args.maxHitDepth);
  counts.matchCount += end - beg;

  size_t qryPos = qryPtr - dis.b;  // coordinate in the query sequence
  size_t maxAlignments = args.maxGaplessAlignmentsPerQueryPosition;

  for (/* noop */; beg < end; ++beg) {
    if (maxAlignments == 0) break;

    // it might be faster to unpack all these refPos values at once:
    size_t refPos = sa.getPosition(beg);  // position in the reference sequence
    size_t diagonal = qryPos - refPos;
    if (dt.isCovered(diagonal, qryPos)) continue;
    ++counts.gaplessExtensionCount;
    int score;

    if (isOverlap) {
      size_t revLen, fwdLen;
      score = dis.gaplessOverlap(refPos, qryPos, revLen, fwdLen);
      if (score < minScoreGapless) continue;
      SegmentPair sp(refPos - revLen, qryPos - revLen, revLen + fwdLen, score);
      dt.addEndpoint(diagonal, sp.end2());
      writeSegmentPair(aligner, qrySeqs, qryData, sp);
    } else {
      int fwdScore, revScore;
      dis.gaplessExtensionScores(refPos, qryPos, fwdScore, revScore);
      score = fwdScore + revScore;
      if (score < minScoreGapless) continue;
      size_t rPos = refPos;
      size_t qPos = qryPos;
      size_t length;
      if (!dis.gaplessEnds(fwdScore, revScore, rPos, qPos, length)) continue;
      SegmentPair sp(rPos, qPos, length, score);
      dt.addEndpoint(diagonal, sp.end2());

      if (args.outputType == 1) {  // we just want gapless alignments
	writeSegmentPair(aligner, qrySeqs, qryData, sp);
      } else {
	gaplessAlns.add(sp);
      }
    }

    --maxAlignments;
    ++counts.gaplessAlignmentCount;

    if (score >= args.minScoreGapped && --counts.maxSignificantAlignments == 0)
      break;
  }
}

// Find query matches to the suffix array, and do gapless extensions
void alignGapless(LastAligner &aligner, SegmentPairPot &gaplessAlns,
		  const MultiSequence &qrySeqs, const SeqData &qryData,
		  const Dispatcher &dis) {
  DiagonalTable dt;  // record already-covered positions on each diagonal
  size_t maxAlignments =
    args.maxAlignmentsPerQueryStrand ? args.maxAlignmentsPerQueryStrand : 1;
  GaplessAlignmentCounts counts = {0, 0, 0, maxAlignments};

  size_t loopBeg = qryData.seqBeg;
  size_t loopEnd = seedSearchEnd(qryData.seqEnd);

  const uchar *querySeq = qryData.seq;
  const uchar *qryBeg = querySeq + loopBeg;
  const uchar *qryEnd = querySeq + loopEnd;

  const unsigned wordLen = wordsFinder.wordLength;

  if (wordLen) {
    unsigned hash = 0;
    qryBeg = wordsFinder.init(qryBeg, qryEnd, &hash);
    while (qryBeg < qryEnd) {
      unsigned c = wordsFinder.baseToCode[*qryBeg];
      ++qryBeg;
      if (c != dnaWordsFinderNull) {
	unsigned w = wordsFinder.next(&hash, c);
	if (w != dnaWordsFinderNull) {
	  alignGapless1(aligner, gaplessAlns, qrySeqs, qryData, dis, dt,
			counts, suffixArrays[0], qryBeg - wordLen, w);
	  if (counts.maxSignificantAlignments == 0) break;
	}
      } else {
	qryBeg = wordsFinder.init(qryBeg, qryEnd, &hash);
      }
    }
  } else {
    std::vector<SubsetMinimizerFinder> minFinders(numOfIndexes);
    for (unsigned x = 0; x < numOfIndexes; ++x) {
      minFinders[x].init(suffixArrays[x].getSeeds()[0], qryBeg, qryEnd);
    }
    const size_t step = args.queryStep;
    const size_t w = args.minimizerWindow;
    for (size_t i = loopBeg; i < loopEnd; i += step) {
      const uchar *qryPtr = querySeq + i;
      for (unsigned x = 0; x < numOfIndexes; ++x) {
	if (w < 2 || minFinders[x].isMinimizer(suffixArrays[x].getSeeds()[0],
					       qryPtr, qryEnd, w)) {
	  alignGapless1(aligner, gaplessAlns, qrySeqs, qryData, dis, dt,
			counts, suffixArrays[x], qryPtr, 0);
	}
      }
      if (counts.maxSignificantAlignments == 0) break;
    }
  }

  LOG2( "initial matches=" << counts.matchCount );
  LOG2( "gapless extensions=" << counts.gaplessExtensionCount );
  LOG2( "gapless alignments=" << counts.gaplessAlignmentCount );
}

// Shrink the SegmentPair to its longest run of identical matches.
// This trims off possibly unreliable parts of the gapless alignment.
// It may not be the best strategy for protein alignment with subset
// seeds: there could be few or no identical matches...
void shrinkToLongestIdenticalRun( SegmentPair& sp, const Dispatcher& dis ){
  const uchar *map2 = scoreMatrix.isCodonCols() ?
    geneticCode.getCodonToAmino() : alph.numbersToUppercase;
  sp.maxIdenticalRun(dis.a, dis.b, alph.numbersToUppercase, map2);
  sp.score = dis.gaplessScore(sp.beg1(), sp.beg2(), sp.size);
}

// Do gapped extensions of the gapless alignments
void alignGapped(LastAligner &aligner, AlignmentPot &gappedAlns,
		 SegmentPairPot &gaplessAlns, const SeqData &qryData,
		 const SubstitutionMatrices &matrices, Phase::Enum phase) {
  Dispatcher dis(phase, qryData, matrices);
  countT gappedExtensionCount = 0, gappedAlignmentCount = 0;

  // Redo the gapless extensions, using gapped score parameters.
  // Without this, if we self-compare a huge sequence, we risk getting
  // huge gapped extensions.
  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    SegmentPair& sp = gaplessAlns.items[i];
    int fwdScore, revScore;
    dis.gaplessExtensionScores(sp.beg1(), sp.beg2(), fwdScore, revScore);
    size_t beg1 = sp.beg1();
    size_t beg2 = sp.beg2();
    size_t length;
    if (dis.gaplessEnds(fwdScore, revScore, beg1, beg2, length)) {
      sp = SegmentPair(beg1, beg2, length, fwdScore + revScore);
    } else {
      SegmentPairPot::mark(sp);
    }
  }

  erase_if( gaplessAlns.items, SegmentPairPot::isMarked );

  gaplessAlns.cull( args.cullingLimitForGaplessAlignments );
  gaplessAlns.sort();  // sort by score descending, and remove duplicates

  LOG2( "redone gapless alignments=" << gaplessAlns.size() );

  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    SegmentPair& sp = gaplessAlns.get(i);

    if( SegmentPairPot::isMarked(sp) ) continue;

    Alignment aln;
    AlignmentExtras extras;  // not used
    aln.seed = sp;

    shrinkToLongestIdenticalRun( aln.seed, dis );

    // do gapped extension from each end of the seed:
    aln.makeXdrop(aligner.engines, args.isGreedy, args.scoreType,
		  dis.a, dis.b, args.globality,
		  dis.m, scoreMatrix.maxScore, scoreMatrix.minScore,
		  dis.r, matrices.stats.lambda(), gapCosts, dis.d,
		  qryData.frameSize, dis.p, dis.t, dis.i, dis.j, alph, extras);
    ++gappedExtensionCount;

    if( aln.score < args.minScoreGapped ) continue;

    if (args.scoreType == 0 &&
	!aln.isOptimal(dis.a, dis.b, args.globality, dis.m, dis.d, gapCosts,
		       qryData.frameSize, dis.p, dis.t, dis.i, dis.j)) {
      // If retained, non-"optimal" alignments can hide "optimal"
      // alignments, e.g. during non-redundantization.
      continue;
    }

    gaplessAlns.markAllOverlaps( aln.blocks );
    gaplessAlns.markTandemRepeats( aln.seed, args.maxRepeatDistance );

    if (phase == Phase::gapped) gappedAlns.add(aln);
    else SegmentPairPot::markAsGood(sp);

    ++gappedAlignmentCount;
    if( gappedAlignmentCount >= args.maxAlignmentsPerQueryStrand ) break;
  }

  LOG2( "gapped extensions=" << gappedExtensionCount );
  LOG2( "gapped alignments=" << gappedAlignmentCount );
}

// Redo gapped extensions, but keep the old alignment scores
static void alignPostgapped(LastAligner &aligner, AlignmentPot &gappedAlns,
			    size_t frameSize, const Dispatcher &dis) {
  AlignmentExtras extras;  // not used
  for (size_t i = 0; i < gappedAlns.size(); ++i) {
    Alignment &aln = gappedAlns.items[i];
    aln.makeXdrop(aligner.engines, args.isGreedy, args.scoreType,
		  dis.a, dis.b, args.globality,
		  dis.m, scoreMatrix.maxScore, scoreMatrix.minScore,
		  0, 0, gapCosts, dis.d,
		  frameSize, dis.p, dis.t, dis.i, dis.j, alph, extras);
  }
}

// Print the gapped alignments, after optionally calculating match
// probabilities and re-aligning using the gamma-centroid algorithm
void alignFinish(LastAligner &aligner, const MultiSequence &qrySeqs,
		 const SeqData &qryData, const AlignmentPot &gappedAlns,
		 const SubstitutionMatrices &matrices, const Dispatcher &dis) {
  for( size_t i = 0; i < gappedAlns.size(); ++i ){
    const Alignment& aln = gappedAlns.items[i];
    AlignmentExtras extras;
    if (args.scoreType != 0) extras.fullScore = -1;  // score is fullScore
    if( args.outputType < 4 ){
      writeAlignment(aligner, qrySeqs, qryData, aln, extras);
    } else {  // calculate match probabilities:
      Alignment probAln;
      probAln.seed = aln.seed;
      probAln.makeXdrop(aligner.engines, args.isGreedy, args.scoreType,
			dis.a, dis.b, args.globality,
			dis.m, scoreMatrix.maxScore, scoreMatrix.minScore,
			dis.r, matrices.stats.lambda(), gapCosts, dis.d,
			qryData.frameSize, dis.p, dis.t, dis.i, dis.j, alph,
			extras, args.gamma, args.outputType);
      assert(aln.score != -INF);
      if (args.maskLowercase == 2 && args.scoreType != 0)
	probAln.score = aln.score;
      writeAlignment(aligner, qrySeqs, qryData, probAln, extras);
    }
  }
}

static void eraseWeakAlignments(AlignmentPot &gappedAlns,
				size_t frameSize, const Dispatcher &dis) {
  for (size_t i = 0; i < gappedAlns.size(); ++i) {
    Alignment &a = gappedAlns.items[i];
    if (!a.hasGoodSegment(dis.a, dis.b, ceil(args.minScoreGapped), dis.m,
			  gapCosts, frameSize, dis.p, dis.t, dis.i, dis.j)) {
      AlignmentPot::mark(a);
    }
  }
  erase_if(gappedAlns.items, AlignmentPot::isMarked);
}

static bool lessForCulling(const AlignmentText &x, const AlignmentText &y) {
  if (x.strandNum != y.strandNum) return x.strandNum < y.strandNum;
  if (x.queryBeg  != y.queryBeg ) return x.queryBeg  < y.queryBeg;
  else                            return x.score     > y.score;
}

// Remove any alignment whose query range overlaps an alignment with
// higher score (and on the same strand):
static void cullOverlappingAlignments(std::vector<AlignmentText> &textAlns,
				      size_t start) {
  size_t end = textAlns.size();
  size_t i = start;
  for (size_t j = start; j < end; ++j) {
    AlignmentText &x = textAlns[j];
    for (size_t k = j + 1; k < end; ++k) {
      AlignmentText &y = textAlns[k];
      if (y.strandNum > x.strandNum || y.queryBeg >= x.queryEnd) break;
      if (x.score > y.score) {
	delete[] y.text;
	y.text = 0;
      }
      if (y.score > x.score) {
	delete[] x.text;
	x.text = 0;
      }
    }
    if (x.text) {
      textAlns[i] = x;
      ++i;
    }
  }
  textAlns.resize(i);
}

// Remove any alignment whose query range lies in LIMIT or more other
// alignments with higher score (and on the same strand):
static void cullFinalAlignments(std::vector<AlignmentText> &textAlns,
				size_t start, size_t limit) {
  if (limit + 1 == 0) return;
  sort(textAlns.begin() + start, textAlns.end(), lessForCulling);
  if (limit == 0) return cullOverlappingAlignments(textAlns, start);
  std::vector<size_t> stash;  // alignments that might dominate subsequent ones
  size_t i = start;  // number of kept alignments so far
  for (size_t j = start; j < textAlns.size(); ++j) {
    AlignmentText &x = textAlns[j];
    size_t numOfDominators = 0;  // number of alignments that dominate x
    size_t a = 0;  // number of kept stash-items so far
    for (size_t b = 0; b < stash.size(); ++b) {
      size_t k = stash[b];
      AlignmentText &y = textAlns[k];
      if (y.strandNum < x.strandNum) break;  // drop the stash
      if (y.queryEnd <= x.queryBeg) continue;  // drop this stash-item
      stash[a++] = k;  // keep this stash-item
      if (y.queryEnd >= x.queryEnd && y.score > x.score) ++numOfDominators;
    }
    stash.resize(a);
    if (numOfDominators >= limit) {
      delete[] x.text;
    } else {
      stash.push_back(i);
      textAlns[i++] = x;  // keep this alignment
    }
  }
  textAlns.resize(i);
}

static void printAlignments(const std::vector<AlignmentText> &textAlns) {
  for (size_t i = 0; i < textAlns.size(); ++i) {
    std::cout << textAlns[i].text;
  }
}

static void clearAlignments(std::vector<AlignmentText> &textAlns) {
  for (size_t i = 0; i < textAlns.size(); ++i) {
    delete[] textAlns[i].text;
  }
  textAlns.clear();
}

void makeQualityPssm(const SeqData &qryData,
		     const SubstitutionMatrices &matrices, bool isMask) {
  int *pssm = qryData.qualityPssm;
  if (!pssm) return;
  const uchar *seqBeg = qryData.seq;
  const uchar *seqEnd = seqBeg + qryData.padLen;

  if (args.inputFormat == sequenceFormat::prb) {
    matrices.maker.make(seqBeg, seqEnd, qryData.qual, pssm, isMask);
  } else {
    const OneQualityScoreMatrix &m =
      isMask ? matrices.oneQualMasked : matrices.oneQual;
    makePositionSpecificScoreMatrix(m, seqBeg, seqEnd, qryData.qual, pssm);
  }
}

static void unmaskLowercase(const SeqData &qryData,
			    const SubstitutionMatrices &matrices) {
  makeQualityPssm(qryData, matrices, false);
  if (scoreMatrix.isCodonCols()) {
    geneticCode.translateWithoutMasking(qryData.seqPadBeg,
					qryData.seqPadEnd, qryData.seq);
  }
}

static void remaskLowercase(const SeqData &qryData,
			    const SubstitutionMatrices &matrices) {
  makeQualityPssm(qryData, matrices, true);
  if (scoreMatrix.isCodonCols()) {
    geneticCode.translate(qryData.seqPadBeg, qryData.seqPadEnd, qryData.seq);
  }
}

// Scan one query sequence against one database volume
void scan(LastAligner &aligner, const MultiSequence &qrySeqs,
	  const SeqData &qryData, const SubstitutionMatrices &matrices) {
  const int maskMode = args.maskLowercase;
  makeQualityPssm(qryData, matrices, maskMode > 0);

  Dispatcher dis0(Phase::gapless, qryData, matrices);
  SegmentPairPot gaplessAlns;
  alignGapless(aligner, gaplessAlns, qrySeqs, qryData, dis0);
  if( args.outputType == 1 ) return;  // we just want gapless alignments
  if( gaplessAlns.size() == 0 ) return;

  if (maskMode == 1 || (maskMode == 2 && args.scoreType == 0))
    unmaskLowercase(qryData, matrices);

  size_t qryLen = qryData.padLen;
  AlignmentPot gappedAlns;
  Centroid &centroid = aligner.engines.centroid;

  if (args.scoreType != 0 && dis0.p) {
    const OneQualityExpMatrix &m =
      (maskMode < 2) ? matrices.oneQualExp : matrices.oneQualExpMasked;
    centroid.setPssm(dis0.p, qryLen, args.temperature, m, dis0.b, dis0.j);
  }

  if (args.maxDropFinal != args.maxDropGapped) {
    alignGapped(aligner, gappedAlns, gaplessAlns, qryData, matrices,
		Phase::pregapped);
    erase_if( gaplessAlns.items, SegmentPairPot::isNotMarkedAsGood );
  }

  alignGapped(aligner, gappedAlns, gaplessAlns, qryData, matrices,
	      Phase::gapped);
  if( gappedAlns.size() == 0 ) return;

  Dispatcher dis3(Phase::postgapped, qryData, matrices);

  if (maskMode == 2 && args.scoreType != 0) {
    unmaskLowercase(qryData, matrices);
    alignPostgapped(aligner, gappedAlns, qryData.frameSize, dis3);
  }

  if (maskMode == 2 && args.scoreType == 0) {
    remaskLowercase(qryData, matrices);
    eraseWeakAlignments(gappedAlns, qryData.frameSize, dis0);
    LOG2("lowercase-filtered alignments=" << gappedAlns.size());
    if (gappedAlns.size() == 0) return;
    if (args.outputType > 3)
      unmaskLowercase(qryData, matrices);
  }

  if( args.outputType > 2 ){  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG2( "nonredundant gapped alignments=" << gappedAlns.size() );
  }

  if (args.outputType > 3 && dis3.p) {
    const OneQualityExpMatrix &m =
      (maskMode < 3) ? matrices.oneQualExp : matrices.oneQualExpMasked;
    centroid.setPssm(dis3.p, qryLen, args.temperature, m, dis3.b, dis3.j);
    if (args.outputType == 7) {
      centroid.setLetterProbsPerPosition(alph.size, qryLen, dis3.b, dis3.j,
					 isUseFastq(args.inputFormat),
					 matrices.maker.qualToProbRight(),
					 matrices.stats.letterProbs2(),
					 alph.numbersToUppercase);
    }
  }

  if (!isCollatedAlignments()) gappedAlns.sort();  // sort by score
  alignFinish(aligner, qrySeqs, qryData, gappedAlns, matrices, dis3);
}

static void tantanMaskOneQuery(const SeqData &qryData) {
  tantanMasker.mask(qryData.seq + qryData.seqBeg, qryData.seq + qryData.seqEnd,
		    queryAlph.numbersToLowercase);
}

static void tantanMaskTranslatedQuery(const SeqData &qryData) {
  size_t frameSize = qryData.padLen / 3;
  size_t dnaBeg = qryData.seqBeg;
  size_t dnaLen = qryData.seqEnd - qryData.seqBeg;
  for (int frame = 0; frame < 3; ++frame) {
    if (dnaLen < 3) break;
    size_t aaBeg = dnaToAa(dnaBeg++, frameSize);
    size_t aaLen = dnaLen-- / 3;
    size_t aaEnd = aaBeg + aaLen;
    tantanMasker.mask(qryData.seq + aaBeg, qryData.seq + aaEnd,
		      alph.numbersToLowercase);
  }
}

// Scan one query sequence strand against one database volume,
// after optionally translating and/or masking the query
void translateAndScan(LastAligner &aligner, MultiSequence &qrySeqs,
		      SeqData &qryData, size_t chunkQryNum,
		      size_t finalCullingLimit,
		      const SubstitutionMatrices &matrices) {
  std::vector<uchar> modifiedQuery;

  if (args.isTranslated()) {
    if (args.tantanSetting && scoreMatrix.isCodonCols()) {
      if (args.isKeepLowercase) {
	err("can't keep lowercase & find simple repeats & use codons");
      }
      tantanMaskOneQuery(qryData);
    }
    modifiedQuery.resize(qryData.padLen);
    qryData.seq = &modifiedQuery[0];
    geneticCode.translate(qryData.seqPadBeg, qryData.seqPadEnd, qryData.seq);
    if (args.tantanSetting && !scoreMatrix.isCodonCols()) {
      tantanMaskTranslatedQuery(qryData);
    }
  } else {
    if (args.tantanSetting) {
      if (args.isKeepLowercase) {
	modifiedQuery.assign(qryData.seqPadBeg, qryData.seqPadEnd);
	qryData.seq = &modifiedQuery[0];
      }
      tantanMaskOneQuery(qryData);
    }
  }

  if (args.outputType == 0) {
    countMatches(aligner.matchCounts[chunkQryNum], qryData);
  } else {
    size_t oldNumOfAlns = aligner.textAlns.size();
    scan(aligner, qrySeqs, qryData, matrices);
    cullFinalAlignments(aligner.textAlns, oldNumOfAlns, finalCullingLimit);
  }

  qryData.seq = qrySeqs.seqWriter() + qrySeqs.padBeg(qryData.seqNum);

  if (args.tantanSetting && !args.isKeepLowercase) {
    queryAlph.makeUppercase(qryData.seq + qryData.seqBeg,
			    qryData.seq + qryData.seqEnd);
  }
}

static void splitAlignments(LastSplitter &splitter,
			    std::vector<AlignmentText> &textAlns,
			    bool isQryQual) {
  if (!args.isSplit) return;

  unsigned linesPerMaf =
    3 + isQryQual + (args.outputType > 3) + (args.outputType > 6);

  splitter.reserve(textAlns.size());
  std::vector<char *> linePtrs((linesPerMaf + 1) * textAlns.size());

  char **beg = linePtrs.empty() ? 0 : &linePtrs[0];
  for (size_t i = 0; i < textAlns.size(); ++i) {
    char **end = beg;
    char *text = textAlns[i].text;
    for (unsigned j = 0; j < linesPerMaf; ++j) {
      *end++ = text;
      text = strchr(text, '\n');
      *text++ = 0;
    }
    *end = text;
    splitter.addMaf(beg, end, false);
    beg = end + 1;
  }

  splitter.split(args.splitOpts, splitParams, false);
  clearAlignments(textAlns);
}

static void alignOneQuery(LastAligner &aligner, MultiSequence &qrySeqs,
			  size_t qryNum, size_t chunkQryNum,
			  size_t finalCullingLimit, bool isFirstVolume) {
  size_t padBeg = qrySeqs.padBeg(qryNum);
  size_t padEnd = qrySeqs.padEnd(qryNum);
  size_t padLen = padEnd - padBeg;

  const uchar *qual = qrySeqs.qualityReader();
  if (qual) qual += padBeg * qrySeqs.qualsPerLetter();

  int *qualityPssm = qualityPssmSpace(aligner, padLen);

  SeqData qryData = {qryNum,
    padLen,
    qrySeqs.seqBeg(qryNum) - padBeg,
    qrySeqs.seqEnd(qryNum) - padBeg,
    args.isFrameshift() ? (padLen / 3) : 0,
    qrySeqs.seqWriter() + padBeg,
    qrySeqs.seqReader() + padBeg,
    qrySeqs.seqReader() + padEnd,
    qual,
    qualityPssm,
    getQueryPssm(qualityPssm, qrySeqs, padBeg)};

  if (isFirstVolume) {
    aligner.numOfNormalLetters +=
      queryAlph.countNormalLetters(qryData.seqPadBeg + qryData.seqBeg,
				   qryData.seqPadBeg + qryData.seqEnd);
    aligner.numOfSequences += 1;
  }

  int qryStrand = args.strand;
  if (qryStrand == 3) qryStrand = 1 - qryNum;
  if (qryStrand == 4) qryStrand = qryNum;

  std::vector<AlignmentText> &textAlns = aligner.textAlns;
  size_t oldNumOfAlns = textAlns.size();

  if (qryStrand == 2 && !isFirstVolume)
    qrySeqs.reverseComplementOneSequence(qryNum, queryAlph.complement);

  if (qryStrand != 0)
    translateAndScan(aligner, qrySeqs, qryData, chunkQryNum, finalCullingLimit,
		     fwdMatrices);

  if (qryStrand == 2 || (qryStrand == 0 && isFirstVolume))
    qrySeqs.reverseComplementOneSequence(qryNum, queryAlph.complement);

  if (qryStrand != 1)
    translateAndScan(aligner, qrySeqs, qryData, chunkQryNum, finalCullingLimit,
		     args.isQueryStrandMatrix ? revMatrices : fwdMatrices);

  if (numOfVolumes < 2) {
    if (isCollatedAlignments()) {
      sort(textAlns.begin() + oldNumOfAlns, textAlns.end());
    }
    splitAlignments(aligner.splitter, textAlns, qrySeqs.qualsPerLetter());
  }
}

static size_t alignSomeQueries(size_t chunkNum, unsigned volume) {
  size_t numOfChunks = aligners.size();
  LastAligner &aligner = aligners[chunkNum];
  size_t beg = firstSequenceInChunk(qrySeqsGlobal, numOfChunks, chunkNum);
  size_t end = firstSequenceInChunk(qrySeqsGlobal, numOfChunks, chunkNum + 1);
  bool isMultiVolume = (numOfVolumes > 1);
  bool isFirstVolume = (volume == 0);
  size_t finalCullingLimit = args.cullingLimitForFinalAlignments ?
    args.cullingLimitForFinalAlignments : isMultiVolume;
  if (args.outputType == 0 && isFirstVolume) {
    aligner.matchCounts.resize(end - beg);
  }
  for (size_t i = beg; i < end; ++i) {
    alignOneQuery(aligner, qrySeqsGlobal, i, i - beg,
		  finalCullingLimit, isFirstVolume);
  }
  if (isMultiVolume && volume + 1 == numOfVolumes) {
    std::vector<AlignmentText> &textAlns = aligner.textAlns;
    cullFinalAlignments(textAlns, 0, args.cullingLimitForFinalAlignments);
    sort(textAlns.begin(), textAlns.end());
    splitAlignments(aligner.splitter, textAlns,
		    qrySeqsGlobal.qualsPerLetter());
  }
  return beg;
}

static void scanOneVolume(unsigned volume, unsigned numOfThreadsLeft) {
  size_t firstSequence = 0;
  if (numOfThreadsLeft > 1) {
#ifdef HAS_CXX_THREADS
    std::thread t(scanOneVolume, volume, numOfThreadsLeft - 1);
    // Exceptions from threads are not handled nicely, but I don't
    // think it matters much.
    firstSequence = alignSomeQueries(numOfThreadsLeft - 1, volume);
    t.join();
#endif
  } else {
    firstSequence = alignSomeQueries(0, volume);
  }
  if (volume + 1 == numOfVolumes) {
    LastAligner &aligner = aligners[numOfThreadsLeft - 1];
    writeCounts(aligner.matchCounts, qrySeqsGlobal, firstSequence);
    aligner.matchCounts.clear();
    printAlignments(aligner.textAlns);
    clearAlignments(aligner.textAlns);
    if (!aligner.splitter.isOutputEmpty()) {
      aligner.splitter.printOutput();
      aligner.splitter.clearOutput();
    }
  }
}

static void openIfFile(mcf::izstream &z, const char *fileName) {
  if (fileName && !isSingleDash(fileName)) openOrThrow(z, fileName);
}

static std::istream &getInput(mcf::izstream &f) {
  return f.is_open() ? f : std::cin;
}

static bool readPairedSequences(MultiSequence &qrySeqs) {
  const bool isMask = (args.maskLowercase > 1);
  std::istream &in1 = getInput(querySequenceFile);
  std::istream &in2 =
    querySequenceFileNames[1] ? getInput(querySequenceFile2) : in1;
  bool x = false;
  bool y = false;
  {
    std::lock_guard<std::mutex> lockGuard(inputMutex);
    if (appendSequence(qrySeqs, in1, -1, args.inputFormat, queryAlph, isMask))
      x = true;
    if (appendSequence(qrySeqs, in2, -1, args.inputFormat, queryAlph, isMask))
      y = true;
  }
  if (x != y) err("unequal numbers of paired sequences");
  if (x) qrySeqs.fixPairedSequenceNames();
  return x;
}

static bool readSequenceData(MultiSequence &qrySeqs) {
  if (args.isPairedQuerySequences) return readPairedSequences(qrySeqs);
  const size_t maxPairedSeqLen = 2000;  // xxx ???
  const bool isMask = (args.maskLowercase > 1);
  std::lock_guard<std::mutex> lockGuard(inputMutex);

  while (*querySequenceFileNames) {
    bool isFile = !isSingleDash(*querySequenceFileNames);
    std::istream &in = isFile ? querySequenceFile : std::cin;
    if (appendSequence(qrySeqs, in, -1, args.inputFormat, queryAlph, isMask)) {
      if (qrySeqs.isFinished() && qrySeqs.seqLen(0) <= maxPairedSeqLen) {
	appendSequence(qrySeqs, in, -1, args.inputFormat, queryAlph, isMask);
      }
      return true;
    }
    if (isFile) querySequenceFile.close();
    ++querySequenceFileNames;
    openIfFile(querySequenceFile, *querySequenceFileNames);
  }

  return false;
}

static void runOneThread(unsigned threadNum) {
  LastAligner &aligner = aligners[threadNum];
  LastSplitter &splitter = aligner.splitter;
  std::vector< std::vector<countT> > &matchCounts = aligner.matchCounts;
  std::vector<AlignmentText> &textAlns = aligner.textAlns;
  MultiSequence qrySeqs;
  initSequences(qrySeqs, queryAlph, args.isTranslated(), false);

  while (readSequenceData(qrySeqs)) {
    if (!qrySeqs.isFinished()) throwSeqTooBig();
    encodeSequences(qrySeqs, args.inputFormat, queryAlph,
		    args.isKeepLowercase, 0);
    if (args.outputType == 0) matchCounts.resize(qrySeqs.finishedSequences());
    for (size_t i = 0; i < qrySeqs.finishedSequences(); ++i) {
      alignOneQuery(aligner, qrySeqs, i, i,
		    args.cullingLimitForFinalAlignments, true);
    }
    if (!splitter.isOutputEmpty()) {
      {
	std::lock_guard<std::mutex> lockGuard(outputMutex);
	splitter.printOutput();
      }
      splitter.clearOutput();
    } else if (!textAlns.empty()) {
      {
	std::lock_guard<std::mutex> lockGuard(outputMutex);
	printAlignments(textAlns);
      }
      clearAlignments(textAlns);
    } else if (!matchCounts.empty()) {
      {
	std::lock_guard<std::mutex> lockGuard(outputMutex);
	writeCounts(matchCounts, qrySeqs, 0);
      }
      matchCounts.clear();
    }
    qrySeqs.reinitForAppending();
  }
}

static void runOneThreadSafely(unsigned threadNum) {
  try {
    runOneThread(threadNum);
  } catch (const std::bad_alloc &e) {
    std::cerr << args.programName << ": out of memory\n";
    raise(SIGTERM);
  } catch (const std::exception &e) {
    std::cerr << args.programName << ": " << e.what() << '\n';
    raise(SIGTERM);  // quick_exit doesn't work on Mac
  }
}

static void runThreads(unsigned numOfThreads) {
  if (numOfThreads > 1) {
#ifdef HAS_CXX_THREADS
    std::thread t(runThreads, numOfThreads - 1);
    runOneThreadSafely(numOfThreads - 1);
    t.join();
#endif
  } else {
    runOneThreadSafely(0);
  }
}

void readIndex(const std::string &baseName, size_t seqCount, int bitsPerBase,
	       int bitsPerInt, bool isCaseSensitive) {
  LOG( "reading " << baseName << "..." );
  refSeqs.fromFiles(baseName, seqCount,
		    referenceFormat != sequenceFormat::fasta,
		    bitsPerBase == 4, bitsPerInt == 32);
  for( unsigned x = 0; x < numOfIndexes; ++x ){
    if( numOfIndexes > 1 ){
      suffixArrays[x].fromFiles(baseName + char('a' + x),
				bitsPerInt, isCaseSensitive,
				alph.encode, alph.letters);
    } else {
      suffixArrays[x].fromFiles(baseName, bitsPerInt, isCaseSensitive,
				alph.encode, alph.letters);
    }
  }

  const std::vector<CyclicSubsetSeed> &seeds = suffixArrays[0].getSeeds();
  assert(!seeds.empty());  // xxx what if numOfIndexes==0 ?
  makeWordsFinder(wordsFinder, &seeds[0], seeds.size(), alph.encode,
		  isCaseSensitive);

  if (scoreMatrix.isCodonCols()) {
    for (unsigned x = 0; x < numOfIndexes; ++x) {
      std::vector<CyclicSubsetSeed> &s = suffixArrays[x].getSeeds();
      for (size_t i = 0; i < s.size(); ++i) {
	s[i].compose(geneticCode.getCodonToAmino());
      }
    }
  }
}

int calcMinScoreGapless(double numLettersInReference) {
  if (args.minScoreGapless >= 0) return args.minScoreGapless;

  // ***** Default setting for minScoreGapless *****

  // This attempts to ensure that the gapped alignment phase will be
  // reasonably fast relative to the gapless alignment phase.

  if (!gaplessEvaluer.isGood()) return ceil(args.minScoreGapped);

  double n = args.maxGaplessAlignmentsPerQueryPosition;
  if (args.maxGaplessAlignmentsPerQueryPosition + 1 == 0) n = 10;  // ?

  // The number of gapless extensions per query position is
  // proportional to: maxGaplessAlignmentsPerQueryPosition * numOfIndexes.
  // We want the number of gapped extensions per query position to be
  // proportional to this:
  double propConst = 0.0004;  // xxx ???
  double e = n * numOfIndexes * propConst;
  // The proportionality constant was guesstimated by some limited
  // trial-and-error.  It should depend on the relative speeds of
  // gapless and gapped extensions.

  double s = gaplessEvaluer.minScore(e, numLettersInReference);
  s = std::max(1.0, s);
  return ceil(std::min(s, args.minScoreGapped));
}

// Read one database volume
void readVolume(unsigned volumeNumber, int bitsPerBase, int bitsPerInt,
		bool isCaseSensitive) {
  std::string baseName = args.lastdbName + stringify(volumeNumber);
  size_t seqCount = -1;
  size_t seqLen = -1;
  readInnerPrj(baseName + ".prj", seqCount, seqLen);
  minScoreGapless = calcMinScoreGapless(seqLen);
  readIndex(baseName, seqCount, bitsPerBase, bitsPerInt, isCaseSensitive);
}

// Scan one batch of query sequences against all database volumes
void scanAllVolumes(int bitsPerBase, int bitsPerInt, bool isCaseSensitive) {
  encodeSequences(qrySeqsGlobal, args.inputFormat, queryAlph,
		  args.isKeepLowercase, 0);

  for (unsigned i = 0; i < numOfVolumes; ++i) {
    if (refSeqs.unfinishedSize() == 0 || numOfVolumes > 1) {
      readVolume(i, bitsPerBase, bitsPerInt, isCaseSensitive);
    }
    scanOneVolume(i, aligners.size());
  }
}

void writeHeader(countT numOfRefSeqs, countT refLetters, std::ostream &out) {
  out << "# LAST version " <<
#include "version.hh"
      << "\n";
  out << "#\n";
  args.writeCommented( out );
  out << "# Reference sequences=" << numOfRefSeqs
      << " normal letters=" << refLetters << "\n";
  if( args.outputType > 0 ) evaluer.writeCommented( out );
  out << "#\n";

  if( args.outputType == 0 ){  // we just want hit counts
    out << "# length\tcount\n"
	<< "#\n";
  } else {  // we want alignments
    if( args.inputFormat != sequenceFormat::pssm || !args.matrixFile.empty() ){
      // we're not reading PSSMs, or we bothered to specify a matrix file
      scoreMatrix.writeCommented( out );
      // Write lambda?
      out << "#\n";
    }

    if( args.outputFormat != 'b' && args.outputFormat != 'B' ) {
      out << "# Coordinates are 0-based.  For - strand matches, coordinates\n";
      out << "# in the reverse complement of the 2nd sequence are used.\n";
      out << "#\n";
    }

    if( args.outputFormat == 't' ){
      out << "# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\t"
	  << "name2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n";
    }
    if( args.outputFormat == 'm' ){
      out << "# name start alnSize strand seqSize alignment\n"
	  << "#\n";
    }
    if( args.outputFormat == 'b' || args.outputFormat == 'B' ){
      out << "# Fields: query id, subject id, % identity, alignment length, "
	  << "mismatches, gap opens, q. start, q. end, s. start, s. end";
      if( evaluer.isGood() ) out << ", evalue, bit score";
      if( args.outputFormat == 'B' )
	out << ", query length, subject length, raw score";
      out << '\n';
    }

    if (args.isSplit) {
      args.splitOpts.print();
      splitParams.print();
      std::cout << "#\n";
    }
  }
}

void lastal(int argc, char **argv) {
  args.fromArgs(argc, argv);
  args.resetCumulativeOptions();  // because we will do fromArgs again
  const LastdbData prj = readOuterPrj(args.lastdbName + ".prj");
  const bool isDna = (alph.letters == alph.dna);
  const bool isProtein = alph.isProtein();
  args.fromArgs(argc, argv);  // command line overrides prj file

  std::string matrixName = args.matrixName( isProtein );
  std::string matrixFile;
  if( !matrixName.empty() ){
    matrixFile = ScoreMatrix::stringFromName( matrixName );
    args.resetCumulativeOptions();
    args.fromString( matrixFile );  // read options from the matrix file
    sequenceFormat::Enum f = args.inputFormat;
    args.fromArgs( argc, argv );  // command line overrides matrix file
    if (isUseQuality(f) != isUseQuality(args.inputFormat)) {
      ERR("option -Q is inconsistent with the matrix file");
    }
  }

  if (prj.minSeedLimit > 1) {
    if (args.outputType == 0)
      ERR("can't use option -j 0: need to re-run lastdb with i <= 1");
    if (prj.minSeedLimit > args.oneHitMultiplicity)
      ERR("can't use option -m < " + stringify(prj.minSeedLimit) +
	  ": need to re-run lastdb with i <= " +
	  stringify(args.oneHitMultiplicity));
    if (args.minHitDepth > 1)
      ERR("can't use option -l > 1: need to re-run lastdb with i <= 1");
  }

  aligners.resize( decideNumberOfThreads( args.numOfThreads,
					  args.programName, args.verbosity ) );
  bool isMultiVolume = (numOfVolumes + 1 > 0 && numOfVolumes > 1);
  args.setDefaultsFromAlphabet(isDna, isProtein, prj.strand,
			       prj.isKeepLowercase, prj.tantanSetting,
			       prj.isCaseSensitive, numOfVolumes + 1 > 0,
			       prj.minimizerWindow);
  makeScoreMatrix(matrixName, matrixFile);
  gapCosts.assign(args.delOpenCosts, args.delGrowCosts,
		  args.insOpenCosts, args.insGrowCosts,
		  args.frameshiftCosts, args.gapPairCost,
		  fwdMatrices.stats.lambda());

  if( args.isTranslated() ){
    if( isDna )  // allow user-defined alphabet
      ERR( "expected protein database, but got DNA" );
    queryAlph.init(queryAlph.dna, false);
    geneticCode.fromString(GeneticCode::stringFromName(args.geneticCodeFile));
    if (scoreMatrix.isCodonCols()) {
      geneticCode.initCodons(queryAlph.encode, alph.encode,
			     args.maskLowercase > 0, args.maskLowercase < 3);
    } else {
      geneticCode.codeTableSet( alph, queryAlph );
    }
  } else {
    queryAlph = alph;
  }

  if (prj.bitsPerBase < CHAR_BIT) {
    if (args.isGreedy) err("can't use option -M with 4-bit lastdb");
    if (isUseFastq(referenceFormat) && isUseFastq(args.inputFormat))
      err("can't do fastq-versus-fastq with 4-bit lastdb");
  }

  if (args.tantanSetting) {
    int t = args.tantanSetting;
    int m = args.maxRepeatUnit;
    if (scoreMatrix.isCodonCols()) {
      if (m < 0) m = 100;
      tantanMasker.init(0, t==2, t==3, m, queryAlph.letters, queryAlph.encode);
    } else {
      if (m < 0) m = prj.maxRepeatUnit ? prj.maxRepeatUnit
		   : isProtein ? 50 : 100;
      tantanMasker.init(isProtein, t==2, t==3, m, alph.letters, alph.encode);
    }
  }

  char defaultInputName[] = "-";
  char* defaultInput[] = { defaultInputName, 0 };
  char** inputBegin = argv + args.inputStart;
  querySequenceFileNames = *inputBegin ? inputBegin : defaultInput;

  if (args.outputType > 0) {
    calculateScoreStatistics(matrixName, prj.numOfLetters, prj.maxSeqLen);
  }

  double minScore = -1;
  double eg2 = -1;
  if (evaluer.isGood()) {
    minScore = (args.maxEvalue > 0) ? evaluer.minScore(args.maxEvalue, 1e18)
      : evaluer.minScore(args.queryLettersPerRandomAlignment);
    if (args.scoreType == 0 || args.outputType < 2)
      minScore = ceil(std::max(1.0, minScore));
    eg2 = 1e18 * evaluer.evaluePerArea(minScore);
  }
  args.setDefaultsFromMatrix(fwdMatrices.stats.lambda(), minScore, eg2);

  minScoreGapless = calcMinScoreGapless(prj.numOfLetters);
  if (!isMultiVolume) args.minScoreGapless = minScoreGapless;
  if (args.outputType > 0) makeQualityScorers();

  if (args.isSplit) {
    setLastSplitParams(splitParams, args.splitOpts, scoreMatrix.cells,
		       scoreMatrix.rowSymbols.c_str(),
		       scoreMatrix.colSymbols.c_str(),
		       args.isQueryStrandMatrix,
		       args.delOpenCosts[0], args.delGrowCosts[0],
		       args.insOpenCosts[0], args.insGrowCosts[0],
		       args.temperature, prj.numOfLetters, args.inputFormat);
  }

  if (numOfVolumes + 1 == 0) {
    readIndex(args.lastdbName, prj.numOfSeqs, prj.bitsPerBase, prj.bitsPerInt,
	      prj.isCaseSensitive);
    numOfVolumes = 1;
  }

  writeHeader(prj.numOfSeqs, prj.numOfLetters, std::cout);
  countT queryBatchCount = 0;

  if (args.batchSize < 1) {
    openIfFile(querySequenceFile, *querySequenceFileNames);
    if (args.isPairedQuerySequences && querySequenceFileNames[1]) {
      if (querySequenceFileNames[2]) err("can't use >2 files with option -2");
      openIfFile(querySequenceFile2, querySequenceFileNames[1]);
    }
    runThreads(aligners.size());
  } else {
    if (args.isPairedQuerySequences) err("can't use option -2 with batches");
    size_t maxSeqLen = -1;
    initSequences(qrySeqsGlobal, queryAlph, args.isTranslated(), false);
    for (char **i = querySequenceFileNames; *i; ++i) {
      mcf::izstream inFileStream;
      std::istream& in = openIn(*i, inFileStream);
      LOG("reading " << *i << "...");
      while (appendSequence(qrySeqsGlobal, in, maxSeqLen, args.inputFormat,
			    queryAlph, args.maskLowercase > 1)) {
	if (qrySeqsGlobal.isFinished()) {
	  maxSeqLen = args.batchSize;
	} else {
	  if (qrySeqsGlobal.finishedSequences() == 0) throwSeqTooBig();
	  // this enables downstream parsers to read one batch at a time:
	  std::cout << "# batch " << queryBatchCount++ << "\n";
	  scanAllVolumes(prj.bitsPerBase, prj.bitsPerInt, prj.isCaseSensitive);
	  qrySeqsGlobal.reinitForAppending();
	  maxSeqLen = -1;
	}
      }
    }
    if (qrySeqsGlobal.finishedSequences() > 0) {
      std::cout << "# batch " << queryBatchCount << "\n";
      scanAllVolumes(prj.bitsPerBase, prj.bitsPerInt, prj.isCaseSensitive);
    }
  }

  countT numOfSequences = 0;
  countT numOfNormalLetters = 0;
  for (size_t i = 0; i < aligners.size(); ++i) {
    numOfSequences += aligners[i].numOfSequences;
    numOfNormalLetters += aligners[i].numOfNormalLetters;
  }
  std::cout << "# Query sequences=" << numOfSequences
	    << " normal letters=" << numOfNormalLetters << "\n";
}

int main( int argc, char** argv )
try{
  lastal( argc, argv );
  if (!flush(std::cout)) ERR( "write error" );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << argv[0] << ": out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << argv[0] << ": " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
