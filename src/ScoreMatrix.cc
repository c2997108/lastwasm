// Copyright 2008, 2009, 2010, 2011, 2014 Martin C. Frith

#include "ScoreMatrix.hh"
#include "ScoreMatrixData.hh"
#include "qualityScoreUtil.hh"
#include "zio.hh"
#include <sstream>
#include <iomanip>
#include <algorithm>  // min, max
#include <stdexcept>
#include <cassert>
#include <cctype>  // toupper, tolower
#include <stddef.h>  // size_t

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

static void makeUppercase(std::string& s) {
  for (size_t i = 0; i < s.size(); ++i) {
    unsigned char c = s[i];
    s[i] = std::toupper(c);
  }
}

namespace cbrc{

typedef std::runtime_error Err;

const char *ScoreMatrix::canonicalName( const std::string& name ){
  for( size_t i = 0; i < COUNTOF(scoreMatrixNicknames); ++i )
    if( name == scoreMatrixNicknames[i].nickname )
      return scoreMatrixNicknames[i].realname;
  return name.c_str();
}

std::string ScoreMatrix::stringFromName( const std::string& name ){
  std::string n = canonicalName( name );

  for( size_t i = 0; i < COUNTOF(scoreMatrices); ++i )
    if( n == scoreMatrices[i].name )
      return scoreMatrices[i].text;

  return slurp( n.c_str() );
}

void ScoreMatrix::setMatchMismatch(int matchScore, int mismatchCost,
				   const std::string& symbols) {
  rowSymbols.assign(symbols.begin(), symbols.end());
  colSymbols.assign(symbols.begin(), symbols.end());

  size_t size = symbols.size();
  cells.resize(size);

  for (size_t i = 0; i < size; ++i) {
    cells[i].assign(size, -mismatchCost);
    cells[i][i] = matchScore;
  }
}

void ScoreMatrix::fromString( const std::string& matString ){
  std::istringstream iss(matString);
  iss >> *this;
  if (!iss) throw Err("can't read the score matrix");
}

static unsigned s2i(const uchar symbolToIndex[], uchar c) {
  return symbolToIndex[c];
}

static void upperAndLowerIndex(unsigned tooBig, const uchar symbolToIndex[],
			       char symbol, unsigned& upper, unsigned& lower) {
  uchar s = symbol;
  upper = symbolToIndex[s];
  lower = symbolToIndex[std::tolower(s)];
  if (upper >= tooBig || lower >= tooBig) {
    throw Err(std::string("bad letter in score matrix: ") + symbol);
  }
}

void ScoreMatrix::init(const uchar symbolToIndex[]) {
  unsigned fastMatrixSize = ALPHABET_CAPACITY;
  assert(!rowSymbols.empty());
  assert(!colSymbols.empty());

  makeUppercase(rowSymbols);
  makeUppercase(colSymbols);

  minScore = cells[0][0];
  maxScore = cells[0][0];

  for( size_t i = 0; i < rowSymbols.size(); ++i ){
    for( size_t j = 0; j < colSymbols.size(); ++j ){
      minScore = std::min( minScore, cells[i][j] );
      maxScore = std::max( maxScore, cells[i][j] );
    }
  }

  // set default score = minScore:
  for (unsigned i = 0; i < fastMatrixSize; ++i) {
    for (unsigned j = 0; j < fastMatrixSize; ++j) {
      caseSensitive[i][j] = minScore;
      caseInsensitive[i][j] = minScore;
    }
  }

  for( size_t i = 0; i < rowSymbols.size(); ++i ){
    for( size_t j = 0; j < colSymbols.size(); ++j ){
      unsigned iu, il, ju, jl;
      upperAndLowerIndex(fastMatrixSize, symbolToIndex, rowSymbols[i], iu, il);
      upperAndLowerIndex(fastMatrixSize, symbolToIndex, colSymbols[j], ju, jl);
      caseSensitive[iu][jl] = std::min( cells[i][j], 0 );
      caseSensitive[il][ju] = std::min( cells[i][j], 0 );
      caseSensitive[il][jl] = std::min( cells[i][j], 0 );
      caseSensitive[iu][ju] = cells[i][j];  // careful: maybe il==iu or jl==ju
      caseInsensitive[iu][ju] = cells[i][j];
      caseInsensitive[iu][jl] = cells[i][j];
      caseInsensitive[il][ju] = cells[i][j];
      caseInsensitive[il][jl] = cells[i][j];
    }
  }

  // set a hugely negative score for the delimiter symbol:
  uchar delimiter = ' ';
  uchar z = symbolToIndex[delimiter];
  assert(z < fastMatrixSize);
  for (unsigned i = 0; i < fastMatrixSize; ++i) {
    caseSensitive[z][i] = -INF;
    caseSensitive[i][z] = -INF;
    caseInsensitive[z][i] = -INF;
    caseInsensitive[i][z] = -INF;    
  }
}

static void calcSomeLetterProbs(const std::string &symbols,
				const std::vector<double> &freqs,
				double probs[], unsigned alphabetSizeForProbs,
				const uchar symbolToIndex[]) {
  double sum = 0;
  for (size_t i = 0; i < freqs.size(); ++i) {
    unsigned j = s2i(symbolToIndex, symbols[i]);
    if (j < alphabetSizeForProbs) {
      if (freqs[i] < 0) throw Err("bad score matrix: letter frequency < 0");
      sum += freqs[i];
    }
  }
  if (sum <= 0) throw Err("bad score matrix: no positive letter frequencies");

  std::fill_n(probs, alphabetSizeForProbs, 0.0);

  for (size_t i = 0; i < freqs.size(); ++i) {
    unsigned j = s2i(symbolToIndex, symbols[i]);
    if (j < alphabetSizeForProbs) {
      probs[j] = freqs[i] / sum;
    }
  }
}

void ScoreMatrix::calcLetterProbs(double rowProbs[], double colProbs[],
				  unsigned alphabetSizeForProbs,
				  const uchar symbolToIndex[]) const {
  calcSomeLetterProbs(rowSymbols, rowFrequencies, rowProbs,
		      alphabetSizeForProbs, symbolToIndex);
  calcSomeLetterProbs(colSymbols, colFrequencies, colProbs,
		      alphabetSizeForProbs, symbolToIndex);
}

void ScoreMatrix::writeCommented( std::ostream& stream ) const{
  int colWidth = colSymbols.size() < 20 ? 3 : 2;

  stream << "# " << ' ';
  for( size_t i = 0; i < colSymbols.size(); ++i ){
    stream << ' ' << std::setw(colWidth) << colSymbols[i];
  }
  stream << '\n';

  for( size_t i = 0; i < rowSymbols.size(); ++i ){
    stream << "# " << rowSymbols[i];
    for( size_t j = 0; j < colSymbols.size(); ++j ){
      stream << ' ' << std::setw(colWidth) << cells[i][j];
    }
    stream << '\n';
  }
}

std::istream& operator>>( std::istream& stream, ScoreMatrix& m ){
  std::string tmpRowSymbols;
  std::string tmpColSymbols;
  std::vector< std::vector<int> > tmpCells;
  std::vector<double> tmpRowFreqs;
  std::vector<double> tmpColFreqs;
  std::string line, word;

  while (stream) {
    if (!getline(stream, line)) {
      if (stream.eof() && !stream.bad()) stream.clear(std::ios::eofbit);
      break;
    }
    std::istringstream iss(line);
    if (!(iss >> word)) continue;  // skip blank lines
    if (tmpColSymbols.empty()) {
      if (word[0] == '#') continue;  // skip comment lines at the top
      do {
	if (word.size() > 1) stream.setstate(std::ios::failbit);
	tmpColSymbols.push_back(word[0]);
      } while (iss >> word);
    } else {
      std::vector<int> row;
      int score;
      double freq;
      for (size_t i = 0; i < tmpColSymbols.size(); ++i) {
	iss >> score;
	row.push_back(score);
      }
      if (word.size() == 1 && iss) {
	tmpRowSymbols.push_back(word[0]);
	tmpCells.push_back(row);
	if (iss >> freq) {
	  tmpRowFreqs.push_back(freq);
	  if (tmpRowFreqs.size() < tmpCells.size()) {
	    stream.setstate(std::ios::failbit);
	  }
	}
      } else {
	std::istringstream iss2(line);
	while (iss2 >> freq) {
	  tmpColFreqs.push_back(freq);
	}
	if (tmpColFreqs.size() > tmpColSymbols.size() || tmpColFreqs.empty()) {
	  stream.setstate(std::ios::failbit);
	}
	break;
      }
    }
  }

  if (tmpCells.empty() || tmpRowFreqs.empty() != tmpColFreqs.empty()) {
    stream.setstate(std::ios::failbit);
  }

  if (stream) {
    m.rowSymbols.swap(tmpRowSymbols);
    m.colSymbols.swap(tmpColSymbols);
    m.cells.swap(tmpCells);
    m.rowFrequencies.swap(tmpRowFreqs);
    m.colFrequencies.swap(tmpColFreqs);
  }

  return stream;
}

const char *ntAmbiguities[] = {
  "M" "AC",
  "S" "CG",
  "K" "GT",
  "W" "TA",
  "R" "AG",
  "Y" "CT",
  "B" "CGT",
  "D" "AGT",
  "H" "ACT",
  "V" "ACG",
  "N" "ACGT"
};

const char *aaAmbiguities[] = {
  "X" "ACDEFGHIKLMNPQRSTVWY"
};

static bool isIn(const std::string& s, char x) {
  return find(s.begin(), s.end(), x) != s.end();
}

static const char *ambiguityList(const char *ambiguities[],
				 size_t numOfAmbiguousSymbols,
				 char symbol, char scratch[]) {
  for (size_t i = 0; i < numOfAmbiguousSymbols; ++i) {
    if (ambiguities[i][0] == symbol) return ambiguities[i] + 1;
  }
  scratch[0] = symbol;
  return scratch;
}

static double symbolProbSum(const uchar symbolToIndex[],
			    const char *symbols, const double probs[]) {
  if (symbols[1] == 0) return 1;
  double p = 0;
  for (size_t i = 0; symbols[i]; ++i) {
    unsigned x = s2i(symbolToIndex, symbols[i]);
    p += probs[x];
  }
  return p;
}

static int jointScore(const uchar symbolToIndex[], int **fastMatrix,
		      double scale,
		      const double rowSymbolProbs[],
		      const double colSymbolProbs[],
		      const char *rSymbols, const char *cSymbols) {
  bool isOneRowSymbol = (rSymbols[1] == 0);
  bool isOneColSymbol = (cSymbols[1] == 0);

  double p = 0;
  for (size_t i = 0; rSymbols[i]; ++i) {
    for (size_t j = 0; cSymbols[j]; ++j) {
      unsigned x = s2i(symbolToIndex, rSymbols[i]);
      unsigned y = s2i(symbolToIndex, cSymbols[j]);
      double r = isOneRowSymbol ? 1 : rowSymbolProbs[x];
      double c = isOneColSymbol ? 1 : colSymbolProbs[y];
      p += r * c * probFromScore(scale, fastMatrix[x][y]);
    }
  }

  double rowProbSum = symbolProbSum(symbolToIndex, rSymbols, rowSymbolProbs);
  double colProbSum = symbolProbSum(symbolToIndex, cSymbols, colSymbolProbs);

  return scoreFromProb(scale, p / (rowProbSum * colProbSum));
}

void ScoreMatrix::addAmbiguousScores(bool isDna, bool isFullyAmbiguousRow,
				     bool isFullyAmbiguousCol,
				     const uchar symbolToIndex[],
				     double scale,
				     const double rowSymbolProbs[],
				     const double colSymbolProbs[]) {
  int *fastMatrix[ALPHABET_CAPACITY];
  std::copy(caseInsensitive, caseInsensitive + ALPHABET_CAPACITY, fastMatrix);

  char scratch[2] = {0};

  const char **ambiguities = isDna ? ntAmbiguities : aaAmbiguities;
  size_t numOfAmbig = isDna ? COUNTOF(ntAmbiguities) : COUNTOF(aaAmbiguities);
  size_t numOfAmbiguousRows = numOfAmbig - 1 + isFullyAmbiguousRow;
  size_t numOfAmbiguousCols = numOfAmbig - 1 + isFullyAmbiguousCol;

  for (size_t k = 0; k < numOfAmbiguousCols; ++k) {
    char ambiguousSymbol = ambiguities[k][0];
    if (isIn(colSymbols, ambiguousSymbol)) continue;
    colSymbols.push_back(ambiguousSymbol);
    for (size_t i = 0; i < rowSymbols.size(); ++i) {
      const char *rSymbols = ambiguityList(ambiguities, numOfAmbiguousRows,
					   rowSymbols[i], scratch);
      const char *cSymbols = ambiguities[k] + 1;
      int s = jointScore(symbolToIndex, fastMatrix, scale,
			 rowSymbolProbs, colSymbolProbs, rSymbols, cSymbols);
      cells[i].push_back(s);
    }
  }

  for (size_t k = 0; k < numOfAmbiguousRows; ++k) {
    char ambiguousSymbol = ambiguities[k][0];
    if (isIn(rowSymbols, ambiguousSymbol)) continue;
    rowSymbols.push_back(ambiguousSymbol);
    cells.resize(cells.size() + 1);
    for (size_t j = 0; j < colSymbols.size(); ++j) {
      const char *rSymbols = ambiguities[k] + 1;
      const char *cSymbols = ambiguityList(ambiguities, numOfAmbiguousCols,
					   colSymbols[j], scratch);
      int s = jointScore(symbolToIndex, fastMatrix, scale,
			 rowSymbolProbs, colSymbolProbs, rSymbols, cSymbols);
      cells.back().push_back(s);
    }
  }
}

}
