// Copyright 2008, 2009, 2010, 2011, 2014 Martin C. Frith

#include "ScoreMatrix.hh"
#include "ScoreMatrixData.hh"
#include "zio.hh"
#include <sstream>
#include <iomanip>
#include <algorithm>  // min, max
#include <stdexcept>
#include <cassert>
#include <cctype>  // toupper, tolower
#include <stddef.h>  // size_t
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

namespace cbrc{

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
  if( !iss ) ERR( "can't read the score matrix" );
}

void ScoreMatrix::init(const uchar symbolToIndex[]) {
  assert( !rowSymbols.empty() && !colSymbols.empty() );

  for (std::string::iterator i = rowSymbols.begin(); i < rowSymbols.end(); ++i)
    *i = std::toupper( *i );

  for (std::string::iterator i = colSymbols.begin(); i < colSymbols.end(); ++i)
    *i = std::toupper( *i );

  minScore = cells[0][0];
  maxScore = cells[0][0];

  for( size_t i = 0; i < rowSymbols.size(); ++i ){
    for( size_t j = 0; j < colSymbols.size(); ++j ){
      minScore = std::min( minScore, cells[i][j] );
      maxScore = std::max( maxScore, cells[i][j] );
    }
  }

  // set default score = minScore:
  for( unsigned i = 0; i < MAT; ++i ){
    for( unsigned j = 0; j < MAT; ++j ){
      caseSensitive[i][j] = minScore;
      caseInsensitive[i][j] = minScore;
    }
  }

  for( size_t i = 0; i < rowSymbols.size(); ++i ){
    for( size_t j = 0; j < colSymbols.size(); ++j ){
      uchar iu = symbolToIndex[ uchar(rowSymbols[i]) ];
      uchar ju = symbolToIndex[ uchar(colSymbols[j]) ];
      uchar il = symbolToIndex[ std::tolower( rowSymbols[i] ) ];
      uchar jl = symbolToIndex[ std::tolower( colSymbols[j] ) ];
      if( il >= MAT )
        ERR( std::string("bad letter in score matrix: ") + rowSymbols[i] );
      if( jl >= MAT )
        ERR( std::string("bad letter in score matrix: ") + colSymbols[j] );
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
  assert( z < MAT );
  for( unsigned i = 0; i < MAT; ++i ){
    caseSensitive[z][i] = -INF;
    caseSensitive[i][z] = -INF;
    caseInsensitive[z][i] = -INF;
    caseInsensitive[i][z] = -INF;    
  }
}

void ScoreMatrix::writeCommented( std::ostream& stream ) const{
  stream << "# " << ' ';
  for( size_t i = 0; i < colSymbols.size(); ++i ){
    stream << ' ' << std::setw(OUTPAD) << colSymbols[i];
  }
  stream << '\n';

  for( size_t i = 0; i < rowSymbols.size(); ++i ){
    stream << "# " << rowSymbols[i];
    for( size_t j = 0; j < colSymbols.size(); ++j ){
      stream << ' ' << std::setw(OUTPAD) << cells[i][j];
    }
    stream << '\n';
  }
}

std::istream& operator>>( std::istream& stream, ScoreMatrix& m ){
  std::string tmpRowSymbols;
  std::string tmpColSymbols;
  std::vector< std::vector<int> > tmpCells;
  std::string line;
  int state = 0;

  while( std::getline( stream, line ) ){
    std::istringstream iss(line);
    char c;
    if( !(iss >> c) ) continue;  // skip blank lines
    if( state == 0 ){
      if( c == '#' ) continue;  // skip comment lines at the top
      do{
	tmpColSymbols.push_back(c);
      }while( iss >> c );
      state = 1;
    }
    else{
      tmpRowSymbols.push_back(c);
      tmpCells.resize( tmpCells.size() + 1 );
      int score;
      while( iss >> score ){
	tmpCells.back().push_back(score);
      }
      if (tmpCells.back().size() != tmpColSymbols.size()) {
	ERR("bad score matrix");
      }
    }
  }

  if( stream.eof() && !stream.bad() && !tmpRowSymbols.empty() ){
    stream.clear();
    m.rowSymbols.swap(tmpRowSymbols);
    m.colSymbols.swap(tmpColSymbols);
    m.cells.swap(tmpCells);
  }

  return stream;
}

}
