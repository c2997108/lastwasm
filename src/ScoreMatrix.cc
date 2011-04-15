// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

#include "ScoreMatrix.hh"
#include <sstream>
#include <iomanip>
#include <algorithm>  // min, max
#include <stdexcept>
#include <cassert>
#include <cctype>  // toupper, tolower
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

namespace cbrc{

// copied from NCBI BLAST:
const char* ScoreMatrix::blosum62 = "\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -1 -4\n\
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1 -2  0 -1 -4\n\
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  4 -3  0 -1 -4\n\
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4 -3  1 -1 -4\n\
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -1 -3 -1 -4\n\
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0 -2  4 -1 -4\n\
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1 -3  4 -1 -4\n\
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -4\n\
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0 -3  0 -1 -4\n\
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3  3 -3 -1 -4\n\
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4  3 -3 -1 -4\n\
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0 -3  1 -1 -4\n\
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3  2 -1 -1 -4\n\
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3  0 -3 -1 -4\n\
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -3 -1 -1 -4\n\
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0 -2  0 -1 -4\n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -1 -4\n\
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -2 -2 -1 -4\n\
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -1 -2 -1 -4\n\
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3  2 -2 -1 -4\n\
B -2 -1  4  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4 -3  0 -1 -4\n\
J -1 -2 -3 -3 -1 -2 -3 -4 -3  3  3 -3  2  0 -3 -2 -1 -2 -1  2 -3  3 -3 -1 -4\n\
Z -1  0  0  1 -3  4  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -2 -2 -2  0 -3  4 -1 -4\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\n\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n\
";

void ScoreMatrix::matchMismatch( int match, int mismatch,
				 const std::string& letters ){
  rows = letters;
  cols = letters;
  std::size_t size = letters.size();

  cells.resize( size );

  for( std::size_t i = 0; i < size; ++i ){
    cells[i].assign( size, -mismatch );
    cells[i][i] = match;
  }
}

void ScoreMatrix::fromString( const std::string& matString ){
  std::istringstream iss(matString);
  iss >> *this;
  assert(iss);
}

void ScoreMatrix::init( const uchar encode[] ){
  assert( !rows.empty() && !cols.empty() );

  for( std::string::iterator i = rows.begin(); i < rows.end(); ++i )
    *i = std::toupper( *i );

  for( std::string::iterator i = cols.begin(); i < cols.end(); ++i )
    *i = std::toupper( *i );

  minScore = cells[0][0];
  maxScore = cells[0][0];

  for( std::size_t i = 0; i < rows.size(); ++i ){
    for( std::size_t j = 0; j < cols.size(); ++j ){
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

  for( std::size_t i = 0; i < rows.size(); ++i ){
    for( std::size_t j = 0; j < cols.size(); ++j ){
      uchar x = encode[ uchar(rows[i]) ];
      uchar y = encode[ uchar(cols[j]) ];
      uchar a = encode[ std::tolower( rows[i] ) ];
      uchar b = encode[ std::tolower( cols[j] ) ];
      if( a >= MAT )
        ERR( std::string("bad letter in score matrix: ") + rows[i] );
      if( b >= MAT )
        ERR( std::string("bad letter in score matrix: ") + cols[j] );
      caseSensitive[x][b] = std::min( cells[i][j], 0 );
      caseSensitive[a][y] = std::min( cells[i][j], 0 );
      caseSensitive[a][b] = std::min( cells[i][j], 0 );
      caseSensitive[x][y] = cells[i][j];  // careful: maybe a==x or b==y
      caseInsensitive[x][y] = cells[i][j];
      caseInsensitive[x][b] = cells[i][j];
      caseInsensitive[a][y] = cells[i][j];
      caseInsensitive[a][b] = cells[i][j];
    }
  }

  // set a hugely negative score for the delimiter symbol:
  uchar z = encode[' '];
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
  for( std::size_t i = 0; i < cols.size(); ++i ){
    stream << ' ' << std::setw(OUTPAD) << cols[i];
  }
  stream << '\n';

  for( std::size_t i = 0; i < rows.size(); ++i ){
    stream << "# " << rows[i];
    for( std::size_t j = 0; j < cols.size(); ++j ){
      stream << ' ' << std::setw(OUTPAD) << cells[i][j];
    }
    stream << '\n';
  }
}

std::istream& operator>>( std::istream& stream, ScoreMatrix& m ){
  std::string tmpRows;
  std::string tmpCols;
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
	tmpCols.push_back(c);
      }while( iss >> c );
      state = 1;
    }
    else{
      tmpRows.push_back(c);
      tmpCells.resize( tmpCells.size() + 1 );
      int score;
      while( iss >> score ){
	tmpCells.back().push_back(score);
      }
      if( tmpCells.back().size() != tmpCols.size() ) ERR( "bad score matrix" );
    }
  }

  if( stream.eof() && !stream.bad() && !tmpRows.empty() ){
    stream.clear();
    m.rows.swap(tmpRows);
    m.cols.swap(tmpCols);
    m.cells.swap(tmpCells);
  }

  return stream;
}

std::ostream& operator<<( std::ostream& stream, const ScoreMatrix& m ){
  for( std::size_t i = 0; i < m.cols.size(); ++i ){
    stream << ' ' << std::setw(m.OUTPAD) << m.cols[i];
  }
  stream << '\n';

  for( std::size_t i = 0; i < m.rows.size(); ++i ){
    stream << m.rows[i];
    for( std::size_t j = 0; j < m.cols.size(); ++j ){
      stream << ' ' << std::setw(m.OUTPAD) << m.cells[i][j];
    }
    stream << '\n';
  }

  return stream;
}

}  // end namespace cbrc
