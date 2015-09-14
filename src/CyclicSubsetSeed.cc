// Copyright 2009, 2010, 2013, 2014 Martin C. Frith

#include "CyclicSubsetSeed.hh"
#include "io.hh"
#include "stringify.hh"
#include <algorithm>  // sort
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <cctype>  // toupper, tolower
#include <stddef.h>  // size_t
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

using namespace cbrc;

std::string CyclicSubsetSeed::stringFromName( const std::string& name ){
  if( name == "BISF" ) return "\
1  CT A G\n\
0  ACGT\n\
1111110101100\n\
";

  if( name == "BISR" ) return "\
1  AG C T\n\
0  ACGT\n\
1111110101100\n\
";

  // From MC Frith & L Noe (2014) Nucleic Acids Research,
  // Supplementary Table 11, row 12.
  if( name == "MAM4" ) return "\
1  A C G T\n\
0  ACGT\n\
T  AG CT\n\
11100TT01T00T10TTTT\n\
TTTT110TT0T001T0T1T1\n\
11TT010T01TT0001T\n\
11TT10T1T101TT\n\
";

  // From MC Frith & L Noe (2014) Nucleic Acids Research,
  // Supplementary Table 12, second-last row.
  if( name == "MAM8" ) return "\
1  A C G T\n\
0  ACGT\n\
T  AG CT\n\
1101T1T0T1T00TT1TT\n\
1TTTTT010TT0TT01011TTT\n\
1TTTT10010T011T0TTTT1\n\
111T011T0T01T100\n\
1T10T100TT01000TT01TT11\n\
111T101TT000T0T10T00T1T\n\
111100T011TTT00T0TT01T\n\
1T1T10T1101101\n\
";

  if( name == "MURPHY10" ) return "\
1  ILMV FWY A C G H P KR ST DENQ\n\
1\n\
";

  return slurp( name );
}

std::string
CyclicSubsetSeed::stringFromPatterns( const std::string& patterns,
				      const std::string& sequenceLetters ){
  std::string spacedLetters;
  for( size_t i = 0; i < sequenceLetters.size(); ++i ){
    if( i > 0 ) spacedLetters += ' ';
    spacedLetters += sequenceLetters[i];
  }

  std::string p = patterns;
  for( size_t i = 0; i < p.size(); ++i )
    switch( p[i] ){
    case ',':
      p[i] = '\n'; break;
    case '#':
      p[i] = '1'; break;
    case '_':
    case '-':
      p[i] = '0'; break;
    case 't':
    case '@':
      p[i] = 'T'; break;
    }

  return "\
1  " + spacedLetters + "\n\
0  " + sequenceLetters + "\n\
T  AG CT\n\
" + p;
}

bool CyclicSubsetSeed::nextPattern( std::istream& in,
				    std::vector< std::string >& seedAlphabet,
				    std::string& pattern ){
  while ( getline( in, pattern ) ){
    std::istringstream iss( pattern );
    std::string x, y;
    iss >> x;
    if( x.empty() || x[0] == '#' ) continue;
    iss >> y;
    if( y.empty() ) return true;
    if( x.size() > 1 ) ERR( "bad seed line: " + pattern );
    seedAlphabet.push_back( pattern );
  }
  return false;
}

static size_t findSeedLetter( const std::vector< std::string >& seedAlphabet,
			      char seedLetter ){
  // go backwards, so that newer definitions override older ones:
  for( size_t j = seedAlphabet.size(); j > 0; ){
    --j;
    const std::string& s = seedAlphabet[j];
    assert( !s.empty() );
    if( s[0] == seedLetter ) return j;
  }
  ERR( "unknown symbol in seed pattern: " + stringify(seedLetter) );
}

void CyclicSubsetSeed::init( const std::vector< std::string >& seedAlphabet,
			     const std::string& pattern,
			     bool isMaskLowercase,
			     const uchar letterCode[] ){
  clear();
  for( size_t i = 0; i < pattern.size(); ++i ){
    char seedLetter = pattern[i];
    if( std::isspace(seedLetter) ) continue;
    size_t j = findSeedLetter( seedAlphabet, seedLetter );
    std::istringstream iss( seedAlphabet[j] );
    iss >> seedLetter;
    appendPosition( iss, isMaskLowercase, letterCode );
  }
}

static void addLetter( uchar numbersToSubsets[], uchar letter, uchar subsetNum,
		       const uchar letterCode[] ){
  uchar number = letterCode[letter];
  if( number >= CyclicSubsetSeed::MAX_LETTERS )
    ERR( "bad symbol in subset-seed: " + stringify(letter) );
  if( numbersToSubsets[number] < CyclicSubsetSeed::DELIMITER )
    ERR( "repeated symbol in subset-seed: " + stringify(letter) );
  numbersToSubsets[number] = subsetNum;
}

void CyclicSubsetSeed::appendPosition( std::istream& inputLine,
				       bool isMaskLowercase,
				       const uchar letterCode[] ){
  std::string inputWord;
  std::vector<std::string> subsetList;
  std::vector<uchar> numbersToSubsets( MAX_LETTERS, DELIMITER );

  for( unsigned subsetNum = 0; inputLine >> inputWord; ++subsetNum ){
    assert( subsetNum < DELIMITER );
    std::string subset;

    for( size_t i = 0; i < inputWord.size(); ++i ){
      uchar upper = std::toupper( inputWord[i] );
      uchar lower = std::tolower( inputWord[i] );
      addLetter( &numbersToSubsets[0], upper, subsetNum, letterCode );
      subset += upper;
      if( !isMaskLowercase && lower != upper ){
	addLetter( &numbersToSubsets[0], lower, subsetNum, letterCode );
      }
    }

    sort( subset.begin(), subset.end() );  // canonicalize
    subsetList.push_back( subset );
  }

  subsetLists.push_back( subsetList );
  subsetMaps.insert( subsetMaps.end(),
		     numbersToSubsets.begin(), numbersToSubsets.end() );
}

void CyclicSubsetSeed::writePosition( std::ostream& out,
				      unsigned position ) const{
  assert( position < subsetLists.size() );
  for( size_t i = 0; i < subsetLists[position].size(); ++i ){
    if( i > 0 ) out << ' ';
    out << subsetLists[position][i];
  }
}
