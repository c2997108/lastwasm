// Copyright 2009 Martin C. Frith

#include "CyclicSubsetSeed.hh"
#include "stringify.hh"
#include <fstream>
#include <sstream>
#include <algorithm>  // sort
#include <stdexcept>
#include <cassert>
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

using namespace cbrc;

// This will surely get replaced by a more interesting seed, after
// some testing:
const char* CyclicSubsetSeed::proteinSeed =
  "A C D E F G H I K L M N P Q R S T V W Y";

void CyclicSubsetSeed::fromFile( const std::string& fileName,
				 bool isMaskLowercase,
				 const uchar letterCode[] ){
  std::ifstream f( fileName.c_str() );
  fromStream( f, isMaskLowercase, letterCode );
  if( f.bad() || !f.eof() ) ERR( "can't read file: " + fileName );
}

void CyclicSubsetSeed::fromString( const std::string& s,
				   bool isMaskLowercase,
				   const uchar letterCode[] ){
  std::istringstream iss(s);
  fromStream( iss, isMaskLowercase, letterCode );
}

static bool isBlankOrComment( std::string line ){
  std::istringstream iss(line);
  char c;
  if( !(iss >> c) ) return true;
  if( c == '#' ) return true;
  return false;
}

void CyclicSubsetSeed::fromStream( std::istream& stream,
				   bool isMaskLowercase,
				   const uchar letterCode[] ){
  std::string line;
  while( std::getline( stream, line ) ){
    if( isBlankOrComment(line) ) continue;
    std::istringstream iss(line);
    appendPosition( iss, isMaskLowercase, letterCode );
  }
}

static std::string exactSeed( const std::string& letters ){
  std::string result;
  for( unsigned i = 0; i < letters.size(); ++i ){
    if( i > 0 ) result += ' ';
    result += letters[i];
  }
  return result;
}

void CyclicSubsetSeed::fromSpacedSeed( const std::string& spacedSeed,
				       const std::string& letters,
				       bool isMaskLowercase,
				       const uchar letterCode[] ){
  std::string es = exactSeed(letters);

  for( unsigned i = 0; i < spacedSeed.size(); ++i ){
    if( spacedSeed[i] == '1' ){
      std::istringstream iss(es);
      appendPosition( iss, isMaskLowercase, letterCode );
    }
    else{
      std::istringstream iss(letters);
      appendPosition( iss, isMaskLowercase, letterCode );
    }
  }
}

void CyclicSubsetSeed::addLetter( std::vector<uchar>& numbersToSubsets,
				  uchar letter, uchar subsetNum,
				  const uchar letterCode[] ){
  uchar number = letterCode[letter];
  if( number >= CyclicSubsetSeed::MAX_LETTERS )
    ERR( "bad symbol in subset-seed: " + stringify(letter) );
  if( numbersToSubsets[number] < CyclicSubsetSeed::DELIMITER )
    ERR( "repeated symbol in subset-seed: "  + stringify(letter) );
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

    for( unsigned i = 0; i < inputWord.size(); ++i ){
      uchar upper = std::toupper( inputWord[i] );
      addLetter( numbersToSubsets, upper, subsetNum, letterCode );
      subset += upper;
      if( !isMaskLowercase ){
	uchar lower = std::tolower( inputWord[i] );
	addLetter( numbersToSubsets, lower, subsetNum, letterCode );
      }
    }

    std::sort( subset.begin(), subset.end() );  // canonicalize
    subsetList.push_back( subset );
  }

  subsetLists.push_back( subsetList );
  subsetMaps.insert( subsetMaps.end(),
		     numbersToSubsets.begin(), numbersToSubsets.end() );
}

void CyclicSubsetSeed::writePosition( std::ostream& out,
				      unsigned position ) const{
  assert( position < subsetLists.size() );
  for( unsigned i = 0; i < subsetLists[position].size(); ++i ){
    if( i > 0 ) out << ' ';
    out << subsetLists[position][i];
  }
}
