// Copyright 2008, 2009, 2010 Martin C. Frith

#include "Alphabet.hh"
#include <istream>
#include <ostream>
#include <algorithm>
#include <stdexcept>

using namespace cbrc;

const char* Alphabet::dna = "ACGT";

// U=selenocysteine, O=pyrrolysine, *=stop?
const char* Alphabet::protein = "ACDEFGHIKLMNPQRSTVWY";

void Alphabet::fromString( const std::string& alphString ){
  letters = alphString;
  init();
}

void Alphabet::tr( uchar* beg, uchar* end ) const{
  for( /* noop */; beg < end; ++beg ){
    uchar code = encode[ *beg ];
    if( code == dummyCode ){
      throw std::runtime_error( std::string("bad symbol: ") + char(*beg) );
    }
    *beg = code;
  }
}

std::string Alphabet::rtString( const uchar* beg, const uchar* end ) const{
  std::string s;
  for( /* noop */; beg < end; ++beg ){
    s += decode[ *beg ];
  }
  return s;
}

void Alphabet::rc( uchar* beg, uchar* end ) const{
  std::reverse( beg, end );

  for( /* noop */; beg < end; ++beg ){
    *beg = complement[ *beg ];
  }
}

void Alphabet::init(){
  // std::toupper doesn't work here!
  std::transform( letters.begin(), letters.end(), letters.begin(), toupper );

  std::fill_n( encode, capacity, dummyCode );

  for( unsigned i = 0; i < capacity; ++i ){
    canonical[i] = i;
  }

  unsigned code = 0;
  addLetters( letters, code );
  size = code;
  addLetters( " ", code );  // add space as a delimiter symbol

  if( code - 1 != letters.size() ){
    throw std::runtime_error( "bad alphabet: " + letters );
  }

  addLetters( "ABCDEFGHIJKLMNOPQRSTUVWXYZ", code );
  addLetters( "*", code );  // sometimes appears in protein sequences
  addLetters( "abcdefghijklmnopqrstuvwxyz", code );

  makeComplement();
}

void Alphabet::addLetters( const std::string& lettersToAdd, unsigned& code ){
  for( unsigned i = 0; i < lettersToAdd.size(); ++i ){
    uchar letter = lettersToAdd[i];
    if( encode[letter] == dummyCode ){
      encode[letter] = code;
      decode[code] = letter;
      canonical[code] = encode[ std::toupper(letter) ];
      ++code;
    }
  }
}

void Alphabet::makeComplement(){
  static const char* x = "ACGTRYSWKMBDHVN";
  static const char* y = "TGCAYRSWMKVHDBN";

  for( unsigned i = 0; i < capacity; ++i ){
    complement[i] = i;
  }

  for( unsigned i = 0; x[i] && y[i]; ++i ){
    uchar xUpper = std::toupper( x[i] );
    uchar yUpper = std::toupper( y[i] );
    complement[ encode[xUpper] ] = encode[yUpper];
    uchar xLower = std::tolower( x[i] );
    uchar yLower = std::tolower( y[i] );
    complement[ encode[xLower] ] = encode[yLower];
  }
}

std::ostream& cbrc::operator<<( std::ostream& s, const Alphabet& a ){
  return s << a.letters;
}

std::istream& cbrc::operator>>( std::istream& s, Alphabet& a ){
  s >> a.letters;
  if( s ) a.init();
  return s;
}
