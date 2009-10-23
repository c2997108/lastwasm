// Copyright 2008, 2009 Martin C. Frith

#include "Alphabet.hh"
#include <istream>
#include <ostream>
#include <algorithm>
#include <stdexcept>

using namespace cbrc;

const char* Alphabet::dna = "ACGT";

// U=selenocysteine, O=pyrrolysine, *=stop?
const char* Alphabet::protein = "ACDEFGHIKLMNPQRSTVWY";

// need to allow "*", because blosum62 includes it:
const char* Alphabet::all = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*";

void Alphabet::fromString( const std::string& alphString ){
  letters = alphString;
  init();
}

void Alphabet::tr( std::vector<uchar>::iterator beg,
		   std::vector<uchar>::iterator end ) const{
  for( /* noop */; beg < end; ++beg ){
    uchar code = encode[ *beg ];
    if( code == 255 ){
      throw std::runtime_error( std::string("bad symbol: ") + char(*beg) );
    }
    *beg = code;
  }
}

std::string Alphabet::rtString( std::vector<uchar>::const_iterator beg,
				std::vector<uchar>::const_iterator end ) const{
  std::string s;
  for( /* noop */; beg < end; ++beg ){
    s += decode[ *beg ];
  }
  return s;
}

void Alphabet::rc( std::vector<uchar>::iterator beg,
		   std::vector<uchar>::iterator end ) const{
  std::reverse( beg, end );

  for( /* noop */; beg < end; ++beg ){
    *beg = complement[ *beg ];
  }
}

void Alphabet::init(){
  size = letters.size();

  // std::toupper doesn't work here!
  std::transform( letters.begin(), letters.end(), letters.begin(), toupper );

  std::fill_n( encode, 256, 255 );  // fill encode with value 255

  for( unsigned i = 0; i < 256; ++i ){
    canonical[i] = i;
  }

  uchar code = 0;

  // add the "real" alphabet letters:
  for( std::string::iterator i = letters.begin(); i < letters.end(); ++i ){
    uchar letter = *i;
    if( letter == ' ' || encode[ letter ] != 255 ){
      throw std::runtime_error( "bad alphabet: " + letters );
    }
    code = newCode( letter, code );
  }

  // add space as a delimiter symbol:
  code = newCode( ' ', code );

  // add any remaining letters (e.g. ambiguous nucleotides):
  for( const char* i = all; *i; ++i ){
    uchar letter = *i;
    if( encode[ letter ] == 255 ) code = newCode( letter, code );
  }

  // add lowercase letters:
  for( const char* i = all; *i; ++i ){
    uchar letter = std::tolower( *i );
    if( encode[ letter ] == 255 ) code = newCode( letter, code );
  }

  makeComplement();
}

Alphabet::uchar Alphabet::newCode( uchar letter, uchar code ){
  encode[letter] = code;
  decode[code] = letter;
  canonical[code] = encode[ std::toupper(letter) ];
  return code + 1;
}

void Alphabet::makeComplement(){
  static const char* x = "ACGTRYSWKMBDHVN";
  static const char* y = "TGCAYRSWMKVHDBN";

  for( unsigned i = 0; i < 256; ++i ){
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
