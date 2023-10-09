// Copyright 2008, 2009, 2010, 2012, 2014 Martin C. Frith

#include "Alphabet.hh"
#include <algorithm>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <cctype>  // toupper, tolower

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

using namespace cbrc;

static const unsigned dummyCode = Alphabet::capacity - 1;

const char* Alphabet::dna = "ACGT";

// U=selenocysteine, O=pyrrolysine, *=stop?
const char* Alphabet::protein = "ACDEFGHIKLMNPQRSTVWY";
const char* Alphabet::proteinWithStop = "ACDEFGHIKLMNPQRSTVWY*";

void Alphabet::tr( uchar* beg, uchar* end, bool isKeepLowercase ) const{
  const uchar* x = isKeepLowercase ? encode : lettersToUppercase;
  for( /* noop */; beg < end; ++beg ){
    uchar code = x[ *beg ];
    if( code == dummyCode ){
      err( std::string("bad symbol in sequence: ") + char(*beg) );
    }
    *beg = code;
  }
}

void Alphabet::init(const std::string &mainLetters, bool is4bit) {
  letters = mainLetters;
  for( std::string::iterator i = letters.begin(); i < letters.end(); ++i )
    *i = std::toupper( *i );

  std::fill_n( encode, capacity, dummyCode );

  unsigned code = 0;
  addLetters( letters, code );
  size = code;
  addLetters( " ", code );  // add space as a delimiter symbol

  if( code - 1 != letters.size() ) err( "bad alphabet: " + letters );

  if (is4bit) addLetters("NRY.acgtnry", code);  // uppercase<=>lowercase: +-8

  addLetters( "ABCDEFGHIJKLMNOPQRSTUVWXYZ", code );
  addLetters( "*", code );  // sometimes appears in protein sequences
  addLetters( "abcdefghijklmnopqrstuvwxyz", code );
  addLetters( ".", code );  // sometimes appears in FASTQ

  initCaseConversions( code );
  makeComplement();
}

void Alphabet::set4bitAmbiguities() {
  static const char a[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*";
  for (const char *i = a; *i; ++i) {
    int up = *i;
    if (encode[up] > 15) {
      int lo = std::tolower(up);
      encode[lo] = encode['n'];
      encode[up] = encode['N'];
      lettersToUppercase[lo] = lettersToUppercase[up] = encode['N'];
    }
  }
}

void Alphabet::addLetters( const std::string& lettersToAdd, unsigned& code ){
  for( unsigned i = 0; i < lettersToAdd.size(); ++i ){
    uchar letter = lettersToAdd[i];
    if( encode[letter] == dummyCode ){
      encode[letter] = code;
      decode[code] = letter;
      ++code;
    }
  }
}

void Alphabet::initCaseConversions( unsigned codeEnd ) {
  for( unsigned i = 0; i < codeEnd; ++i ){
    uchar letter = decode[i];
    numbersToUppercase[i] = encode[ std::toupper(letter) ];
    numbersToLowercase[i] = encode[ std::tolower(letter) ];
  }

  for( unsigned i = codeEnd; i < capacity; ++i ){
    numbersToUppercase[i] = i;
    numbersToLowercase[i] = i;
  }

  if (letters == dna) {  // kludge for RNA:
    encode['U'] = encode['T'];
    encode['u'] = encode['t'];
  }

  for( unsigned i = 0; i < capacity; ++i ){
    lettersToUppercase[i] = numbersToUppercase[ encode[i] ];
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
