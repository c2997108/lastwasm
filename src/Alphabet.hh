// Copyright 2008, 2009, 2010 Martin C. Frith

// This struct maps characters to codes (small integers) and back.

// We allow for both "proper" letters (e.g. ACGT for DNA), which get
// the lowest codes, and "improper" letters.  This is because real
// sequence data includes additional letters (e.g. ambiguous bases),
// so we have to handle them.  In addition, the space character
// represents a special delimiter.

#ifndef ALPHABET_HH
#define ALPHABET_HH

#include <string>
#include <iosfwd>

namespace cbrc{

struct Alphabet{
  typedef unsigned char uchar;

  static const char* dna;
  static const char* protein;

  static const unsigned capacity = 256;

  // make an Alphabet from a string containing the "proper" letters
  void fromString( const std::string& alphString );

  // translate (encode) a sequence of letters to numbers, in place
  void tr( uchar* beg, uchar* end ) const;

  // reverse-translate (decode) a sequence of numbers to letters
  std::string rtString( const uchar* beg, const uchar* end ) const;

  // reverse and complement a sequence of numbers, in place
  void rc( uchar* beg, uchar* end ) const;

  std::string letters;    // the "proper" letters, e.g. ACGT for DNA
  unsigned size;          // same as letters.size(): excludes delimiters
  uchar encode[capacity];  // translate ASCII letters to codes (small integers)
  uchar decode[capacity];  // translate codes to ASCII letters
  uchar canonical[capacity];   // translate lowercase codes to uppercase codes
  uchar complement[capacity];  // translate DNA codes to their complements

  void init();
  void addLetters( const std::string& lettersToAdd, unsigned& code );
  void makeComplement();
};

std::ostream& operator<<( std::ostream& s, const Alphabet& a );
std::istream& operator>>( std::istream& s, Alphabet& a );

}  // end namespace
#endif
