// Copyright 2008 Martin C. Frith

// This struct maps characters to codes (small integers) and back.
// The mapping can be either case-sensitive or case-insensitive.

// We allow for both "proper" letters (e.g. ACGT for DNA), which get
// the lowest codes, and "improper" letters.  This is because real
// sequence data includes additional letters (e.g. ambiguous bases),
// so we have to handle them.  In addition, the space character
// represents a special delimiter.

#ifndef ALPHABET_HH
#define ALPHABET_HH
#include <vector>
#include <string>
#include <iosfwd>

namespace cbrc{

struct Alphabet{
  typedef unsigned char uchar;

  static const char* dna;
  static const char* protein;
  static const char* all;

  // make an Alphabet from a string containing the "proper" letters
  void fromString( const std::string& alphString );

  // make lowercase letters have the same codes as uppercase letters
  void makeCaseInsensitive();

  void tr( std::vector<uchar>::iterator beg,  // translate (encode) in place
	   std::vector<uchar>::iterator end ) const;

  std::string rtString( std::vector<uchar>::const_iterator beg,  // decode
			std::vector<uchar>::const_iterator end ) const;

  void rc( std::vector<uchar>::iterator beg,  // reverse & complement in place
           std::vector<uchar>::iterator end ) const;

  // reduce the alphabet by merging codes for "improper" letters
  void mergeImproperCodes( std::vector<uchar>::iterator beg,
			   std::vector<uchar>::iterator end ) const;

  std::string letters;    // the "proper" letters, e.g. ACGT for DNA
  unsigned size;          // same as letters.size(): excludes delimiters
  uchar encode[256];      // translate ASCII letters to codes (small integers)
  uchar decode[256];      // translate codes to ASCII letters
  uchar canonical[256];   // translate lowercase codes to uppercase codes
  uchar complement[256];  // translate DNA codes to their complements

  void init();
  uchar newCode( uchar letter, uchar code );
  void makeComplement();
};

std::ostream& operator<<( std::ostream& s, const Alphabet& a );
std::istream& operator>>( std::istream& s, Alphabet& a );

}  // end namespace cbrc
#endif  // ALPHABET_HH
