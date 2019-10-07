// Copyright 2009 Toshiyuki Sato

#include "GeneticCode.hh"
#include "GeneticCodeData.hh"
#include "Alphabet.hh"
#include "zio.hh"

#include <cctype>  // toupper, tolower, islower
#include <fstream>
#include <sstream>
#include <stdexcept>
//#include <iostream>  // for debugging

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

namespace cbrc{

std::string GeneticCode::stringFromName(const std::string &name) {
  for (size_t i = 0; i < COUNTOF(geneticCodes); ++i)
    if (name == geneticCodes[i].name)
      return geneticCodes[i].text;
  return slurp(name.c_str());
}

void GeneticCode::fromString( const std::string& s ){
  std::istringstream iss(s);
  iss >> *this;
}

//
void GeneticCode::codeTableSet( const Alphabet& aaAlph, const Alphabet& dnaAlph )
{
  uchar codon[3];

  genome2residue.assign( UNKNOWN, 'X' );

  for( size_t i = 0 ; i < AAs.size() ; i++ ){
    char aminoAcid = std::toupper( AAs[i] );

    for( int x = 0; x < 2; ++x ){
      codon[0] = x ? std::tolower( Base[0][i] ) : std::toupper( Base[0][i] );

      for( int y = 0; y < 2; ++y ){
	codon[1] = y ? std::tolower( Base[1][i] ) : std::toupper( Base[1][i] );

	for( int z = 0; z < 2; ++z ){
	  codon[2] = z ? std::tolower( Base[2][i] ) : std::toupper( Base[2][i] );

	  int c = codon2number2( codon, dnaAlph );
	  genome2residue[c] = aminoAcid;
	}
      }
    }
  }

  // codons containing DNA delimiters, or lowercase bases
  for( int i = 0 ; i < NumMember ; i++ ){
    codon[0]= dnaAlph.decode[i];

    for( int j = 0 ; j < NumMember ; j++ ){
      codon[1]= dnaAlph.decode[j];

      for( int k = 0 ; k < NumMember ; k++ ){
	codon[2]= dnaAlph.decode[k];

	int c = codon2number2( codon, dnaAlph );

	if( codon[0] == ' ' || codon[1] == ' ' || codon[2] == ' ' ){
	  genome2residue[c] = ' ';  // delimiter
	}
	else if( std::islower( codon[0] ) ||
		 std::islower( codon[1] ) ||
		 std::islower( codon[2] ) ){
	  genome2residue[c] = std::tolower( genome2residue[c] );
	}
      }
    }
  }

  aaAlph.tr( &genome2residue.front(), &genome2residue.back() + 1 );

  return;
}

//
void GeneticCode::translate( const uchar* beg, const uchar* end,
			     uchar* dest ) const{
  size_t size = end - beg;

  for( size_t i = 0 ; i < 3 ; i++ ){
    for( size_t j = i ; j+2 < size ; j+=3 ){
      *dest++ = translation( beg + j );
    }

    // this ensures that each reading frame has exactly the same size:
    if( i > size % 3 ) *dest++ = 20;
  }

  // this ensures that the size of the translated sequence is exactly "size":
  if( size % 3 > 0 ) *dest++ = 20;
  if( size % 3 > 1 ) *dest++ = 20;

  return;
}

//
int GeneticCode::codon2number2( const uchar* codon, const Alphabet& dnaAlph ){
  uchar c[3] = { codon[0], codon[1], codon[2] };
  dnaAlph.tr( c, c + 3 );
  return codon2number( c );
}

//
std::istream& operator>>( std::istream& stream, GeneticCode& codon  ){
  std::string	readbuf, label, dummy;

  while( getline(stream,readbuf) ){
    std::istringstream readline(readbuf.c_str());
    readline >> label;

    if( label == "AAs" ){
      readline >> dummy;
      readline >> codon.AAs;
    }
    else if( label == "Base1" ){
      readline >> dummy;
      readline >> codon.Base[0];
    }
    else if( label == "Base2" ){
      readline >> dummy;
      readline >> codon.Base[1];
    }
    else if( label == "Base3" ){
      readline >> dummy;
      readline >> codon.Base[2];
    }
  }

  if( codon.AAs.size() != codon.Base[0].size() ||
      codon.AAs.size() != codon.Base[1].size() ||
      codon.AAs.size() != codon.Base[2].size() ){
    throw std::runtime_error( "bad genetic code table" );
  }

  return stream;
}

} // end namespace cbrc
