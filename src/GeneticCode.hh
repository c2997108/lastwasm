// Copyright 2009 Toshiyuki Sato

#ifndef GENETICCODE_HH
#define GENETICCODE_HH

#include <string>
#include <vector>
#include <iosfwd>
#include <cassert>

namespace cbrc{

class Alphabet;

class GeneticCode
{
 private:
  typedef unsigned char uchar;

  GeneticCode( GeneticCode &c );
  const GeneticCode & operator=( const GeneticCode &c );
  std::string			AAs;
  std::string			Base[3];
  static const int 		NumMember = 54;			// DNA member
  static const int 		UNKNOWN = NumMember*NumMember*NumMember;	// unknown residue
  std::vector<uchar>		genome2residue;
 protected:
  virtual int			codon2number( const uchar* codon );
  virtual int			codon2number2( std::vector<uchar> codon, const Alphabet& dnaAlph );
  friend std::istream& operator>>( std::istream& stream, GeneticCode& codon  );
 public:
  GeneticCode(){
    //    std::cout << "Constructing GeneticCode.\n";
  }
  virtual			~GeneticCode(){
    //    std::cout << "Destructing GeneticCode.\n";
  }
  virtual void 			fromFile( const std::string& codeTable );
  virtual void 			fromString( const std::string& s );
  virtual void			codeTableSet( const Alphabet& aaAlph, const Alphabet& dnaAlph );
  virtual void 			translate( const uchar* beg,
                                           const uchar* end,
					   uchar* dest );

  static const char* standard;  // the standard genetic code
};

// Convert an amino-acid (translated) coordinate to a DNA coordinate
inline unsigned aaToDna( unsigned aaCoordinate, unsigned frameSize ){
  if( frameSize == 0 ) return aaCoordinate;  // for non-translated sequences
  unsigned frame = aaCoordinate / frameSize;
  unsigned offset = aaCoordinate % frameSize;
  return frame + offset * 3;
}

// Convert a DNA coordinate to an amino-acid (translated) coordinate
inline unsigned dnaToAa( unsigned dnaCoordinate, unsigned frameSize ){
  if( frameSize == 0 ) return dnaCoordinate;  // for non-translated sequences
  unsigned frame = dnaCoordinate % 3;
  unsigned offset = dnaCoordinate / 3;
  return frame * frameSize + offset;
}

// Convert begin and end coordinates to a size and a frameshift
inline void sizeAndFrameshift( unsigned beg, unsigned end,
			       unsigned frameSize,  // 0 means not translated
			       unsigned& size, unsigned& frameshift ){
  if( frameSize ){  // if it's a translated sequence:
    unsigned dnaBeg = aaToDna( beg, frameSize );
    unsigned dnaEnd = aaToDna( end, frameSize );
    unsigned dnaSize = dnaEnd - dnaBeg;
    assert( dnaBeg <= dnaEnd + 1 );  // allow a -1 frameshift
    size = ( dnaSize + 1 ) / 3;
    frameshift = ( dnaSize + 3 ) % 3;
  }
  else{  // if it's not a translated sequence:
    assert( beg <= end );
    size = end - beg;
    frameshift = 0;
  }
}

} // end namespace cbrc

#endif
