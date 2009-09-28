// Copyright 2009 Toshiyuki Sato

#ifndef GENETICCODE_HH
#define GENETICCODE_HH

#include <string>
#include <vector>
#include <iosfwd>

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
  virtual void 			translate( std::vector<uchar>::const_iterator beg,
					   std::vector<uchar>::const_iterator end,
					   std::vector<uchar>::iterator dest );

  static const char* standard;  // the standard genetic code
};

// Convert an amino-acid (translated) coordinate to a DNA coordinate
inline unsigned aaToDna( unsigned aaCoordinate, unsigned frameSize ){
  unsigned frame = aaCoordinate / frameSize;
  unsigned offset = aaCoordinate % frameSize;
  return frame + offset * 3;
}

// Convert a DNA coordinate to an amino-acid (translated) coordinate
inline unsigned dnaToAa( unsigned dnaCoordinate, unsigned frameSize ){
  unsigned frame = dnaCoordinate % 3;
  unsigned offset = dnaCoordinate / 3;
  return frame * frameSize + offset;
}

} // end namespace cbrc

#endif
