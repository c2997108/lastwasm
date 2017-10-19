// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

#include "MultiSequence.hh"
#include "io.hh"
#include <sstream>
#include <algorithm>  // upper_bound
#include <cassert>
#include <iterator>  // istreambuf_iterator

using namespace cbrc;

void MultiSequence::initForAppending( indexT padSizeIn ){
  padSize = padSizeIn;
  seq.v.assign( padSize, ' ' );
  ends.v.assign( 1, padSize );
  names.v.clear();
  nameEnds.v.assign( 1, 0 );
  qualityScoresPerLetter = 0;
}

void MultiSequence::reinitForAppending(){
  size_t n = finishedSequences();
  size_t s = padBeg(n);

  seq.v.erase(seq.v.begin(), seq.v.begin() + s);
  names.v.erase(names.v.begin(), names.v.begin() + nameEnds.v[n]);
  ends.v.resize(1);
  nameEnds.v.resize(1);
  if( !names.v.empty() ) nameEnds.v.push_back( names.v.size() );

  qualityScores.v.erase(qualityScores.v.begin(),
			qualityScores.v.begin() + s * qualsPerLetter());

  if (!pssm.empty()) {
    pssm.erase(pssm.begin(), pssm.begin() + s * scoreMatrixRowSize);
  }
}

void MultiSequence::fromFiles( const std::string& baseName, indexT seqCount,
                               size_t qualitiesPerLetter ){
  ends.m.open( baseName + ".ssp", seqCount + 1 );
  seq.m.open( baseName + ".tis", ends.m.back() );
  nameEnds.m.open( baseName + ".sds", seqCount + 1 );
  names.m.open( baseName + ".des", nameEnds.m.back() );
  padSize = ends.m[0];

  qualityScores.m.open( baseName + ".qua",
                        ends.m.back() * qualitiesPerLetter );
  qualityScoresPerLetter = qualitiesPerLetter;
}

void MultiSequence::toFiles( const std::string& baseName ) const{
  memoryToBinaryFile( ends.begin(), ends.end(), baseName + ".ssp" );

  memoryToBinaryFile( seq.begin(), seq.begin() + ends.back(),
		      baseName + ".tis" );

  memoryToBinaryFile( nameEnds.begin(), nameEnds.begin() + ends.size(),
		      baseName + ".sds" );

  memoryToBinaryFile( names.begin(),
		      names.begin() + nameEnds[ finishedSequences() ],
		      baseName + ".des" );

  memoryToBinaryFile( qualityScores.begin(),
                      qualityScores.begin() + ends.back() * qualsPerLetter(),
                      baseName + ".qua" );
}

std::istream& MultiSequence::readFastaName( std::istream& stream ){
  std::string line, word;
  getline( stream, line );
  std::istringstream iss(line);
  iss >> word;
  if( !stream ) return stream;
  addName(word);
  return stream;
}

std::istream&
MultiSequence::appendFromFasta( std::istream& stream, indexT maxSeqLen ){
  if( isFinished() ){
    char c = '>';
    stream >> c;
    if( c != '>' )
      throw std::runtime_error("bad FASTA sequence data: missing '>'");
    readFastaName(stream);
    if( !stream ) return stream;
  }

  std::istreambuf_iterator<char> inpos(stream);
  std::istreambuf_iterator<char> endpos;
  while( inpos != endpos ){
    uchar c = *inpos;
    if( c > ' ' ){  // faster than isspace
      if( c == '>' ) break;  // we have hit the next FASTA sequence
      if( seq.v.size() >= maxSeqLen ) break;
      seq.v.push_back(c);
    }
    ++inpos;
  }

  if( isFinishable(maxSeqLen) ) finish();
  return stream;
}

void MultiSequence::finish(){
  assert( !isFinished() );
  seq.v.insert( seq.v.end(), padSize, ' ' );
  ends.v.push_back( seq.v.size() );
  assert( ends.v.back() == seq.v.size() );
}

MultiSequence::indexT MultiSequence::whichSequence( indexT coordinate ) const{
  const indexT* u = std::upper_bound( ends.begin(), ends.end(), coordinate );
  assert( u != ends.begin() && u != ends.end() );
  return u - ends.begin() - 1;
}

static void reverseComplementPssm(int *beg, int *end,
				  const uchar *complement) {
  while (beg < end) {
    end -= scoreMatrixRowSize;
    for (unsigned i = 0; i < scoreMatrixRowSize; ++i) {
      unsigned j = complement[i];
      if (beg < end || i < j) std::swap(beg[i], end[j]);
    }
    beg += scoreMatrixRowSize;
  }
}

void MultiSequence::reverseComplementOneSequence(indexT seqNum,
						 const uchar *complement) {
  size_t b = seqBeg(seqNum);
  size_t e = seqEnd(seqNum);

  uchar *s = seqWriter();
  std::reverse(s + b, s + e);
  for (size_t i = b; i < e; ++i) {
    s[i] = complement[s[i]];
  }

  reverse(qualityScores.v.begin() + b * qualsPerLetter(),
	  qualityScores.v.begin() + e * qualsPerLetter());

  if (!pssm.empty()) {
    int *p = &pssm[0];
    reverseComplementPssm(p + b * scoreMatrixRowSize,
			  p + e * scoreMatrixRowSize, complement);
  }
}
