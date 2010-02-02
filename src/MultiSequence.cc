// Copyright 2008, 2009, 2010 Martin C. Frith

#include "MultiSequence.hh"
#include "io.hh"
#include <sstream>
#include <algorithm>  // upper_bound
#include <cassert>

using namespace cbrc;

void MultiSequence::initForAppending( indexT padSizeIn ){
  padSize = padSizeIn;
  seq.v.assign( padSize, ' ' );
  ends.v.assign( 1, padSize );
  names.v.clear();
  nameEnds.v.assign( 1, 0 );
}

void MultiSequence::reinitForAppending(){
  seq.v.erase( seq.v.begin(), seq.v.begin() + ends.v.back() - padSize );
  names.v.erase( names.v.begin(),
                 names.v.begin() + nameEnds.v[ finishedSequences() ] );
  ends.v.resize(1);
  nameEnds.v.resize(1);
  if( !names.v.empty() ) nameEnds.v.push_back( names.v.size() );
}

void MultiSequence::fromFiles( const std::string& baseName, indexT seqCount ){
  ends.m.open( baseName + ".ssp", seqCount + 1 );
  seq.m.open( baseName + ".tis", ends.m.back() );
  nameEnds.m.open( baseName + ".sds", seqCount + 1 );
  names.m.open( baseName + ".des", nameEnds.m.back() );
  padSize = ends.m[0];
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
}

std::istream& MultiSequence::readFastaName( std::istream& stream ){
  char c;
  stream >> c;  // don't check that it's '>': works for FASTQ too
  std::string line, word;
  getline( stream, line );
  std::istringstream iss(line);
  iss >> word;
  if( !stream ) return stream;
  names.v.insert( names.v.end(), word.begin(), word.end() );
  nameEnds.v.push_back( names.v.size() );
  return stream;
}

// probably slower than it could be:
std::istream&
MultiSequence::appendFromFasta( std::istream& stream, std::size_t maxBytes ){
  if( isFinished() ){
    readFastaName(stream);
    if( !stream ) return stream;
  }

  uchar c;
  while( isFinishable(maxBytes) && stream >> c ){  // skips whitespace
    if( c == '>' ){
      stream.unget();
      break;
    }
    seq.v.push_back(c);
  }

  if( isFinishable(maxBytes) ) finish();
  if( stream.eof() && !stream.bad() ) stream.clear();
  return stream;
}

void MultiSequence::finish(){
  assert( !isFinished() );
  seq.v.insert( seq.v.end(), padSize, ' ' );
  ends.v.push_back( seq.v.size() );
}

void MultiSequence::unfinish(){
  assert( isFinished() );
  ends.v.pop_back();
  seq.v.erase( seq.v.end() - padSize, seq.v.end() );
}

bool MultiSequence::isFinishable( std::size_t maxBytes ) const{
  return seq.v.size() + padSize <= maxBytes;
}

MultiSequence::indexT MultiSequence::whichSequence( indexT coordinate ) const{
  const indexT* u = std::upper_bound( ends.begin(), ends.end(), coordinate );
  assert( u != ends.begin() && u != ends.end() );
  return u - ends.begin() - 1;
}

std::string MultiSequence::seqName( indexT seqNum ) const{
  return std::string( names.begin() + nameEnds[ seqNum ],
		      names.begin() + nameEnds[ seqNum + 1 ] );
}
