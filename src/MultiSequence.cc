// Copyright 2008, 2009, 2010 Martin C. Frith

#include "MultiSequence.hh"
#include "io.hh"
#include <sstream>
#include <algorithm>  // upper_bound
#include <cassert>

using namespace cbrc;

void MultiSequence::initForAppending( indexT padSizeIn ){
  padSize = padSizeIn;
  seq.assign( padSize, ' ' );
  ends.assign( 1, padSize );
  names.clear();
  nameEnds.assign( 1, 0 );
}

void MultiSequence::reinitForAppending(){
  seq.erase( seq.begin(), seq.begin() + ends.back() - padSize );
  names.erase( names.begin(),
	       names.begin() + nameEnds[ finishedSequences() ] );
  ends.resize(1);
  nameEnds.resize(1);
  if( !names.empty() ) nameEnds.push_back( names.size() );
}

void MultiSequence::fromFiles( const std::string& baseName, indexT seqCount ){
  ends.resize( seqCount + 1 );  // unwanted zero-fill
  vectorFromBinaryFile( ends, baseName + ".ssp" );

  seq.resize( ends.back() );  // unwanted zero-fill
  vectorFromBinaryFile( seq, baseName + ".tis" );

  nameEnds.resize( seqCount + 1 );  // unwanted zero-fill
  vectorFromBinaryFile( nameEnds, baseName + ".sds" );

  names.resize( nameEnds.back() );  // unwanted zero-fill
  vectorFromBinaryFile( names, baseName + ".des" );

  padSize = ends[0];
}

void MultiSequence::toFiles( const std::string& baseName ) const{
  vectorToBinaryFile( ends, baseName + ".ssp" );

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
  names.insert( names.end(), word.begin(), word.end() );
  nameEnds.push_back( names.size() );
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
    seq.push_back(c);
  }

  if( isFinishable(maxBytes) ) finish();
  if( stream.eof() && !stream.bad() ) stream.clear();
  return stream;
}

void MultiSequence::finish(){
  assert( !isFinished() );
  seq.insert( seq.end(), padSize, ' ' );
  ends.push_back( seq.size() );
}

void MultiSequence::unfinish(){
  assert( isFinished() );
  ends.pop_back();
  seq.erase( seq.end() - padSize, seq.end() );
}

bool MultiSequence::isFinishable( std::size_t maxBytes ) const{
  return seq.size() + padSize <= maxBytes;
}

MultiSequence::indexT MultiSequence::whichSequence( indexT coordinate ) const{
  std::vector<indexT>::const_iterator u =
    std::upper_bound( ends.begin(), ends.end(), coordinate );
  assert( u != ends.begin() && u != ends.end() );
  return u - ends.begin() - 1;
}

std::string MultiSequence::seqName( indexT seqNum ) const{
  return std::string( names.begin() + nameEnds[ seqNum ],
		      names.begin() + nameEnds[ seqNum + 1 ] );
}
