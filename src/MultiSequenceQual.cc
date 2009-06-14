// Copyright 2009 Martin C. Frith

#include "MultiSequence.hh"
#include "stringify.hh"
#include <algorithm>  // max_element
#include <limits>  // numeric_limits

// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

std::istream&
MultiSequence::appendFromFastq( std::istream& stream, std::size_t maxBytes ){
  const uchar padQualityScore = 0;  // dummy value: should never get used

  // initForAppending:
  if( qualityScores.empty() )
    qualityScores.insert( qualityScores.end(), padSize, padQualityScore );

  // reinitForAppending:
  if( qualityScores.size() > seq.size() )
    qualityScores.erase( qualityScores.begin(),
			 qualityScores.end() - seq.size() );

  if( isFinished() ){
    readFastaName(stream);
    if( !stream ) return stream;

    uchar c;

    // don't bother to obey maxBytes exactly: harmless for short sequences
    while( stream >> c && c != '+' ){  // skips whitespace
      seq.push_back(c);
    }

    stream.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    while( qualityScores.size() < seq.size() && stream >> c ){  // skips WS
      qualityScores.push_back(c);
    }

    if( seq.size() != qualityScores.size() )
      throw std::runtime_error("bad FASTQ data");
  }

  if( isFinishable(maxBytes) ){
    finish();
    qualityScores.insert( qualityScores.end(), padSize, padQualityScore );
  }

  return stream;
}

std::istream&
MultiSequence::appendFromPrb( std::istream& stream, std::size_t maxBytes,
			      unsigned alphSize, const uchar decode[] ){
  const uchar padQualityScore = 0;  // dummy value: should never get used
  indexT qualPadSize = padSize * alphSize;
  indexT qualSize = seq.size() * alphSize;

  // initForAppending:
  if( qualityScores.empty() )
    qualityScores.insert( qualityScores.end(), qualPadSize, padQualityScore );

  // reinitForAppending:
  if( qualityScores.size() > qualSize )
    qualityScores.erase( qualityScores.begin(),
			 qualityScores.end() - qualSize );

  if( isFinished() ){
    std::string line;
    getline( stream, line );  // slow but simple
    if( !stream ) return stream;

    // give the sequence a boring name:
    static std::size_t lineCount = 0;
    std::string name = stringify( ++lineCount );
    names.insert( names.end(), name.begin(), name.end() );
    nameEnds.push_back( names.size() );

    std::istringstream iss(line);
    int q;
    while( iss >> q ){
      if( q < -64 || q > 62 )
	throw std::runtime_error( "quality score too large: " + stringify(q) );
      qualityScores.push_back( q + 64 );  // ASCII-encode the quality score
    }

    if( qualityScores.size() % alphSize != 0 )
      throw std::runtime_error("bad PRB data");

    for( CI(uchar) i = qualityScores.begin() + qualSize;
	 i < qualityScores.end(); i += alphSize ){
      unsigned maxIndex = std::max_element( i, i + alphSize ) - i;
      seq.push_back( decode[ maxIndex ] );
    }
  }

  if( isFinishable(maxBytes) ){
    finish();
    qualityScores.insert( qualityScores.end(), qualPadSize, padQualityScore );
  }

  return stream;
}
