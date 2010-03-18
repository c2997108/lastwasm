// Copyright 2009, 2010 Martin C. Frith

#include "MultiSequence.hh"
#include "stringify.hh"
#include <algorithm>  // max_element
#include <limits>  // numeric_limits

// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

std::istream&
MultiSequence::appendFromFastq( std::istream& stream, indexT maxSeqLen ){
  const uchar padQualityScore = 0;  // dummy value: should never get used

  // initForAppending:
  if( qualityScores.empty() )
    qualityScores.insert( qualityScores.end(), padSize, padQualityScore );

  // reinitForAppending:
  if( qualityScores.size() > seq.v.size() )
    qualityScores.erase( qualityScores.begin(),
			 qualityScores.end() - seq.v.size() );

  if( isFinished() ){
    readFastaName(stream);
    if( !stream ) return stream;

    uchar c;

    // don't bother to obey maxSeqLen exactly: harmless for short sequences
    while( stream >> c && c != '+' ){  // skips whitespace
      seq.v.push_back(c);
    }

    stream.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    while( qualityScores.size() < seq.v.size() && stream >> c ){  // skips WS
      qualityScores.push_back(c);
    }

    if( seq.v.size() != qualityScores.size() )
      throw std::runtime_error("bad FASTQ data");
  }

  if( isFinishable(maxSeqLen) ){
    finish();
    qualityScores.insert( qualityScores.end(), padSize, padQualityScore );
  }

  return stream;
}

std::istream&
MultiSequence::appendFromPrb( std::istream& stream, indexT maxSeqLen,
			      unsigned alphSize, const uchar decode[] ){
  const uchar padQualityScore = 0;  // dummy value: should never get used
  std::size_t qualPadSize = padSize * alphSize;
  std::size_t qualSize = seq.v.size() * alphSize;

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
    addName(name);

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
      seq.v.push_back( decode[ maxIndex ] );
    }
  }

  if( isFinishable(maxSeqLen) ){
    finish();
    qualityScores.insert( qualityScores.end(), qualPadSize, padQualityScore );
  }

  return stream;
}
