// Copyright 2008, 2009 Martin C. Frith

// Read fasta-format sequences; construct a suffix array of them; and
// write the results to files.

#include "LastdbArguments.hh"
#include "SuffixArray.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "PeriodicSpacedSeed.hh"
#include "io.hh"
#include "stringify.hh"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE

#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastdb: " << x << '\n'

using namespace cbrc;

typedef unsigned indexT;

// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet( Alphabet& alph, const LastdbArguments& args ){
  if( !args.userAlphabet.empty() )  alph.fromString( args.userAlphabet );
  else if( args.isProtein )         alph.fromString( alph.protein );
  else                              alph.fromString( alph.dna );
  if( !args.isCaseSensitive )       alph.makeCaseInsensitive();
}

// Write the .prj file for the whole database
void writeOuterPrj( const std::string& fileName, const Alphabet& alph,
		    const PeriodicSpacedSeed& seed, unsigned volumes ){
  std::ofstream f( fileName.c_str() );
  f << "version=" <<
#include "version.hh"
    << '\n';
  f << "alphabet=" << alph << '\n';
  f << "spacedseed=" << seed << '\n';
  f << "volumes=" << volumes << '\n';
  if( !f ) throw std::runtime_error("can't write file: " + fileName);
}

// Write a per-volume .prj file, with info about a database volume
void writeInnerPrj( const std::string& fileName,
		    const MultiSequence& multi, const SuffixArray& sa ){
  std::ofstream f( fileName.c_str() );
  f << "totallength=" << multi.ends.back() << '\n';
  f << "specialcharacters=" << multi.ends.back() - sa.indexSize() << '\n';
  f << "numofsequences=" << multi.finishedSequences() << '\n';
  f << "prefixlength=" << sa.maxBucketPrefix() << '\n';
  if( !f ) throw std::runtime_error("can't write file: " + fileName);
}

// Make one database volume, from one batch of sequences
void makeVolume( MultiSequence& multi, SuffixArray& sa,
		 const LastdbArguments& args, const Alphabet& alph,
		 unsigned volumeNumber ){
  std::string baseName = args.lastdbName + stringify(volumeNumber);

  LOG( "writing tis, des, ssp, sds..." );
  multi.toFiles( baseName );

  LOG( "recoding..." );
  alph.mergeImproperCodes( multi.seq.begin(), multi.seq.end() );

  LOG( "sorting..." );
  sa.sortIndex();

  LOG( "bucketing..." );
  sa.makeBuckets( args.bucketDepth );

  LOG( "writing suf, bck..." );
  sa.toFiles( baseName );

  LOG( "writing prj..." );
  writeInnerPrj( baseName + ".prj", multi, sa );

  LOG( "done!" );
}

// Read the next sequence, adding it to the MultiSequence and the SuffixArray
std::istream&
appendFromFasta( MultiSequence& multi, SuffixArray& sa,
		 const LastdbArguments& args, const Alphabet& alph,
		 std::istream& in ){
  std::size_t maxSeqBytes = args.volumeSize - sa.indexBytes();
  if( args.volumeSize < sa.indexBytes() ) maxSeqBytes = 0;
  if( multi.finishedSequences() == 0 ) maxSeqBytes = std::size_t(-1);

  indexT oldSeqSize = multi.seq.size();

  multi.appendFromFasta( in, maxSeqBytes );

  // encode the newly-read sequence
  alph.tr( multi.seq.begin() + oldSeqSize, multi.seq.end() );

  if( in && multi.isFinished() ){
    std::size_t maxIndexBytes = args.volumeSize - multi.seq.size();
    if( args.volumeSize < multi.seq.size() ) maxIndexBytes = 0;
    if( multi.finishedSequences() == 1 ) maxIndexBytes = std::size_t(-1);

    if( !sa.addIndices( *(multi.ends.end() - 2), *(multi.ends.end() - 1),
			args.indexStep, maxIndexBytes ) ){
      multi.unfinish();
    }
  }

  return in;
}

void lastdb( int argc, char** argv ){
  LastdbArguments args;
  args.fromArgs( argc, argv );
  Alphabet alph;
  PeriodicSpacedSeed seed;
  MultiSequence multi;
  makeAlphabet( alph, args );
  seed.fromString( args.spacedSeed );
  multi.initForAppending( seed.maxOffset );
  alph.tr( multi.seq.begin(), multi.seq.end() );
  SuffixArray sa( multi.seq, seed.offsets, alph.size );
  unsigned volumeNumber = 0;

  for( char** i = argv + args.inputStart; i < argv + argc; ++i ){
    LOG( "reading " << *i << "..." );
    std::ifstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );

    while( appendFromFasta( multi, sa, args, alph, in ) ){
      if( !multi.isFinished() ){
	makeVolume( multi, sa, args, alph, volumeNumber++ );
	sa.clear();
	multi.reinitForAppending();
      }
    }
  }

  if( multi.finishedSequences() > 0 ){
    makeVolume( multi, sa, args, alph, volumeNumber++ );
  }

  writeOuterPrj( args.lastdbName + ".prj", alph, seed, volumeNumber );
}

int main( int argc, char** argv )
try{
  lastdb( argc, argv );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << "lastdb: out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << "lastdb: " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
