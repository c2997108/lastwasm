// Copyright 2008, 2009, 2010 Martin C. Frith

// Read fasta-format sequences; construct a suffix array of them; and
// write the results to files.

#include "LastdbArguments.hh"
#include "SubsetSuffixArray.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "CyclicSubsetSeed.hh"
#include "io.hh"
#include "stringify.hh"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <cstring>  // strchr
#include <numeric>  // accumulate

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastdb: " << x << '\n'

using namespace cbrc;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet( Alphabet& alph, const LastdbArguments& args ){
  if( !args.userAlphabet.empty() )  alph.fromString( args.userAlphabet );
  else if( args.isProtein )         alph.fromString( alph.protein );
  else                              alph.fromString( alph.dna );
}

// Does the first sequence look like it isn't really DNA?
bool isDubiousDna( const Alphabet& alph, const MultiSequence& multi ){
  typedef unsigned char uchar;
  const uchar* seq = multi.seqReader() + multi.seqBeg(0);
  unsigned dnaCount = 0;

  for( indexT i = 0; i < 100; ++i ){  // look at the first 100 letters
    uchar c = alph.canonical[ seq[i] ];
    if( c == alph.size ) return false;  // we hit the end of the sequence early
    if( c < alph.size || c == alph.encode[ 'N' ] ) ++dnaCount;
  }

  if( dnaCount < 90 ) return true;  // more than 10% unexpected letters
  else return false;
}

// Set up a subset seed, based on the user options
void makeSubsetSeed( CyclicSubsetSeed& seed, const LastdbArguments& args,
		     const Alphabet& alph ){
  if( !args.subsetSeedFile.empty() ){
    seed.fromFile( args.subsetSeedFile, args.isCaseSensitive, alph.encode );
  }
  else if( !args.spacedSeed.empty() ){
    seed.fromSpacedSeed( args.spacedSeed, alph.letters,
			 args.isCaseSensitive, alph.encode );
  }
  else if( alph.letters == alph.dna ){
    seed.fromString( seed.yassSeed, args.isCaseSensitive, alph.encode );
  }
  else if( alph.letters == alph.protein ){
    seed.fromString( seed.proteinSeed, args.isCaseSensitive, alph.encode );
  }
  else{
    seed.fromSpacedSeed( "1", alph.letters,
                         args.isCaseSensitive, alph.encode );
  }
}

// Write the .prj file for the whole database
void writeOuterPrj( const std::string& fileName, const LastdbArguments& args,
                    const Alphabet& alph, countT sequenceCount,
                    const std::vector<countT>& letterCounts,
		    const CyclicSubsetSeed& seed, unsigned volumes ){
  countT letterTotal = std::accumulate( letterCounts.begin(),
                                        letterCounts.end(), countT(0) );

  std::ofstream f( fileName.c_str() );
  f << "version=" <<
#include "version.hh"
    << '\n';
  f << "alphabet=" << alph << '\n';
  f << "numofsequences=" << sequenceCount << '\n';
  f << "numofletters=" << letterTotal << '\n';
  f << "letterfreqs=";
  for( unsigned i = 0; i < letterCounts.size(); ++i ){
    if( i > 0 ) f << ' ';
    f << letterCounts[i];
  }
  f << '\n';

  if( !args.isCountsOnly ){
    f << "masklowercase=" << args.isCaseSensitive << '\n';
    f << "volumes=" << volumes << '\n';
    for( unsigned i = 0; i < seed.span(); ++i ){
      f << "subsetseed=";
      seed.writePosition( f, i );
      f << '\n';
    }
  }

  if( !f ) ERR( "can't write file: " + fileName );
}

// Write a per-volume .prj file, with info about a database volume
void writeInnerPrj( const std::string& fileName,
		    const MultiSequence& multi, const SubsetSuffixArray& sa ){
  std::ofstream f( fileName.c_str() );
  f << "totallength=" << multi.finishedSize() << '\n';
  f << "specialcharacters=" << multi.finishedSize() - sa.indexSize() << '\n';
  f << "numofsequences=" << multi.finishedSequences() << '\n';
  f << "prefixlength=" << sa.maxBucketPrefix() << '\n';
  if( !f ) ERR( "can't write file: " + fileName );
}

// Make one database volume, from one batch of sequences
void makeVolume( SubsetSuffixArray& sa, const MultiSequence& multi,
		 const LastdbArguments& args, const CyclicSubsetSeed& seed,
		 unsigned volumeNumber ){
  if( args.isCountsOnly ) return;

  std::string baseName = args.lastdbName + stringify(volumeNumber);

  LOG( "sorting..." );
  sa.sortIndex( multi.seqReader(), seed );

  LOG( "bucketing..." );
  sa.makeBuckets( multi.seqReader(), seed, args.bucketDepth );

  LOG( "writing..." );
  writeInnerPrj( baseName + ".prj", multi, sa );
  multi.toFiles( baseName );
  sa.toFiles( baseName );

  LOG( "done!" );
}

// Read the next sequence, adding it to the MultiSequence and the SuffixArray
std::istream&
appendFromFasta( MultiSequence& multi, SubsetSuffixArray& sa,
		 const LastdbArguments& args, const Alphabet& alph,
		 const CyclicSubsetSeed& seed, std::istream& in ){
  indexT maxSeqLen = args.volumeSize - sa.indexBytes();
  if( maxSeqLen < args.volumeSize - sa.indexBytes() ) maxSeqLen = indexT(-1);
  if( args.volumeSize < sa.indexBytes() ) maxSeqLen = 0;
  if( multi.finishedSequences() == 0 ) maxSeqLen = indexT(-1);

  indexT oldUnfinishedSize = multi.unfinishedSize();
  indexT oldFinishedSize = multi.finishedSize();

  multi.appendFromFasta( in, maxSeqLen );

  if( !multi.isFinished() && multi.finishedSequences() == 0 )
    ERR( "encountered a sequence that's too long" );

  // encode the newly-read sequence
  alph.tr( multi.seqWriter() + oldUnfinishedSize,
           multi.seqWriter() + multi.unfinishedSize() );

  if( in && multi.isFinished() && !args.isCountsOnly ){
    std::size_t maxIndexBytes = args.volumeSize - multi.unfinishedSize();
    if( args.volumeSize < multi.unfinishedSize() ) maxIndexBytes = 0;
    if( multi.finishedSequences() == 1 ) maxIndexBytes = std::size_t(-1);

    if( !sa.addIndices( multi.seqReader(),
                        oldFinishedSize, multi.finishedSize(),
                        args.indexStep, seed, maxIndexBytes ) ){
      multi.unfinish();
    }
  }

  return in;
}

void lastdb( int argc, char** argv ){
  LastdbArguments args;
  args.fromArgs( argc, argv );
  Alphabet alph;
  CyclicSubsetSeed seed;
  MultiSequence multi;
  SubsetSuffixArray sa;
  makeAlphabet( alph, args );
  makeSubsetSeed( seed, args, alph );
  multi.initForAppending(1);
  alph.tr( multi.seqWriter(), multi.seqWriter() + multi.unfinishedSize() );
  unsigned volumeNumber = 0;
  countT sequenceCount = 0;
  std::vector<countT> letterCounts( alph.size );

  for( char** i = argv + args.inputStart; i < argv + argc; ++i ){
    std::ifstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );
    LOG( "reading " << *i << "..." );

    while( appendFromFasta( multi, sa, args, alph, seed, in ) ){
      if( !args.isProtein && args.userAlphabet.empty() &&
          sequenceCount == 0 && isDubiousDna( alph, multi ) ){
        std::cerr << "lastdb: that's some funny-lookin DNA\n";
      }

      if( multi.isFinished() ){
        ++sequenceCount;
        indexT lastSeq = multi.finishedSequences() - 1;
        alph.count( multi.seqReader() + multi.seqBeg(lastSeq),
                    multi.seqReader() + multi.seqEnd(lastSeq),
                    &letterCounts[0] );
      }
      else{
	makeVolume( sa, multi, args, seed, volumeNumber++ );
	sa.clear();
	multi.reinitForAppending();
      }
    }
  }

  if( multi.finishedSequences() > 0 ){
    makeVolume( sa, multi, args, seed, volumeNumber++ );
  }

  writeOuterPrj( args.lastdbName + ".prj", args, alph,
                 sequenceCount, letterCounts, seed, volumeNumber );
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
