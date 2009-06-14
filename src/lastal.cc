// Copyright 2008, 2009 Martin C. Frith

// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "LastalArguments.hh"
#include "QualityScoreCalculator.hh"
#include "LambdaCalculator.hh"
#include "SuffixArray.hh"
#include "Centroid.hh"
#include "XdropAligner.hh"
#include "AlignmentPot.hh"
#include "Alignment.hh"
#include "SegmentPairPot.hh"
#include "SegmentPair.hh"
#include "ScoreMatrix.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "PeriodicSpacedSeed.hh"
#include "DiagonalTable.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "io.hh"
#include "stringify.hh"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <cassert>

#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastal: " << x << '\n'

using namespace cbrc;

namespace {
  typedef unsigned indexT;
  typedef unsigned char uchar;
  typedef unsigned long long countT;

  LastalArguments args;
  Alphabet alph;
  PeriodicSpacedSeed spacedSeed;
  ScoreMatrix scoreMatrix;
  GeneralizedAffineGapCosts gapCosts;
  XdropAligner xdropAligner;
  MultiSequence query;  // sequence that hasn't been indexed by lastdb
  MultiSequence text;  // sequence that has been indexed by lastdb
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  enum { MAT = 64 };
  const int (*matGapless)[MAT];  // score matrix for gapless alignment
  const int (*matGapped)[MAT];  // score matrix for gapped alignment
  QualityScoreCalculator qualityScoreCalculator;
  int (*pssm)[MAT];  // Position Specific Score Matrix
}

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( const std::string& matrixFile ){
  if( !matrixFile.empty() ){
    scoreMatrix.fromString( matrixFile );
  }
  else if( args.matchScore < 0 && args.mismatchCost < 0 &&
	   alph.letters == alph.protein ){
    scoreMatrix.fromString( scoreMatrix.blosum62 );
  }
  else{
    scoreMatrix.matchMismatch( args.matchScore, args.mismatchCost,
			       alph.letters );
  }

  scoreMatrix.init( alph.encode );

  matGapless = args.maskLowercase < 2 ?
    scoreMatrix.caseInsensitive : scoreMatrix.caseSensitive;

  matGapped = args.maskLowercase < 3 ?
    scoreMatrix.caseInsensitive : scoreMatrix.caseSensitive;
}

// Read the .prj file for the whole database
void readOuterPrj( const std::string& fileName, unsigned& volumes ){
  std::ifstream f( fileName.c_str() );
  if( !f ) throw std::runtime_error("can't open file: " + fileName );

  std::string word;
  while( std::getline( f >> std::ws, word, '=' ) ){  // eat leading whitespace
    /**/ if( word == "alphabet" ) f >> alph;
    else if( word == "spacedseed" ) f >> spacedSeed;
    else if( word == "volumes" ) f >> volumes;
    else f >> word;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( alph.letters.empty() || spacedSeed.pattern.empty() || volumes == -1u ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) throw std::runtime_error("can't read file: " + fileName);
}

// Read a per-volume .prj file, with info about a database volume
void readInnerPrj( const std::string& fileName, indexT& seqCount,
		   indexT& delimiterNum, indexT& bucketDepth ){
  std::ifstream f( fileName.c_str() );
  if( !f ) throw std::runtime_error("can't open file: " + fileName );

  std::string word;
  while( std::getline( f >> std::ws, word, '=' ) ){  // eat leading whitespace
    /**/ if( word == "numofsequences" ) f >> seqCount;
    else if( word == "specialcharacters" ) f >> delimiterNum;
    else if( word == "prefixlength" ) f >> bucketDepth;
    else f >> word;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( seqCount == -1u || delimiterNum == -1u || bucketDepth == -1u ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) throw std::runtime_error("can't read file: " + fileName);
}

// Write match counts for each query sequence
void writeCounts( std::ostream& out ){
  LOG( "writing..." );

  for( indexT i = 0; i < matchCounts.size(); ++i ){
    out << query.seqName(i) << '\n';

    for( indexT j = args.minHitDepth - 1; j < matchCounts[i].size(); ++j ){
      out << j+1 << '\t' << matchCounts[i][j] << '\n';
    }

    out << '\n';  // blank line afterwards
  }
}

// Count all matches, of all sizes, of a query batch against a suffix array
void countMatches( const SuffixArray& suffixArray, char strand ){
  LOG( "counting..." );
  indexT seqNum = strand == '+' ? 0 : query.finishedSequences() - 1;

  for( indexT i = 0; i < query.ends.back(); i += args.queryStep ){
    if( strand == '+' ){
      for( ;; ){
	if( seqNum == query.finishedSequences() ) return;
	if( query.seqEnd(seqNum) > i ) break;
	++seqNum;
      }
      // speed-up (for spaced seeds it's not optimal):
      if( args.minHitDepth > query.seqEnd(seqNum) - i ) continue;
    }
    else{
      indexT j = query.ends.back() - i;
      for( ;; ){
	if( seqNum == -1u ) return;
	if( query.seqBeg(seqNum) < j ) break;
	--seqNum;
      }
      // speed-up (for spaced seeds it's not optimal):
      if( args.minHitDepth > j - query.seqBeg(seqNum) ) continue;
    }

    suffixArray.countMatches( matchCounts[seqNum], &query.seq[i] );
  }
}

// Find query matches to the suffix array, and do gapless extensions
void alignGapless( SegmentPairPot& gaplessAlns, const SuffixArray& suffixArray,
		   char strand, std::ostream& out ){
  const uchar* qseq = &query.seq[0];
  const uchar* tseq = &text.seq[0];
  DiagonalTable dt;  // record already-covered positions on each diagonal
  countT gaplessExtensionCount = 0, gaplessAlignmentCount = 0;

  for( indexT i = 0; i < query.ends.back(); i += args.queryStep ){
    const indexT* beg;
    const indexT* end;
    suffixArray.match( beg, end, qseq + i,
		       args.oneHitMultiplicity, args.minHitDepth );

    // Tried: if we hit a delimiter when using contiguous seeds, then
    // increase "i" to the delimiter position.  This gave a speed-up
    // of only 3%, with 34-nt tags.

    for( /* noop */; beg < end; ++beg ){  // loop over suffix-array matches
      if( dt.isCovered( i, *beg ) ) continue;

      // tried: first get the score only, not the endpoints: slower!
      SegmentPair sp;  // do the gapless extension:
      sp.makeXdrop( *beg, i, tseq, qseq,
		    matGapless, args.maxDropGapless, pssm );
      ++gaplessExtensionCount;

      // Tried checking the score after isOptimal & addEndpoint, but
      // the number of extensions decreased by < 10%, and it was
      // slower overall.
      if( sp.score < args.minScoreGapless ) continue;

      if( !sp.isOptimal( tseq, qseq, matGapless, args.maxDropGapless, pssm ) ){
	continue;  // ignore sucky gapless extensions
      }

      if( args.outputType == 1 ){  // we just want gapless alignments
	Alignment aln;
	aln.fromSegmentPair(sp);
	aln.write( text, query, strand, alph, args.outputFormat, out );
      }
      else{
	// Redo gapless extension, using gapped score parameters.  Without
	// this, if we self-compare a huge sequence, we risk getting a
	// huge gapped extension.
	sp.makeXdrop( *beg, i, tseq, qseq,
		      matGapped, args.maxDropGapped, pssm );
	if( !sp.isOptimal( tseq, qseq, matGapped, args.maxDropGapped, pssm ) ){
	  continue;
	}
	gaplessAlns.add(sp);  // add the gapless alignment to the pot
      }

      ++gaplessAlignmentCount;
      dt.addEndpoint( sp.end2(), sp.end1() );
    }
  }

  LOG( "gapless extensions=" << gaplessExtensionCount );
  LOG( "gapless alignments=" << gaplessAlignmentCount );
}

// Do gapped extensions of the gapless alignments
void alignGapped( AlignmentPot& gappedAlns, SegmentPairPot& gaplessAlns,
		  Centroid& centroid ){
  const uchar* qseq = &query.seq[0];
  const uchar* tseq = &text.seq[0];
  countT gappedExtensionCount = 0;

  gaplessAlns.sort();  // sort the gapless alignments by score, highest first

  for( std::size_t i = 0; i < gaplessAlns.size(); ++i ){
    const SegmentPair& sp = gaplessAlns.get(i);

    if( sp.score == 0 ) continue;  // it has been marked as redundant

    Alignment aln;
    aln.seed = sp;

    // Shrink the seed to its longest run of identical matches.  This
    // trims off possibly unreliable parts of the gapless alignment.
    aln.seed.maxIdenticalRun( tseq, qseq, alph.canonical, matGapped, pssm );

    // do gapped extension from each end of the seed:
    aln.makeXdrop( xdropAligner, centroid, tseq, qseq, matGapped,
		   scoreMatrix.maxScore, gapCosts, args.maxDropGapped, pssm );
    ++gappedExtensionCount;

    if( aln.score < args.minScoreGapped ) continue;

    if( !aln.isOptimal( tseq, qseq,
			matGapped, args.maxDropGapped, gapCosts, pssm ) ){
      // If retained, non-"optimal" alignments can hide "optimal"
      // alignments, e.g. during non-reduntantization.
      continue;
    }

    gaplessAlns.markAllOverlaps( aln.blocks );
    gaplessAlns.markTandemRepeats( aln.seed, args.maxRepeatDistance ); 
    gappedAlns.add(aln);  // add the gapped alignment to the pot
 }

  LOG( "gapped extensions=" << gappedExtensionCount );
  LOG( "gapped alignments=" << gappedAlns.size() );
}

// Print the gapped alignments, after optionally calculating match
// probabilities and re-aligning using the gamma-centroid algorithm
void alignFinish( const AlignmentPot& gappedAlns,
		  Centroid& centroid, char strand, std::ostream& out ){
  for( std::size_t i = 0; i < gappedAlns.size(); ++i ){
    const Alignment& aln = gappedAlns.items[i];
    if( args.outputType < 4 ){
      aln.write( text, query, strand, alph, args.outputFormat, out );
    }
    else{  // calculate match probabilities:
      Alignment probAln;
      probAln.seed = aln.seed;
      probAln.makeXdrop( xdropAligner, centroid, &text.seq[0], &query.seq[0],
			 matGapped, scoreMatrix.maxScore, gapCosts,
			 args.maxDropGapped, pssm, args.gamma, args.outputType );
      probAln.write( text, query, strand, alph, args.outputFormat, out );
    }
  }
}

// Scan one batch of query sequences against one database volume
void scan( const SuffixArray& suffixArray, char strand, std::ostream& out ){
  if( args.outputType == 0 ){  // we just want match counts
    countMatches( suffixArray, strand );
    return;
  }

  if( args.inputFormat > 0 ){
    LOG( "making PSSM..." );
    qualityScoreCalculator.makePssm( pssm, &query.qualityScores[0],
				     &query.seq[0], query.ends.back(),
				     args.inputFormat < 3 );
  }

  LOG( "scanning..." );

  SegmentPairPot gaplessAlns;
  alignGapless( gaplessAlns, suffixArray, strand, out );
  if( args.outputType == 1 ) return;  // we just want gapless alignments

  Centroid centroid( xdropAligner, matGapped, args.temperature );  // slow?

  AlignmentPot gappedAlns;
  alignGapped( gappedAlns, gaplessAlns, centroid );

  if( args.outputType > 2 ){  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG( "nonredundant gapped alignments=" << gappedAlns.size() );
  }

  gappedAlns.sort();  // sort by score
  alignFinish( gappedAlns, centroid, strand, out );
}

// Read one database volume
void readVolume( SuffixArray& suffixArray, unsigned volumeNumber ){
  std::string baseName = args.lastdbName + stringify(volumeNumber);
  LOG( "reading " << baseName << "..." );
  indexT seqCount = -1u, delimiterNum = -1u, bucketDepth = -1u;
  readInnerPrj( baseName + ".prj", seqCount, delimiterNum, bucketDepth );
  text.fromFiles( baseName, seqCount );
  suffixArray.fromFiles( baseName,
			 text.ends.back() - delimiterNum, bucketDepth );
}

void reverseComplementQuery(){
  LOG( "reverse complementing..." );
  alph.rc( query.seq.begin(), query.seq.begin() + query.ends.back() );
  if( args.inputFormat > 0 ){
    std::reverse( query.qualityScores.begin(),  // XXX ugly
		  query.qualityScores.begin() + query.ends.back() *
		  (query.qualityScores.size() / query.seq.size()) );
  }
}

// Scan one batch of query sequences against all database volumes
void scanAllVolumes( SuffixArray& suffixArray,
		     unsigned volumes, std::ostream& out ){
  if( args.outputType == 0 ) matchCounts.resize( query.finishedSequences() );
  else if( args.inputFormat > 0 ) pssm = new int[ query.ends.back() ][MAT];

  for( unsigned i = 0; i < volumes; ++i ){
    if( text.seq.empty() || volumes > 1 ) readVolume( suffixArray, i );

    if( args.strand == 2 && i > 0 ) reverseComplementQuery();

    if( args.strand != 0 ) scan( suffixArray, '+', out );

    if( args.strand == 2 || (args.strand == 0 && i == 0) )
      reverseComplementQuery();

    if( args.strand != 1 ) scan( suffixArray, '-', out );
  }

  if( args.outputType == 0 ) writeCounts( out );
  else if( args.inputFormat > 0 ) delete[] pssm;

  LOG( "query batch done!" );
}

void writeHeader( std::ostream& out ){
  out << "# LAST version " <<
#include "version.hh"
      << "\n";
  out << "#\n";
  args.writeCommented( out );
  out << "#\n";

  if( args.outputType == 0 ){  // we just want hit counts
    out << "# depth\tcount\n";
  }
  else{  // we want alignments
    scoreMatrix.writeCommented( out );
    out << "#\n";
    out << "# Coordinates are 0-based.  For - strand matches, coordinates\n";
    out << "# in the reverse complement of the 2nd sequence are used.\n";
    out << "#\n";

    if( args.outputFormat == 0 ){  // tabular format
      out << "# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\t"
	  << "name2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n";
    }
    else{  // MAF format
      out << "# name start alnSize strand seqSize alignment\n";
    }
  }

  out << "#\n";
}

// Read the next sequence, adding it to the MultiSequence
std::istream& appendFromFasta( std::istream& in ){
  std::size_t maxSeqBytes = args.batchSize;
  if( query.finishedSequences() == 0 ) maxSeqBytes = std::size_t(-1);

  indexT oldSeqSize = query.seq.size();

  /**/ if( args.inputFormat < 1 ) query.appendFromFasta( in, maxSeqBytes );
  else if( args.inputFormat < 3 ) query.appendFromFastq( in, maxSeqBytes );
  else query.appendFromPrb( in, maxSeqBytes, alph.size, alph.decode );

  // encode the newly-read sequence
  alph.tr( query.seq.begin() + oldSeqSize, query.seq.end() );

  return in;
}

void lastal( int argc, char** argv ){
  std::clock_t startTime = std::clock();

  args.fromArgs( argc, argv );

  std::string matrixFile;
  if( !args.matrixFile.empty() ){
    matrixFile = slurp( args.matrixFile );
    args.fromString( matrixFile );  // read options from the matrix file
    args.fromArgs( argc, argv );  // command line overrides matrix file
  }

  unsigned volumes = -1u;  // initialize it to an "error" value
  readOuterPrj( args.lastdbName + ".prj", volumes );

  args.setDefaultsFromAlphabet( alph.letters == alph.dna,
				alph.letters == alph.protein );
  makeScoreMatrix( matrixFile );  // before alph.makeCaseInsensitive
  if( args.maskLowercase < 1 ) alph.makeCaseInsensitive();

  double lambda = -1;
  if( args.temperature < 0 && (args.outputType > 3 || args.inputFormat > 0) ){
    // it makes no difference whether we use matGapped or matGapless here:
    lambda = LambdaCalculator::calculate( matGapped, alph.size );
    if( lambda < 0 )
      throw std::runtime_error("can't calculate lambda for this score matrix");
  }

  args.setDefaultsFromMatrix( scoreMatrix.maxScore, lambda );

  gapCosts.assign( args.gapExistCost, args.gapExtendCost, args.gapPairCost );

  SuffixArray suffixArray( text.seq, spacedSeed.offsets, alph.size );

  if( args.inputFormat > 0 ){
    assert( matGapless == matGapped );
    int asciiOffset = (args.inputFormat < 2) ? 33 : 64;
    bool isMatchMismatch = args.matrixFile.empty() && args.matchScore > 0;
    qualityScoreCalculator.init( matGapped, alph.size, args.temperature,
				 args.maskLowercase > 2, isMatchMismatch,
				 args.matchScore, -args.mismatchCost,
				 alph.canonical,
				 args.inputFormat < 2, asciiOffset );
  }

  std::ofstream outFileStream;
  std::ostream& out = openOut( args.outFile, outFileStream );
  writeHeader( out );
  out.precision(3);  // print non-integers more compactly

  query.initForAppending( spacedSeed.maxOffset );
  alph.tr( query.seq.begin(), query.seq.end() );

  for( char** i = argv + args.inputStart; i < argv + argc; ++i ){
    LOG( "reading " << *i << "..." );
    std::ifstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );

    while( appendFromFasta( in ) ){
      if( !query.isFinished() ){
	scanAllVolumes( suffixArray, volumes, out );
	query.reinitForAppending();
      }
    }
  }

  if( query.finishedSequences() > 0 ){
    scanAllVolumes( suffixArray, volumes, out );
  }

  out.precision(6);  // reset the precision to the default value
  out << "# CPU time: " << (clock() - startTime + 0.0) / CLOCKS_PER_SEC
      << " seconds\n";
}

int main( int argc, char** argv )
try{
  lastal( argc, argv );
  return EXIT_SUCCESS;
}
catch( const std::exception& e ) {
  std::cerr << "lastal: " << e.what() << '\n';
  return EXIT_FAILURE;
}
