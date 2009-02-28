// Copyright 2008, 2009 Martin C. Frith

// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "LastalArguments.hh"
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

namespace cbrc{

typedef unsigned indexT;
typedef unsigned char uchar;
typedef unsigned long long countT;

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( ScoreMatrix& sm, LastalArguments& args,
		      const Alphabet& alph ){
  if( !args.matrixFile.empty() ){
    sm.fromFile( args.matrixFile );
  }
  else if( args.matchScore < 0 && args.mismatchCost < 0 &&
	   alph.letters == alph.protein ){
    sm.fromString( sm.blosum62 );
  }
  else{
    // if matchScore or mismatchCost weren't set, use default values:
    if( args.matchScore < 0 )    args.matchScore = 1;
    if( args.mismatchCost < 0 )  args.mismatchCost = 1;
    sm.matchMismatch( args.matchScore, args.mismatchCost, alph.letters );
  }

  sm.init( alph.encode );
}

// Read the .prj file for the whole database
void readOuterPrj( const std::string& fileName, Alphabet& alph,
		   PeriodicSpacedSeed& mask, unsigned& volumes ){
  std::ifstream f( fileName.c_str() );
  if( !f ) throw std::runtime_error("can't open file: " + fileName );

  std::string word;
  while( std::getline( f >> std::ws, word, '=' ) ){  // eat leading whitespace
    /**/ if( word == "alphabet" ) f >> alph;
    else if( word == "spacedseed" ) f >> mask;
    else if( word == "volumes" ) f >> volumes;
    else f >> word;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( alph.letters.empty() || mask.pattern.empty() || volumes == -1u ){
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
void writeCounts( const std::vector< std::vector<countT> >& matchCounts,
		  const MultiSequence& query, indexT minDepth,
		  std::ostream& out ){
  for( indexT i = 0; i < matchCounts.size(); ++i ){
    out << query.seqName(i) << '\n';

    for( indexT j = minDepth-1; j < matchCounts[i].size(); ++j ){
      out << j+1 << '\t' << matchCounts[i][j] << '\n';
    }

    out << '\n';  // blank line afterwards
  }
}

// Count all matches, of all sizes, of a query batch against a suffix array
void countMatches( std::vector< std::vector<countT> >& matchCounts,
		   const MultiSequence& query, const SuffixArray& sa,
		   const LastalArguments& args, char strand ){
  indexT seqNum = strand == '+' ? 0 : query.finishedSequences() - 1;

  for( indexT i = 0; i < query.ends.back(); i += args.queryStep ){
    if( strand == '+' ){
      for( ;; ){
	if( seqNum == query.finishedSequences() ) return;
	if( query.seqEnd(seqNum) > i ) break;
	++seqNum;
      }
    }
    else{
      for( ;; ){
	if( seqNum == -1u ) return;
	if( query.seqBeg(seqNum) < query.ends.back() - i ) break;
	--seqNum;
      }
    }

    sa.countMatches( matchCounts[seqNum], &query.seq[i] );
  }
}

// Find query matches to the suffix array, and do gapless extensions
void alignGapless( SegmentPairPot& gaplessAlns,
		   const MultiSequence& query, const MultiSequence& text,
		   const SuffixArray& sa,
		   const LastalArguments& args, const Alphabet& alph,
		   const int matGapless[64][64], const int matGapped[64][64],
		   char strand, std::ostream& out ){
  const uchar* qseq = &query.seq[0];
  const uchar* tseq = &text.seq[0];
  DiagonalTable dt;  // record already-covered positions on each diagonal
  countT gaplessExtensionCount = 0, gaplessAlignmentCount = 0;

  for( indexT i = 0; i < query.ends.back(); i += args.queryStep ){
    const indexT* beg;
    const indexT* end;
    sa.match( beg, end, qseq + i, args.oneHitMultiplicity, args.minHitDepth );

    for( /* noop */; beg < end; ++beg ){  // loop over suffix-array matches
      if( dt.isCovered( i, *beg ) ) continue;

      // tried: first get the score only, not the endpoints: slower!
      SegmentPair sp;  // do the gapless extension:
      sp.makeXdrop( *beg, i, tseq, qseq, matGapless, args.maxDropGapless );
      ++gaplessExtensionCount;

      // Tried checking the score after isOptimal & addEndpoint, but
      // the number of extensions decreased by < 10%, and it was
      // slower overall.
      if( sp.score < args.minScoreGapless ) continue;

      if( !sp.isOptimal( tseq, qseq, matGapless, args.maxDropGapless ) ){
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
	sp.makeXdrop( *beg, i, tseq, qseq, matGapped, args.maxDropGapped );
	if( !sp.isOptimal( tseq, qseq, matGapped, args.maxDropGapped ) ){
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
		  const std::vector<uchar>& text,
		  const std::vector<uchar>& query,
		  const LastalArguments& args, const Alphabet& alph,
		  const int mat[64][64], int matMax,
		  const GeneralizedAffineGapCosts& gap,
		  XdropAligner& aligner, Centroid& centroid ){
  countT gappedExtensionCount = 0;

  gaplessAlns.sort();  // sort the gapless alignments by score, highest first

  for( std::size_t i = 0; i < gaplessAlns.size(); ++i ){
    const SegmentPair& sp = gaplessAlns.get(i);

    if( sp.score == 0 ) continue;  // it has been marked as redundant

    Alignment aln;
    aln.seed = sp;

    // Shrink the seed to its longest run of identical matches.  This
    // trims off possibly unreliable parts of the gapless alignment.
    aln.seed.maxIdenticalRun( &text[0], &query[0], alph.canonical, mat );

    // do gapped extension from each end of the seed:
    aln.makeXdrop( aligner, centroid, &text[0], &query[0],
		   mat, matMax, gap, args.maxDropGapped );
    ++gappedExtensionCount;

    if( aln.score < args.minScoreGapped ) continue;

    if( !aln.isOptimal( &text[0], &query[0], mat, args.maxDropGapped, gap ) ){
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
		  const MultiSequence& query, const MultiSequence& text,
		  const LastalArguments& args, const Alphabet& alph,
		  const int mat[64][64], int matMax,
		  const GeneralizedAffineGapCosts& gap,
		  XdropAligner& aligner, Centroid& centroid,
		  char strand, std::ostream& out ){
  for( std::size_t i = 0; i < gappedAlns.size(); ++i ){
    const Alignment& aln = gappedAlns.items[i];
    if( args.outputType < 4 ){
      aln.write( text, query, strand, alph, args.outputFormat, out );
    }
    else{  // calculate match probabilities:
      Alignment probAln;
      probAln.seed = aln.seed;
      probAln.makeXdrop( aligner, centroid, &text.seq[0], &query.seq[0],
			 mat, matMax, gap, args.maxDropGapped,
			 args.gamma, args.outputType );
      probAln.write( text, query, strand, alph, args.outputFormat, out );
    }
  }
}

// Scan one batch of query sequences against one database volume
void scan( const MultiSequence& query, const MultiSequence& text,
	   const SuffixArray& sa, const LastalArguments& args,
	   const Alphabet& alph, const ScoreMatrix& sm,
	   std::vector< std::vector<countT> >& matchCounts,
	   char strand, std::ostream& out ){
  LOG( "scanning..." );

  if( args.outputType == 0 ){  // we just want match counts
    countMatches( matchCounts, query, sa, args, strand );
    return;
  }

  const int (*matGapless)[64] =  // score matrix for gapless alignment
    args.maskLowercase < 2 ? sm.caseInsensitive : sm.caseSensitive;

  const int (*matGapped)[64] =   // score matrix for gapped alignment
    args.maskLowercase < 3 ? sm.caseInsensitive : sm.caseSensitive;

  SegmentPairPot gaplessAlns;
  alignGapless( gaplessAlns, query, text, sa, args, alph,
		matGapless, matGapped, strand, out );
  if( args.outputType == 1 ) return;  // we just want gapless alignments

  XdropAligner aligner;
  Centroid centroid( aligner, matGapped, args.temperature );  // slow?
  GeneralizedAffineGapCosts gap( args.gapExistCost, args.gapExtendCost,
				 args.gapPairCost );

  AlignmentPot gappedAlns;
  alignGapped( gappedAlns, gaplessAlns, text.seq, query.seq, args, alph,
	       matGapped, sm.maxScore, gap, aligner, centroid );

  if( args.outputType > 2 ){  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG( "nonredundant gapped alignments=" << gappedAlns.size() );
  }

  gappedAlns.sort();  // sort by score
  alignFinish( gappedAlns, query, text, args, alph,
	       matGapped, sm.maxScore, gap, aligner, centroid, strand, out );
}

// Scan one batch of query sequences against all database volumes
void scanAllVolumes( MultiSequence& query,
		     const LastalArguments& args, const Alphabet& alph,
		     const PeriodicSpacedSeed& mask, const ScoreMatrix& sm,
		     unsigned volumes, std::ostream& out ){
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  if( args.outputType == 0 ) matchCounts.resize( query.finishedSequences() );

  for( unsigned i = 0; i < volumes; ++i ){
    std::string baseName = args.lastdbName + stringify(i);
    LOG( "reading " << baseName << "..." );

    indexT seqCount = -1u, delimiterNum = -1u, bucketDepth = -1u;
    readInnerPrj( baseName + ".prj", seqCount, delimiterNum, bucketDepth );

    MultiSequence text;
    text.fromFiles( baseName, seqCount );

    SuffixArray sa( text.seq, mask.offsets, alph.size );
    sa.fromFiles( baseName, text.ends.back() - delimiterNum, bucketDepth );

    if( args.strand == 2 && i > 0 ){
      LOG( "re-reverse complementing..." );
      alph.rc( query.seq.begin(), query.seq.begin() + query.ends.back() );
    }

    if( args.strand != 0 ){
      scan( query, text, sa, args, alph, sm, matchCounts, '+', out );
    }

    if( args.strand == 2 || (args.strand == 0 && i == 0) ){
      LOG( "reverse complementing..." );
      alph.rc( query.seq.begin(), query.seq.begin() + query.ends.back() );
    }

    if( args.strand != 1 ){
      scan( query, text, sa, args, alph, sm, matchCounts, '-', out );
    }
  }

  if( args.outputType == 0 ){
    LOG( "writing..." );
    writeCounts( matchCounts, query, args.minHitDepth, out );
  }

  LOG( "done!" );
}

void writeHeader( std::ostream& out, const LastalArguments& args,
		  const ScoreMatrix sm ){
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
    sm.writeCommented( out );
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
std::istream&
appendFromFasta( MultiSequence& query,
		 const LastalArguments& args, const Alphabet& alph,
		 std::istream& in ){
  std::size_t maxSeqBytes = args.batchSize;
  if( query.finishedSequences() == 0 ) maxSeqBytes = std::size_t(-1);

  indexT oldSeqSize = query.seq.size();

  query.appendFromFasta( in, maxSeqBytes );

  // encode the newly-read sequence
  alph.tr( query.seq.begin() + oldSeqSize, query.seq.end() );

  return in;
}

void lastal( int argc, char** argv ){
  std::clock_t startTime = std::clock();

  LastalArguments args;
  args.fromArgs( argc, argv );

  Alphabet alph;
  PeriodicSpacedSeed mask;
  unsigned volumes = -1u;  // initialize it to an "error" value
  readOuterPrj( args.lastdbName + ".prj", alph, mask, volumes );

  ScoreMatrix sm;
  makeScoreMatrix( sm, args, alph );  // before alph.makeCaseInsensitive!
  args.setDefaults( alph.letters == alph.dna, alph.letters == alph.protein,
		    sm.maxScore );
  if( args.maskLowercase < 1 ) alph.makeCaseInsensitive();

  if( args.temperature == -1 && args.outputType > 3 ){
    args.temperature =
      1 / LambdaCalculator::calculate( sm.caseSensitive, alph.size );
    if( args.temperature < 0 )
      throw std::runtime_error("can't calculate lambda for this score matrix");
  }

  std::ofstream outFileStream;
  std::ostream& out = openOut( args.outFile, outFileStream );
  writeHeader( out, args, sm );
  MultiSequence query;
  query.initForAppending( mask.maxOffset );
  alph.tr( query.seq.begin(), query.seq.end() );
  out.precision(3);  // print non-integers more compactly

  for( char** i = argv + args.inputStart; i < argv + argc; ++i ){
    LOG( "reading " << *i << "..." );
    std::ifstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );

    while( appendFromFasta( query, args, alph, in ) ){
      if( !query.isFinished() ){
	scanAllVolumes( query, args, alph, mask, sm, volumes, out );
	query.reinitForAppending();
      }
    }
  }

  if( query.finishedSequences() > 0 ){
    scanAllVolumes( query, args, alph, mask, sm, volumes, out );
  }

  out.precision(6);  // reset the precision to the default value
  out << "# CPU time: " << (clock() - startTime + 0.0) / CLOCKS_PER_SEC
      << " seconds\n";
}

}

int main( int argc, char** argv )
try{
  cbrc::lastal( argc, argv );
  return EXIT_SUCCESS;
}
catch( const std::exception& e ) {
  std::cerr << "lastal: " << e.what() << '\n';
  return EXIT_FAILURE;
}
