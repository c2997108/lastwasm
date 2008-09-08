// Copyright 2008 Martin C. Frith

// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "LastalArguments.hh"
#include "SuffixArray.hh"
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
#include <cassert>

#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastal: " << x << '\n'

typedef unsigned indexT;
typedef unsigned char uchar;
typedef unsigned long long countT;

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( cbrc::ScoreMatrix& sm, cbrc::LastalArguments& args,
		      const cbrc::Alphabet& alph ){
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
void readOuterPrj( const std::string& fileName, cbrc::Alphabet& alph,
		   cbrc::PeriodicSpacedSeed& mask, unsigned& volumes ){
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
		  const cbrc::MultiSequence& query, indexT minDepth,
		  std::ostream& out ){
  for( indexT i = 0; i < matchCounts.size(); ++i ){
    out << query.names[i] << '\n';

    for( indexT j = minDepth-1; j < matchCounts[i].size(); ++j ){
      out << j+1 << '\t' << matchCounts[i][j] << '\n';
    }
  }
}

// Count all matches, of all sizes, of a query batch against a suffix array
void countMatches( std::vector< std::vector<countT> >& matchCounts,
		   const cbrc::MultiSequence& query,
		   const cbrc::SuffixArray& sa,
		   const cbrc::LastalArguments& args, char strand ){
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
void alignGapless( cbrc::SegmentPairPot& gaplessAlns,
		   const cbrc::MultiSequence& query,
		   const cbrc::MultiSequence& text,
		   const cbrc::SuffixArray& sa,
		   const cbrc::LastalArguments& args,
		   const cbrc::ScoreMatrix& sm,
		   const cbrc::Alphabet& alph,
		   char strand, std::ostream& out ){
  const std::vector<uchar>& qseq = query.seq;
  const std::vector<uchar>& tseq = text.seq;
  const int (*mat)[64] =
    args.maskLowercase < 2 ? sm.caseInsensitive : sm.caseSensitive;
  cbrc::DiagonalTable dt;  // record already-covered positions on each diagonal
  countT gaplessExtensionCount = 0, gaplessAlignmentCount = 0;

  for( indexT i = 0; i < query.ends.back(); i += args.queryStep ){
    const indexT* beg;
    const indexT* end;
    sa.match( beg, end, &qseq[i], args.oneHitMultiplicity, args.minHitDepth );

    for( /* noop */; beg < end; ++beg ){  // loop over suffix-array matches
      if( dt.isCovered( i, *beg ) ) continue;

      // tried: first get the score only, not the endpoints: slower!
      cbrc::SegmentPair sp;  // do the gapless extension:
      sp.makeXdrop( *beg, i, &tseq[0], &qseq[0], mat, args.maxDropGapless );
      ++gaplessExtensionCount;

      // Tried checking the score after isOptimal & addEndpoint, but
      // the number of extensions decreased by < 10%, and it was
      // slower overall.
      if( sp.score < args.minScoreGapless ) continue;

      if( !sp.isOptimal( &tseq[0], &qseq[0], mat, args.maxDropGapless ) ){
	continue;  // ignore sucky gapless extensions
      }

      dt.addEndpoint( sp.end2(), sp.end1() );

      if( args.outputType == 1 ){  // we just want gapless alignments
	cbrc::Alignment aln;
	aln.fromSegmentPair(sp);
	aln.write( text, query, strand, alph, args.outputFormat, out );
      }
      else gaplessAlns.add(sp);  // add the gapless alignment to the pot
      ++gaplessAlignmentCount;
    }
  }

  LOG( "gapless extensions=" << gaplessExtensionCount );
  LOG( "gapless alignments=" << gaplessAlignmentCount );
}

// Do gapped extensions of the gapless alignments
void alignGapped( cbrc::AlignmentPot& gappedAlns,
		  cbrc::SegmentPairPot& gaplessAlns,
		  const std::vector<uchar>& text,
		  const std::vector<uchar>& query,
		  const cbrc::LastalArguments& args,
		  const cbrc::ScoreMatrix& sm,
		  const cbrc::Alphabet& alph ){
  cbrc::XdropAligner aligner;
  cbrc::GeneralizedAffineGapCosts gap( args.gapExistCost, args.gapExtendCost,
				       args.gapPairCost );
  const int (*mat)[64] =
    args.maskLowercase < 3 ? sm.caseInsensitive : sm.caseSensitive;
  countT gappedExtensionCount = 0;

  gaplessAlns.sort();  // sort the gapless alignments by score, highest first

  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    const cbrc::SegmentPair& sp = gaplessAlns.get(i);

    if( sp.score == 0 ) continue;  // it has been marked as redundant

    cbrc::SegmentPair seed;

    // Redo gapless extension, using gapped score parameters.  Without
    // this, if we self-compare a huge sequence, we risk getting a
    // huge gapped extension.  (Do this before sorting them?)
    seed.makeXdrop( sp.beg1(), sp.beg2(), &text[0], &query[0],
		    mat, args.maxDropGapped );

    // Shrink the seed to its longest run of identical matches.  This
    // trims off possibly unreliable parts of the gapless alignment.
    seed.maxIdenticalRun( &text[0], &query[0], alph.canonical, mat );

    cbrc::Alignment aln;  // do gapped extension from each end of the seed:
    aln.makeXdrop( aligner, seed, &text[0], &query[0],
		   mat, args.maxDropGapped, gap );
    ++gappedExtensionCount;

    if( aln.score < args.minScoreGapped ) continue;

    if( !aln.isOptimal( &text[0], &query[0], mat, args.maxDropGapped, gap ) ){
      // If retained, non-"optimal" alignments can hide "optimal"
      // alignments, e.g. during non-reduntantization.
      continue;
    }

    gaplessAlns.markAllOverlaps( aln.blocks );
    gaplessAlns.markTandemRepeats( seed, args.maxRepeatDistance ); 
    gappedAlns.add(aln);  // add the gapped alignment to the pot
 }

  LOG( "gapped extensions=" << gappedExtensionCount );
  LOG( "gapped alignments=" << gappedAlns.size() );
}

// Scan one batch of query sequences against one database volume
void scan( const cbrc::MultiSequence& query, const cbrc::MultiSequence& text,
	   const cbrc::SuffixArray& sa, const cbrc::LastalArguments& args,
	   const cbrc::Alphabet& alph, const cbrc::ScoreMatrix& sm,
	   std::vector< std::vector<countT> >& matchCounts,
	   char strand, std::ostream& out ){
  LOG( "scanning..." );

  if( args.outputType == 0 ){  // we just want match counts
    countMatches( matchCounts, query, sa, args, strand );
    return;
  }

  cbrc::SegmentPairPot gaplessAlns;
  alignGapless( gaplessAlns, query, text, sa, args, sm, alph, strand, out );
  if( args.outputType == 1 ) return;  // we just want gapless alignments

  cbrc::AlignmentPot gappedAlns;
  alignGapped( gappedAlns, gaplessAlns, text.seq, query.seq, args, sm, alph );

  if( args.outputType == 3 ){  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG( "nonredundant gapped alignments=" << gappedAlns.size() );
  }

  gappedAlns.sort();  // sort by score
  for( size_t i = 0; i < gappedAlns.size(); ++i ){
    gappedAlns.items[i].write( text, query, strand, alph,
			       args.outputFormat, out );
  }
}

// Scan one batch of query sequences against all database volumes
void scanAllVolumes( cbrc::MultiSequence& query,
		     const cbrc::LastalArguments& args,
		     const cbrc::Alphabet& alph,
		     const cbrc::PeriodicSpacedSeed& mask,
		     const cbrc::ScoreMatrix& sm,
		     unsigned volumes, std::ostream& out ){
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  if( args.outputType == 0 ) matchCounts.resize( query.finishedSequences() );

  for( unsigned i = 0; i < volumes; ++i ){
    std::string baseName = args.lastdbName + cbrc::stringify(i);
    LOG( "reading " << baseName << "..." );

    indexT seqCount = -1u, delimiterNum = -1u, bucketDepth = -1u;
    readInnerPrj( baseName + ".prj", seqCount, delimiterNum, bucketDepth );

    cbrc::MultiSequence text;
    text.fromFiles( baseName, seqCount );

    cbrc::SuffixArray sa( text.seq, mask.offsets, alph.size );
    sa.fromFiles( baseName, text.ends.back() - delimiterNum, bucketDepth );

    if( args.strand == 2 && i > 0 ){
      LOG( "re-reverse complementing..." );
      alph.rc( query.seq.begin(), query.seq.begin() + query.ends.back() );
    }

    if( args.strand != 0 ){
      scan( query, text, sa, args, alph, sm, matchCounts, '+', out );
    }

    if( args.strand == 2 || args.strand == 0 && i == 0 ){
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

void writeHeader( std::ostream& out, const cbrc::LastalArguments& args,
		  const cbrc::ScoreMatrix sm ){
  out << "# LAST\n";
  out << "# Version " <<
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
appendFromFasta( cbrc::MultiSequence& query,
		 const cbrc::LastalArguments& args, const cbrc::Alphabet& alph,
		 std::istream& in ){
  std::size_t maxSeqBytes = args.batchSize;
  if( query.finishedSequences() == 0 ) maxSeqBytes = std::size_t(-1);

  indexT oldSeqSize = query.seq.size();

  query.appendFromFasta( in, maxSeqBytes );

  // encode the newly-read sequence
  alph.tr( query.seq.begin() + oldSeqSize, query.seq.end() );

  return in;
}

int main( int argc, char** argv )
try{
  std::clock_t startTime = std::clock();
  cbrc::LastalArguments args( argc, argv );
  cbrc::Alphabet alph;
  cbrc::PeriodicSpacedSeed mask;
  cbrc::ScoreMatrix sm;
  unsigned volumes = -1u;  // initialize it to an "error" value
  readOuterPrj( args.lastdbName + ".prj", alph, mask, volumes );
  makeScoreMatrix( sm, args, alph );  // before alph.makeCaseInsensitive!
  args.setDefaults( alph.letters == alph.dna, alph.letters == alph.protein,
		    sm.maxScore );
  if( args.maskLowercase < 1 ) alph.makeCaseInsensitive();

  std::ofstream outFileStream;
  std::ostream& out = cbrc::openOut( args.outFile, outFileStream );
  writeHeader( out, args, sm );
  cbrc::MultiSequence query;
  query.initForAppending( mask.maxOffset );
  alph.tr( query.seq.begin(), query.seq.end() );

  for( char** i = argv + args.inputStart; i < argv + argc; ++i ){
    LOG( "reading " << *i << "..." );
    std::ifstream inFileStream;
    std::istream& in = cbrc::openIn( *i, inFileStream );

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

  out << "# CPU time: " << (clock() - startTime + 0.0) / CLOCKS_PER_SEC
      << " seconds\n";

  return EXIT_SUCCESS;
}
catch( const std::exception& e ) {
  std::cerr << "lastal: " << e.what() << '\n';
  return EXIT_FAILURE;
}
