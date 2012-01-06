// Copyright 2008, 2009, 2011, 2012 Martin C. Frith

#include "Alignment.hh"
#include "Centroid.hh"
#include "GappedXdropAligner.hh"
#include "GeneticCode.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "TwoQualityScoreMatrix.hh"
#include <cassert>

// make C++ tolerable:
#define IT(type) std::vector<type>::iterator
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

void Alignment::fromSegmentPair( const SegmentPair& sp ){
  blocks.assign( 1, sp );
  score = sp.score;
}

// Does x precede and touch y in both sequences?
static bool isNext( const SegmentPair& x, const SegmentPair& y ){
  return x.end1() == y.beg1() && x.end2() == y.beg2();
}

void Alignment::makeXdrop( GappedXdropAligner& aligner, Centroid& centroid,
			   const uchar* seq1, const uchar* seq2,
			   const int scoreMatrix[MAT][MAT], int smMax,
			   const GeneralizedAffineGapCosts& gap, int maxDrop,
			   int frameshiftCost, indexT frameSize,
			   const int pssm2[][MAT],
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2,
			   double gamma, int outputType ){
  score = seed.score;

  // extend a gapped alignment in the left/reverse direction from the seed:
  extend( blocks, columnAmbiguityCodes, aligner, centroid, seq1, seq2,
	  seed.beg1(), seed.beg2(), false,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, sm2qual, qual1, qual2, gamma, outputType );

  // convert left-extension coordinates to sequence coordinates:
  indexT seedBeg1 = seed.beg1();
  indexT seedBeg2 = aaToDna( seed.beg2(), frameSize );
  for( IT(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    i->start1 = seedBeg1 - i->start1 - i->size;
    i->start2 = dnaToAa( seedBeg2 - i->start2, frameSize ) - i->size;
  }

  // extend a gapped alignment in the right/forward direction from the seed:
  std::vector<SegmentPair> forwardBlocks;
  std::vector<uchar> forwardAmbiguities;
  extend( forwardBlocks, forwardAmbiguities, aligner, centroid, seq1, seq2,
	  seed.end1(), seed.end2(), true,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, sm2qual, qual1, qual2, gamma, outputType );

  // convert right-extension coordinates to sequence coordinates:
  indexT seedEnd1 = seed.end1();
  indexT seedEnd2 = aaToDna( seed.end2(), frameSize );
  for( IT(SegmentPair) i = forwardBlocks.begin(); i < forwardBlocks.end();
       ++i ){
    i->start1 = seedEnd1 + i->start1;
    i->start2 = dnaToAa( seedEnd2 + i->start2, frameSize );
  }

  // check that the seed isn't very bizarre and dubious:
  assert( (seed.size > 0) ||
	  (!blocks.empty() && isNext( blocks.back(), seed )) ||
	  (!forwardBlocks.empty() && isNext( seed, forwardBlocks.back() )) );

  // splice together the two extensions and the seed (a bit messy):

  if( !blocks.empty() && isNext( blocks.back(), seed ) ){
    blocks.back().size += seed.size;
  }
  else blocks.push_back(seed);

  if( !forwardBlocks.empty() && isNext( seed, forwardBlocks.back() ) ){
    blocks.back().size += forwardBlocks.back().size;
    forwardBlocks.pop_back();
  }

  blocks.insert( blocks.end(), forwardBlocks.rbegin(), forwardBlocks.rend() );

  if( outputType > 3 ){  // set the un-ambiguity of the core to a max value:
    columnAmbiguityCodes.insert( columnAmbiguityCodes.end(), seed.size, 126 );
  }

  columnAmbiguityCodes.insert( columnAmbiguityCodes.end(),
                               forwardAmbiguities.rbegin(),
                               forwardAmbiguities.rend() );
}

bool Alignment::isOptimal( const uchar* seq1, const uchar* seq2,
			   const int scoreMatrix[MAT][MAT], int maxDrop,
			   const GeneralizedAffineGapCosts& gap,
			   int frameshiftCost, indexT frameSize,
			   const int pssm2[][MAT],
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2 ){
  int maxScore = 0;
  int runningScore = 0;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      indexT gapBeg1 = (i-1)->end1();
      indexT gapEnd1 = i->beg1();
      indexT gapSize1 = gapEnd1 - gapBeg1;

      indexT gapBeg2 = (i-1)->end2();
      indexT gapEnd2 = i->beg2();
      indexT gapSize2, frameshift2;
      sizeAndFrameshift( gapBeg2, gapEnd2, frameSize, gapSize2, frameshift2 );
      if( frameshift2 ) runningScore -= frameshiftCost;

      runningScore -= gap.cost( gapSize1, gapSize2 );
      if( runningScore <= 0 || runningScore < maxScore - maxDrop ){
	return false;
      }
    }

    const uchar* s1 = seq1 + i->beg1();
    const uchar* s2 = seq2 + i->beg2();
    const uchar* e1 = seq1 + i->end1();
    const int (*p2)[MAT] = pssm2 ? pssm2 + i->beg2() : 0;
    const uchar* q1 = qual1 ? qual1 + i->beg1() : 0;
    const uchar* q2 = qual2 ? qual2 + i->beg2() : 0;

    while( s1 < e1 ){
      /**/ if( sm2qual ) runningScore += sm2qual( *s1++, *s2++, *q1++, *q2++ );
      else if( pssm2 )   runningScore += ( *p2++ )[ *s1++ ];
      else               runningScore += scoreMatrix[ *s1++ ][ *s2++ ];

      if( runningScore > maxScore ) maxScore = runningScore;
      else if( runningScore <= 0 ||                  // non-optimal prefix
	       (s1 == e1 && i+1 == blocks.end()) ||  // non-optimal suffix
	       runningScore < maxScore - maxDrop ){  // excessive score drop
	return false;
      }
    }
  }

  return true;
}

void Alignment::extend( std::vector< SegmentPair >& chunks,
			std::vector< uchar >& ambiguityCodes,
			GappedXdropAligner& aligner, Centroid& centroid,
			const uchar* seq1, const uchar* seq2,
			indexT start1, indexT start2, bool isForward,
			const int sm[MAT][MAT], int smMax, int maxDrop,
			const GeneralizedAffineGapCosts& gap,
			int frameshiftCost, indexT frameSize,
			const int pssm2[][MAT],
                        const TwoQualityScoreMatrix& sm2qual,
                        const uchar* qual1, const uchar* qual2,
			double gamma, int outputType ){
  if( frameSize ){
    assert( outputType < 4 );
    assert( !pssm2 );
    assert( !sm2qual );

    indexT f = aaToDna( start2, frameSize ) + 1;
    indexT r = aaToDna( start2, frameSize ) - 1;

    const uchar* frame0 = seq2 + start2;
    const uchar* frame1 = seq2 + dnaToAa( isForward ? f : r, frameSize );
    const uchar* frame2 = seq2 + dnaToAa( isForward ? r : f, frameSize );

    score += aligner.align3( seq1 + start1, frame0, frame1, frame2, isForward,
                             sm, gap.exist, gap.extend, gap.extendPair,
                             frameshiftCost, maxDrop, smMax );

    std::size_t end1, end2, size;
    // This should be OK even if end2 < size * 3:
    while( aligner.getNextChunk3( end1, end2, size,
                                  gap.exist, gap.extend, gap.extendPair,
                                  frameshiftCost ) )
      chunks.push_back( SegmentPair( end1 - size, end2 - size * 3, size ) );

    return;
  }

  score +=
      sm2qual ? aligner.align2qual( seq1 + start1, qual1 + start1,
                                    seq2 + start2, qual2 + start2,
                                    isForward, sm2qual,
                                    gap.exist, gap.extend, gap.extendPair,
                                    maxDrop, smMax )
      : pssm2 ? aligner.alignPssm( seq1 + start1, pssm2 + start2, isForward,
                                   gap.exist, gap.extend, gap.extendPair,
                                   maxDrop, smMax )
      :         aligner.align( seq1 + start1, seq2 + start2, isForward,
                               sm, gap.exist, gap.extend, gap.extendPair,
                               maxDrop, smMax );

  if( outputType < 5 ){  // ordinary alignment, not gamma-centroid
    std::size_t end1, end2, size;
    while( aligner.getNextChunk( end1, end2, size,
                                 gap.exist, gap.extend, gap.extendPair ) )
      chunks.push_back( SegmentPair( end1 - size, end2 - size, size ) );
  }

  if( outputType > 3 ){  // calculate match probabilities
    assert( !sm2qual );
    centroid.reset();
    centroid.forward( seq1, seq2, start1, start2, isForward, gap );
    centroid.backward( seq1, seq2, start1, start2, isForward, gap );

    if( outputType > 4 ){  // do gamma-centroid alignment
      centroid.dp( gamma );
      centroid.traceback( chunks, gamma );
    }

    centroid.getColumnAmbiguities( ambiguityCodes, chunks, isForward );
  }
}
