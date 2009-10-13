// Copyright 2008, 2009 Martin C. Frith

#include "Alignment.hh"
#include "Centroid.hh"
#include "GeneticCode.hh"
#include "GeneralizedAffineGapCosts.hh"
#include <cassert>

// make C++ tolerable:
#define IT(type) std::vector<type>::iterator
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

void Alignment::fromSegmentPair( const SegmentPair& sp ){
  blocks.assign( 1, sp );
  score = sp.score;
  centroidScore = -1;
}

void Alignment::makeXdrop( Xdrop3FrameAligner& aligner, Centroid& centroid,
			   const uchar* seq1, const uchar* seq2,
			   const int scoreMatrix[MAT][MAT], int smMax,
			   const GeneralizedAffineGapCosts& gap, int maxDrop,
			   int frameshiftCost, indexT frameSize,
			   const int pssm2[][MAT],
			   double gamma, int outputType ){
  assert( seed.size > 0 );  // relax this requirement?
  score = seed.score;
  centroidScore = (outputType < 5 ? -1 : gamma * seed.size);

  extend( blocks, matchProbabilities, aligner, centroid, seq1, seq2,
	  seed.beg1(), seed.beg2(), XdropAligner::REVERSE,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, gamma, outputType );

  // convert left-extension coordinates to sequence coordinates:
  indexT seedBeg1 = seed.beg1();
  indexT seedBeg2 = aaToDna( seed.beg2(), frameSize );
  for( IT(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    i->start1 = seedBeg1 - i->start1 - i->size;
    i->start2 = dnaToAa( seedBeg2 - i->start2, frameSize ) - i->size;
  }

  if( !blocks.empty() &&
      blocks.back().end1() == seed.beg1() &&
      blocks.back().end2() == seed.beg2() ){
    blocks.back().size += seed.size;
  }
  else blocks.push_back(seed);

  if( outputType > 3 ){  // set match probabilities in the seed to 1.0 (?)
    matchProbabilities.insert( matchProbabilities.end(), seed.size, 1.0 );
  }

  std::vector<SegmentPair> forwardBlocks;
  std::vector<double> forwardProbs;

  extend( forwardBlocks, forwardProbs, aligner, centroid, seq1, seq2,
	  seed.end1(), seed.end2(), XdropAligner::FORWARD,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, gamma, outputType );

  // convert right-extension coordinates to sequence coordinates:
  indexT seedEnd1 = seed.end1();
  indexT seedEnd2 = aaToDna( seed.end2(), frameSize );
  for( IT(SegmentPair) i = forwardBlocks.begin(); i < forwardBlocks.end();
       ++i ){
    i->start1 = seedEnd1 + i->start1;
    i->start2 = dnaToAa( seedEnd2 + i->start2, frameSize );
  }

  SegmentPair& b = blocks.back();
  if( !forwardBlocks.empty() &&
      b.end1() == forwardBlocks.back().beg1() &&
      b.end2() == forwardBlocks.back().beg2() ){
    b.size += forwardBlocks.back().size;
    forwardBlocks.pop_back();
  }

  blocks.insert( blocks.end(), forwardBlocks.rbegin(), forwardBlocks.rend() );
  matchProbabilities.insert( matchProbabilities.end(),
			     forwardProbs.rbegin(), forwardProbs.rend() );
}

bool Alignment::isOptimal( const uchar* seq1, const uchar* seq2,
			   const int scoreMatrix[MAT][MAT], int maxDrop,
			   const GeneralizedAffineGapCosts& gap,
			   int frameshiftCost, indexT frameSize,
			   const int pssm2[][MAT] ){
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

    while( s1 < e1 ){
      if( pssm2 ) runningScore += ( *p2++ )[ *s1++ ];
      else runningScore += scoreMatrix[ *s1++ ][ *s2++ ];

      if( runningScore > maxScore ) maxScore = runningScore;
      else if( runningScore <= 0 ||                  // non-optimal prefix
	       (s1 == e1 && i+1 == blocks.end()) ||   // non-optimal suffix
	       runningScore < maxScore - maxDrop ){  // excessive score drop
	return false;
      }
    }
  }

  return true;
}

void Alignment::extend( std::vector< SegmentPair >& chunks,
			std::vector< double >& probs,
			Xdrop3FrameAligner& aligner, Centroid& centroid,
			const uchar* seq1, const uchar* seq2,
			indexT start1, indexT start2,
			XdropAligner::direction dir,
			const int sm[MAT][MAT], int smMax, int maxDrop,
			const GeneralizedAffineGapCosts& gap,
			int frameshiftCost, indexT frameSize,
			const int pssm2[][MAT],
			double gamma, int outputType ){
  if( frameSize ){
    assert( outputType < 4 );
    assert( !pssm2 );
    score += aligner.fillThreeFrame( seq1, seq2, start1, start2, dir, sm,
				     maxDrop, gap, frameshiftCost, frameSize );
    aligner.traceThreeFrame( chunks, seq1, seq2, start1, start2, dir,
			     sm, gap, frameshiftCost, frameSize );
    return;
  }

  score += aligner.fill( seq1, seq2, start1, start2, dir,
			 sm, smMax, maxDrop, gap, pssm2 );

  if( outputType < 5 ){  // ordinary alignment, not gamma-centroid
    aligner.traceback( chunks, seq1, seq2, start1, start2, dir,
		       sm, pssm2, gap );
  }

  if( outputType > 3 ){  // calculate match probabilities
    assert( !pssm2 );
    centroid.reset();
    centroid.forward( seq1, seq2, start1, start2, dir, gap );
    centroid.backward( seq1, seq2, start1, start2, dir, gap );

    if( outputType > 4 ){  // do gamma-centroid alignment
      centroidScore += centroid.dp( gamma );
      centroid.traceback( chunks, gamma );
    }

    centroid.chunkProbabilities( probs, chunks );
  }
}
