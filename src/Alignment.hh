// Copyright 2008, 2009, 2010, 2011, 2012 Martin C. Frith

// This struct holds a gapped, pair-wise alignment.

#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH
#include "SegmentPair.hh"
#include <cstddef>  // size_t
#include <string>
#include <vector>
#include <iosfwd>

namespace cbrc{

typedef unsigned char uchar;

class GappedXdropAligner;
class GeneralizedAffineGapCosts;
class MultiSequence;
class Alphabet;
class Centroid;
class TwoQualityScoreMatrix;

struct Alignment{
  typedef SegmentPair::indexT indexT;

  enum { MAT = 64 };

  // make a single-block alignment:
  void fromSegmentPair( const SegmentPair& sp );

  // Make an Alignment by doing gapped X-drop extension in both
  // directions starting from a seed SegmentPair.  The resulting
  // Alignment might not be "optimal" (see below).
  // If outputType > 3: calculates match probabilities.
  // If outputType > 4: does gamma-centroid alignment.
  void makeXdrop( GappedXdropAligner& aligner, Centroid& centroid,
		  const uchar* seq1, const uchar* seq2,
		  const int scoreMatrix[MAT][MAT], int smMax,
		  const GeneralizedAffineGapCosts& gap, int maxDrop,
		  int frameshiftCost, indexT frameSize,
		  const int pssm2[][MAT],
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2,
		  const Alphabet& alph, double gamma = 0, int outputType = 0 );

  // Check that the Alignment has no prefix with score <= 0, no suffix
  // with score <= 0, and no sub-segment with score < -maxDrop.
  // Alignments that pass this test may be non-optimal in other ways.
  bool isOptimal( const uchar* seq1, const uchar* seq2,
                  const int scoreMatrix[MAT][MAT], int maxDrop,
                  const GeneralizedAffineGapCosts& gap,
		  int frameshiftCost, indexT frameSize,
		  const int pssm2[][MAT],
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2 );

  void write( const MultiSequence& seq1, const MultiSequence& seq2,
	      char strand, bool isTranslated, const Alphabet& alph,
	      int format, std::ostream& os ) const;

  // data:
  std::vector<SegmentPair> blocks;  // the gapless blocks of the alignment
  int score;
  SegmentPair seed;  // the alignment remembers its seed
  std::vector<uchar> columnAmbiguityCodes;  // char or uchar?
  std::vector<double> expectedCounts;  // expected emission & transition counts

  indexT beg1() const{ return blocks.front().beg1(); }
  indexT beg2() const{ return blocks.front().beg2(); }
  indexT end1() const{ return blocks.back().end1(); }
  indexT end2() const{ return blocks.back().end2(); }

  void extend( std::vector< SegmentPair >& chunks,
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
	       const Alphabet& alph, double gamma, int outputType );

  void writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
		 char strand, bool isTranslated, std::ostream& os ) const;

  void writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
		 char strand, bool isTranslated, const Alphabet& alph,
		 std::ostream& os ) const;

  std::size_t numColumns( indexT frameSize ) const;

  char* writeTopSeq( const uchar* seq, const Alphabet& alph,
		     indexT frameSize, char* dest ) const;

  char* writeBotSeq( const uchar* seq, const Alphabet& alph,
		     indexT frameSize, char* dest ) const;

  char* writeTopQual( const uchar* qualities,
		      std::size_t qualsPerBase, char* dest ) const;

  char* writeBotQual( const uchar* qualities,
		      std::size_t qualsPerBase, char* dest ) const;
};

}  // end namespace cbrc
#endif  // ALIGNMENT_HH
