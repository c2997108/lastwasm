// Copyright 2008 Martin C. Frith

// This struct holds a gapped, pair-wise alignment.

#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH
#include "SegmentPair.hh"
#include <string>
#include <vector>
#include <iosfwd>

namespace cbrc{

class XdropAligner;
class GeneralizedAffineGapCosts;
class MultiSequence;
class Alphabet;

struct Alignment{
  typedef unsigned indexT;
  typedef unsigned char uchar;

  enum { MAT = 64 };

  // make a single-block alignment:
  void fromSegmentPair( const SegmentPair& sp );

  // Make an Alignment by doing gapped X-drop alignment in both
  // directions starting from a seed SegmentPair.  The resulting
  // Alignment might not be "optimal" (see below).
  void makeXdrop( XdropAligner& aligner,
		  const uchar* seq1, const uchar* seq2,
		  const int scoreMatrix[MAT][MAT], int maxDrop,
		  const GeneralizedAffineGapCosts& gap );

  // Check that the Alignment has no prefix with score <= 0, no suffix
  // with score <= 0, and no sub-segment with score < -maxDrop.
  // Alignments that pass this test may be non-optimal in other ways.
  bool isOptimal( const uchar* seq1, const uchar* seq2,
                  const int scoreMatrix[MAT][MAT], int maxDrop,
                  const GeneralizedAffineGapCosts& gap );

  void write( const MultiSequence& seq1, const MultiSequence& seq2,
	      char strand, const Alphabet& alph, int format,
	      std::ostream& os ) const;

  // data:
  std::vector<SegmentPair> blocks;  // the gapless blocks of the alignment
  int score;
  SegmentPair seed;  // the alignment remembers its seed

  indexT beg1() const{ return blocks.front().beg1(); }
  indexT beg2() const{ return blocks.front().beg2(); }
  indexT end1() const{ return blocks.back().end1(); }
  indexT end2() const{ return blocks.back().end2(); }
  indexT range1() const{ return end1() - beg1(); }
  indexT range2() const{ return end2() - beg2(); }

  void writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
		 char strand, std::ostream& os ) const;

  void writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
		 char strand, const Alphabet& alph, std::ostream& os ) const;

  std::string topString( const std::vector<uchar>& seq,
			 const Alphabet& alph ) const;

  std::string botString( const std::vector<uchar>& seq,
			 const Alphabet& alph ) const;
};

}  // end namespace cbrc
#endif  // ALIGNMENT_HH
