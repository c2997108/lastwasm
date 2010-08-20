// Copyright 2008, 2009, 2010 Martin C. Frith

// This struct holds a pair of equal-length segments, in a pair of
// sequences.  In other words, it holds a gapless alignment.

#ifndef SEGMENTPAIR_HH
#define SEGMENTPAIR_HH

namespace cbrc{

struct SegmentPair{
  typedef unsigned indexT;
  typedef unsigned char uchar;

  enum { MAT = 64 };

  SegmentPair(){}

  SegmentPair( indexT s1, indexT s2, indexT sz, int sc = 0 )
    : start1(s1), start2(s2), size(sz), score(sc){}

  // Make a SegmentPair by doing gapless X-drop alignment in both
  // directions starting from a seed point.  The resulting SegmentPair
  // might not be "optimal" (see below).
  void makeXdrop( indexT seed1, indexT seed2,
		  const uchar* seq1, const uchar* seq2,
                  const int scoreMatrix[MAT][MAT], int maxDrop,
		  const int pssm2[][MAT] );

  // Check that the SegmentPair has no prefix with score <= 0, no
  // suffix with score <= 0, and no sub-segment with score < -maxDrop.
  bool isOptimal( const uchar* seq1, const uchar* seq2,
		  const int scoreMatrix[MAT][MAT], int maxDrop,
		  const int pssm2[][MAT] ) const;

  // Shrink the SegmentPair to the longest run of identical letters
  // within it.  Allow (upper/lower)case to differ, using "canonical".
  void maxIdenticalRun( const uchar* seq1, const uchar* seq2,
			const uchar* canonical,
			const int scoreMatrix[MAT][MAT],
			const int pssm2[][MAT] );

  indexT beg1() const{ return start1; }         // start in sequence 1
  indexT beg2() const{ return start2; }         // start in sequence 2
  indexT end1() const{ return start1 + size; }  // end in sequence 1
  indexT end2() const{ return start2 + size; }  // end in sequence 2
  indexT diagonal() const{ return start1 - start2; }  // may wrap around!

  bool operator==( const SegmentPair& s ) const{
    return start1 == s.start1
        && start2 == s.start2
        && size == s.size
        && score == s.score; }

  indexT start1;
  indexT start2;
  indexT size;
  int score;

  // Since makeXdrop() is one of the most speed-critical routines, we
  // have two (almost-identical) versions: one using an ordinary score
  // matrix, the other using a Position-Specific-Score-Matrix.

  void makeXdrop( indexT seed1, indexT seed2,
		  const uchar* seq1, const uchar* seq2,
                  const int scoreMatrix[MAT][MAT], int maxDrop );

  void makeXdropPssm( indexT seed1, indexT seed2,
		      const uchar* seq1,
		      const int pssm2[][MAT], int maxDrop );

  // Used by maxIdenticalRun().
  void calculateScore( const uchar* seq1, const uchar* seq2,
		       const int scoreMatrix[MAT][MAT],
		       const int pssm2[][MAT] );
};

}  // end namespace cbrc
#endif  // SEGMENTPAIR_HH
