// Copyright 2008, 2009 Martin C. Frith

// This struct holds multiple sequences and their names.  The
// sequences are concatenated, with delimiters between them.

// The final sequence may be "unfinished".  This happens if the
// sequence data hits a memory limit before we finish reading it.

#ifndef MULTISEQUENCE_HH
#define MULTISEQUENCE_HH
#include <vector>
#include <string>
#include <iosfwd>

namespace cbrc{

struct MultiSequence{
  typedef unsigned indexT;
  typedef unsigned char uchar;

  // initialize with leftmost delimiter pad, ready for appending sequences
  void initForAppending( indexT padSizeIn );

  // re-initialize, but keep the last sequence if it is unfinished
  void reinitForAppending();

  // when reading from files, assume no. of sequences is known from elsewhere:
  void fromFiles( const std::string& baseName, indexT seqCount );
  void toFiles( const std::string& baseName ) const;

  // Append a sequence with delimiters.  Don't let the total size of
  // the concatenated sequences plus pads exceed maxBytes: thus it may
  // not finish reading the sequence.
  std::istream& appendFromFasta( std::istream& stream, std::size_t maxBytes );

  // read a FASTA header: read the whole line but store just the first word
  std::istream& readFastaName( std::istream& stream );

  // finish the last sequence: add final pad and end coordinate
  void finish();

  // unfinish the last sequence: remove final pad and end coordinate
  void unfinish();

  // did we finish reading the last sequence?
  bool isFinished() const{ return ends.size() == nameEnds.size(); }

  // how many finished sequences are there?
  indexT finishedSequences() const{ return ends.size() - 1; }

  // can we finish the last sequence and stay within the memory limit?
  bool isFinishable( std::size_t maxBytes ) const;

  // which sequence is the coordinate in?
  indexT whichSequence( indexT coordinate ) const;

  indexT seqBeg( indexT seqNum ) const{ return ends[seqNum]; }
  indexT seqEnd( indexT seqNum ) const{ return ends[seqNum+1] - padSize; }
  indexT seqLen( indexT seqNum ) const{ return seqEnd(seqNum) - ends[seqNum]; }
  std::string seqName( indexT seqNum ) const;

  // data:
  indexT padSize;  // number of delimiter chars between sequences
  std::vector<uchar> seq;  // concatenated sequences
  std::vector<indexT> ends;  // coordinates of ends of delimiter pads
  std::vector<char> names;  // concatenated sequence names (to save memory)
  std::vector<indexT> nameEnds;  // endpoints of the names
};

}  // end namespace cbrc
#endif  // MULTISEQUENCE_HH
