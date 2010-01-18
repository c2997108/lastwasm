// Copyright 2008, 2009, 2010 Martin C. Frith

// This class holds multiple sequences and their names.  The sequences
// are concatenated, with delimiters between them.

// The final sequence may be "unfinished".  This happens if the
// sequence data hits a memory limit before we finish reading it.

#ifndef MULTISEQUENCE_HH
#define MULTISEQUENCE_HH

#include <vector>
#include <string>
#include <iosfwd>

namespace cbrc{

class MultiSequence{
 public:
  typedef unsigned indexT;
  typedef unsigned char uchar;

  // initialize with leftmost delimiter pad, ready for appending sequences
  void initForAppending( indexT padSizeIn );

  // re-initialize, but keep the last sequence if it is unfinished
  void reinitForAppending();

  // read seqCount finished sequences, and their names, from binary files
  void fromFiles( const std::string& baseName, indexT seqCount );

  // write all the finished sequences and their names to binary files
  void toFiles( const std::string& baseName ) const;

  // Append a sequence with delimiters.  Don't let the total size of
  // the concatenated sequences plus pads exceed maxBytes: thus it may
  // not finish reading the sequence.
  std::istream& appendFromFasta( std::istream& stream, std::size_t maxBytes );

  // As above, but read quality scores too.
  std::istream& appendFromFastq( std::istream& stream, std::size_t maxBytes );

  // As above, but read quality scores too.
  std::istream& appendFromPrb( std::istream& stream, std::size_t maxBytes,
			       unsigned alphSize, const uchar decode[] );

  // finish the last sequence: add final pad and end coordinate
  void finish();

  // unfinish the last sequence: remove final pad and end coordinate
  void unfinish();

  // did we finish reading the last sequence?
  bool isFinished() const{ return ends.size() == nameEnds.size(); }

  // how many finished sequences are there?
  indexT finishedSequences() const{ return ends.size() - 1; }

  // total length of finished sequences plus delimiters
  indexT finishedSize() const{ return ends.back(); }

  // total length of finished and unfinished sequences plus delimiters
  indexT unfinishedSize() const{ return seq.size(); }

  // which sequence is the coordinate in?
  indexT whichSequence( indexT coordinate ) const;

  indexT seqBeg( indexT seqNum ) const{ return ends[seqNum]; }
  indexT seqEnd( indexT seqNum ) const{ return ends[seqNum+1] - padSize; }
  indexT seqLen( indexT seqNum ) const{ return seqEnd(seqNum) - ends[seqNum]; }
  std::string seqName( indexT seqNum ) const;

  // get a pointer to the start of the sequence data
  const uchar* seqReader() const{ return &seq[0]; }
  /***/ uchar* seqWriter()      { return &seq[0]; }

  // swap the sequence data with some other sequence data
  void swapSeq( std::vector<uchar>& otherSeq ){ seq.swap(otherSeq); }

  // get a pointer to the start of the quality data
  const uchar* qualityReader() const{ return &qualityScores[0]; }
  /***/ uchar* qualityWriter()      { return &qualityScores[0]; }

  // How many quality scores are there per letter?  There might be
  // none at all, one per letter, or several (e.g. 4) per letter.
  unsigned qualsPerLetter() const
  { return qualityScores.size() / seq.size(); }

 private:
  indexT padSize;  // number of delimiter chars between sequences
  std::vector<uchar> seq;  // concatenated sequences
  std::vector<indexT> ends;  // coordinates of ends of delimiter pads
  std::vector<char> names;  // concatenated sequence names (to save memory)
  std::vector<indexT> nameEnds;  // endpoints of the names

  // The quality scores may be ASCII-coded: to get the real scores,
  // subtract e.g. 33 or 64.  The real scores might be related to
  // error probabilities in one of these ways:
  // Qphred = -10*log10(p)
  // Qsolexa = -10*log10(p/(1-p))
  std::vector<uchar> qualityScores;

  // read a FASTA header: read the whole line but store just the first word
  std::istream& readFastaName( std::istream& stream );

  // can we finish the last sequence and stay within the memory limit?
  bool isFinishable( std::size_t maxBytes ) const;
};

}  // end namespace cbrc
#endif  // MULTISEQUENCE_HH
