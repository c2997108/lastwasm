// Copyright 2008, 2009 Martin C. Frith

#include "Alignment.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include "stringify.hh"
#include <iomanip>
#include <algorithm>
#include <cassert>

// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

void Alignment::write( const MultiSequence& seq1, const MultiSequence& seq2,
		       char strand, const Alphabet& alph, int format,
		       std::ostream& os ) const{
  /**/ if( format == 0 ) writeTab( seq1, seq2, strand, os );
  else if( format == 1 ) writeMaf( seq1, seq2, strand, alph, os );
}

void Alignment::writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, std::ostream& os ) const{
  indexT size2 = seq2.ends.back();
  indexT w1 = seq1.whichSequence( beg1() );
  indexT w2 = seq2.whichSequence( strand == '+' ? beg2() : size2 - beg2() );
  indexT seqStart1 = seq1.seqBeg(w1);
  indexT seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  if( centroidScore < 0 ) os << score << '\t';
  else                    os << centroidScore << '\t';

  os << seq1.seqName(w1) << '\t'
     << beg1() - seqStart1 << '\t'
     << range1() << '\t'
     << '+' << '\t'
     << seq1.seqLen(w1) << '\t'
     << seq2.seqName(w2) << '\t'
     << beg2() - seqStart2 << '\t'
     << range2() << '\t'
     << strand << '\t'
     << seq2.seqLen(w2) << '\t';

  for( unsigned i = 0; i < blocks.size(); ++i ){
    if( i > 0 ){
      os << blocks[i].beg1() - blocks[i-1].end1() << ':'
	 << blocks[i].beg2() - blocks[i-1].end2() << ',';
    }
    os << blocks[i].size << ( i+1 < blocks.size() ? ',' : '\n' );
  }
}

void Alignment::writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, const Alphabet& alph,
			  std::ostream& os ) const{
  indexT size2 = seq2.ends.back();
  indexT w1 = seq1.whichSequence( beg1() );
  indexT w2 = seq2.whichSequence( strand == '+' ? beg2() : size2 - beg2() );
  indexT seqStart1 = seq1.seqBeg(w1);
  indexT seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);
  const std::string b1 = stringify( beg1() - seqStart1 );
  const std::string b2 = stringify( beg2() - seqStart2 );
  const std::string r1 = stringify( range1() );
  const std::string r2 = stringify( range2() );
  const std::string s1 = stringify( seq1.seqLen(w1) );
  const std::string s2 = stringify( seq2.seqLen(w2) );

  const std::size_t nw = std::max( n1.size(), n2.size() );
  const std::size_t bw = std::max( b1.size(), b2.size() );
  const std::size_t rw = std::max( r1.size(), r2.size() );
  const std::size_t sw = std::max( s1.size(), s2.size() );

  if( centroidScore < 0 ) os << "a score=" << score << '\n';
  else                    os << "a score=" << centroidScore << '\n';

  os << "s "
     << std::setw( nw ) << std::left << n1 << std::right << ' '
     << std::setw( bw ) << b1 << ' '
     << std::setw( rw ) << r1 << ' ' << '+' << ' '
     << std::setw( sw ) << s1 << ' ' << topString( seq1.seq, alph ) << '\n'
     << "s "
     << std::setw( nw ) << std::left << n2 << std::right << ' '
     << std::setw( bw ) << b2 << ' '
     << std::setw( rw ) << r2 << ' ' << strand << ' '
     << std::setw( sw ) << s2 << ' ' << botString( seq2.seq, alph ) << '\n';

  if( matchProbabilities.size() > 0 ){
    os << 'p';
    CI(double) p = matchProbabilities.begin();
    for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
      if( i > blocks.begin() ){
	for( indexT j = (i-1)->end1(); j < i->beg1(); ++j ) os << " -";
	for( indexT j = (i-1)->end2(); j < i->beg2(); ++j ) os << " -";
      }
      for( indexT j = 0; j < i->size; ++j ) os << ' ' << *p++;
    }
    os << '\n';
  }

  os << '\n';  // blank line afterwards
}

std::string Alignment::topString( const std::vector<uchar>& seq,
				  const Alphabet& alph ) const{
  std::string s;
  for( unsigned i = 0; i < blocks.size(); ++i ){
    if( i > 0 ){
      s.append( alph.rtString( seq.begin() + blocks[i-1].end1(),
			       seq.begin() + blocks[i].beg1() ) );
      s.append( blocks[i].beg2() - blocks[i-1].end2(), '-' );
    }
    s.append( alph.rtString( seq.begin() + blocks[i].beg1(),
			     seq.begin() + blocks[i].end1() ) );
  }
  return s;
}

std::string Alignment::botString( const std::vector<uchar>& seq,
				  const Alphabet& alph ) const{
  std::string s;
  for( unsigned i = 0; i < blocks.size(); ++i ){
    if( i > 0 ){
      s.append( blocks[i].beg1() - blocks[i-1].end1(), '-' );
      s.append( alph.rtString( seq.begin() + blocks[i-1].end2(),
			       seq.begin() + blocks[i].beg2() ) );
    }
    s.append( alph.rtString( seq.begin() + blocks[i].beg2(),
			     seq.begin() + blocks[i].end2() ) );
  }
  return s;
}
