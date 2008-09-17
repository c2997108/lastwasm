// Copyright 2008 Martin C. Frith

#include "Alignment.hh"
#include "XdropAligner.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "stringify.hh"
#include <iomanip>
#include <algorithm>
#include <cassert>

// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

namespace cbrc{

void Alignment::fromSegmentPair( const SegmentPair& sp ){
  blocks.assign( 1, sp );
  score = sp.score;
}

void Alignment::makeXdrop( XdropAligner& aligner,
			   const uchar* seq1, const uchar* seq2,
			   const int scoreMatrix[MAT][MAT], int maxDrop,
			   const GeneralizedAffineGapCosts& gap ){
  assert( seed.size > 0 );  // relax this requirement?
  score = seed.score;

  score += aligner.fill( seq1, seq2, seed.beg1(), seed.beg2(),
			 XdropAligner::REVERSE, scoreMatrix, maxDrop, gap );
  aligner.traceback( blocks, seq1, seq2, seed.beg1(), seed.beg2(),
		     XdropAligner::REVERSE, scoreMatrix, gap );

  for( unsigned i = 0; i < blocks.size(); ++i ){
    blocks[i].start1 = seed.beg1() - blocks[i].end1();
    blocks[i].start2 = seed.beg2() - blocks[i].end2();
  }

  if( !blocks.empty() &&
      blocks.back().end1() == seed.beg1() &&
      blocks.back().end2() == seed.beg2() ){
    blocks.back().size += seed.size;
  }
  else blocks.push_back(seed);

  std::vector<SegmentPair> forwardBlocks;

  score += aligner.fill( seq1, seq2, seed.end1(), seed.end2(),
			 XdropAligner::FORWARD, scoreMatrix, maxDrop, gap );
  aligner.traceback( forwardBlocks, seq1, seq2, seed.end1(), seed.end2(),
		     XdropAligner::FORWARD, scoreMatrix, gap );

  for( unsigned i = 0; i < forwardBlocks.size(); ++i ){
    forwardBlocks[i].start1 += seed.end1();
    forwardBlocks[i].start2 += seed.end2();
  }

  SegmentPair& b = blocks.back();
  if( !forwardBlocks.empty() &&
      b.end1() == forwardBlocks.back().beg1() &&
      b.end2() == forwardBlocks.back().beg2() ){
    b.size += forwardBlocks.back().size;
    forwardBlocks.pop_back();
  }

  blocks.insert( blocks.end(), forwardBlocks.rbegin(), forwardBlocks.rend() );
}

bool Alignment::isOptimal( const uchar* seq1, const uchar* seq2,
			   const int scoreMatrix[MAT][MAT], int maxDrop,
			   const GeneralizedAffineGapCosts& gap ){
  int maxScore = 0;
  int runningScore = 0;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){
      indexT gapSize1 = i->beg1() - (i-1)->end1();
      indexT gapSize2 = i->beg2() - (i-1)->end2();
      runningScore -= gap.cost( gapSize1, gapSize2 );
      if( runningScore <= 0 || runningScore < maxScore - maxDrop ){
	return false;
      }
    }

    const uchar* s1 = seq1 + i->beg1();
    const uchar* s2 = seq2 + i->beg2();
    const uchar* e1 = seq1 + i->end1();

    while( s1 < e1 ){
      runningScore += scoreMatrix[ *s1++ ][ *s2++ ];
      if( runningScore > maxScore ) maxScore = runningScore;
      else if( runningScore <= 0 ||                  // non-optimal prefix
	       s1 == e1 && i+1 == blocks.end() ||    // non-optimal suffix
	       runningScore < maxScore - maxDrop ){  // excessive score drop
	return false;
      }
    }
  }

  return true;
}

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

  os << score << '\t'
     << seq1.seqName(w1) << '\t'
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

  os << "a score=" << score << '\n'
     << "s "
     << std::setw( nw ) << std::left << n1 << std::right << ' '
     << std::setw( bw ) << b1 << ' '
     << std::setw( rw ) << r1 << ' ' << '+' << ' '
     << std::setw( sw ) << s1 << ' ' << topString( seq1.seq, alph ) << '\n'
     << "s "
     << std::setw( nw ) << std::left << n2 << std::right << ' '
     << std::setw( bw ) << b2 << ' '
     << std::setw( rw ) << r2 << ' ' << strand << ' '
     << std::setw( sw ) << s2 << ' ' << botString( seq2.seq, alph ) << '\n'
     << '\n';  // blank line afterwards
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

}  // end namespace cbrc
