// Copyright 2008, 2009 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include "stringify.hh"
#include <iomanip>
#include <algorithm>
#include <cassert>

// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

// write x - y as a signed integer
static void writeSignedDifference( unsigned x, unsigned y, std::ostream& os ){
  if( x >= y )  os << x - y;
  else          os << '-' << y - x;
}

void Alignment::write( const MultiSequence& seq1, const MultiSequence& seq2,
		       char strand, bool isTranslated, const Alphabet& alph,
		       int format, std::ostream& os ) const{
  assert( !blocks.empty() );
  if( format == 0 ) writeTab( seq1, seq2, strand, isTranslated, os );
  else              writeMaf( seq1, seq2, strand, isTranslated, alph, os );
}

void Alignment::writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated,
			  std::ostream& os ) const{
  indexT alnBeg1 = beg1();
  indexT alnEnd1 = end1();
  indexT w1 = seq1.whichSequence(alnBeg1);
  indexT seqStart1 = seq1.seqBeg(w1);

  indexT size2 = seq2.ends.back();
  indexT frameSize2 = isTranslated ? (size2 / 3) : 0;
  indexT alnBeg2 = aaToDna( beg2(), frameSize2 );
  indexT alnEnd2 = aaToDna( end2(), frameSize2 );
  indexT w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  indexT seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  if( centroidScore < 0 ) os << score << '\t';
  else                    os << centroidScore << '\t';

  os << seq1.seqName(w1) << '\t'
     << alnBeg1 - seqStart1 << '\t'
     << alnEnd1 - alnBeg1 << '\t'
     << '+' << '\t'
     << seq1.seqLen(w1) << '\t';

  os << seq2.seqName(w2) << '\t'
     << alnBeg2 - seqStart2 << '\t'
     << alnEnd2 - alnBeg2 << '\t'
     << strand << '\t'
     << seq2.seqLen(w2) << '\t';

  for( unsigned i = 0; i < blocks.size(); ++i ){
    if( i > 0 ){  // between each pair of aligned blocks:
      indexT gapBeg1 = blocks[i-1].end1();
      indexT gapEnd1 = blocks[i].beg1();
      writeSignedDifference( gapEnd1, gapBeg1, os );  // allow -1 frameshift
      os << ':';

      indexT gapBeg2 = aaToDna( blocks[i-1].end2(), frameSize2 );
      indexT gapEnd2 = aaToDna( blocks[i].beg2(), frameSize2 );
      writeSignedDifference( gapEnd2, gapBeg2, os );  // allow -1 frameshift
      os << ',';
    }
    os << blocks[i].size << ( i+1 < blocks.size() ? ',' : '\n' );
  }
}

void Alignment::writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated, const Alphabet& alph,
			  std::ostream& os ) const{
  indexT alnBeg1 = beg1();
  indexT alnEnd1 = end1();
  indexT w1 = seq1.whichSequence(alnBeg1);
  indexT seqStart1 = seq1.seqBeg(w1);

  indexT size2 = seq2.ends.back();
  indexT frameSize2 = isTranslated ? (size2 / 3) : 0;
  indexT alnBeg2 = aaToDna( beg2(), frameSize2 );
  indexT alnEnd2 = aaToDna( end2(), frameSize2 );
  indexT w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  indexT seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);
  const std::string b1 = stringify( alnBeg1 - seqStart1 );
  const std::string b2 = stringify( alnBeg2 - seqStart2 );
  const std::string r1 = stringify( alnEnd1 - alnBeg1 );
  const std::string r2 = stringify( alnEnd2 - alnBeg2 );
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
     << std::setw( sw ) << s1 << ' '
     << topString( seq1.seqBase(), alph, frameSize2 ) << '\n';

  os << "s "
     << std::setw( nw ) << std::left << n2 << std::right << ' '
     << std::setw( bw ) << b2 << ' '
     << std::setw( rw ) << r2 << ' ' << strand << ' '
     << std::setw( sw ) << s2 << ' '
     << botString( seq2.seqBase(), alph, frameSize2 ) << '\n';

  if( seq2.qualsPerLetter() > 0 ){
    os << "q "
       << std::setw( nw ) << std::left << n2 << std::right << ' '
       << std::setw( bw + rw + sw + 5 ) << ""
       << qualityString( seq2.qualityBase(), seq2.qualsPerLetter() ) << '\n';
  }

  if( matchProbabilities.size() > 0 ){
    os << 'p';
    CI(double) p = matchProbabilities.begin();
    for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
      if( i > blocks.begin() ){  // between each pair of aligned blocks:
	// assume we're not doing translated alignment
	for( indexT j = (i-1)->end1(); j < i->beg1(); ++j ) os << " -";
	for( indexT j = (i-1)->end2(); j < i->beg2(); ++j ) os << " -";
      }
      for( indexT j = 0; j < i->size; ++j ) os << ' ' << *p++;
    }
    os << '\n';
  }

  os << '\n';  // blank line afterwards
}

std::string Alignment::topString( const uchar* seq, const Alphabet& alph,
				  indexT frameSize ) const{
  std::string s;

  for( unsigned i = 0; i < blocks.size(); ++i ){
    if( i > 0 ){  // between each pair of aligned blocks:
      indexT gapSize, frameshift;

      // append unaligned chunk of top sequence:
      indexT gapBeg1 = blocks[i-1].end1();
      indexT gapEnd1 = blocks[i].beg1();
      s.append( alph.rtString( seq + gapBeg1, seq + gapEnd1 ) );

      // append gaps for unaligned chunk of bottom sequence:
      indexT gapBeg2 = blocks[i-1].end2();
      indexT gapEnd2 = blocks[i].beg2();
      sizeAndFrameshift( gapBeg2, gapEnd2, frameSize, gapSize, frameshift );
      if( frameshift ) s.push_back( '-' );
      s.append( gapSize, '-' );
    }

    // append aligned chunk of top sequence:
    s.append( alph.rtString( seq + blocks[i].beg1(),
			     seq + blocks[i].end1() ) );
  }

  return s;
}

std::string Alignment::botString( const uchar* seq, const Alphabet& alph,
				  indexT frameSize ) const{
  std::string s;

  for( unsigned i = 0; i < blocks.size(); ++i ){
    if( i > 0 ){  // between each pair of aligned blocks:
      indexT gapSize, frameshift;

      // append gaps for unaligned chunk of top sequence:
      indexT gapBeg1 = blocks[i-1].end1();
      indexT gapEnd1 = blocks[i].beg1();
      s.append( gapEnd1 - gapBeg1, '-' );

      //append unaligned chunk of bottom sequence:
      indexT gapBeg2 = blocks[i-1].end2();
      indexT gapEnd2 = blocks[i].beg2();
      sizeAndFrameshift( gapBeg2, gapEnd2, frameSize, gapSize, frameshift );
      if( frameshift == 1 ) s.push_back( '\\' );
      if( frameshift == 2 ) s.push_back( '/' );
      s.append( alph.rtString( seq + gapEnd2 - gapSize, seq + gapEnd2 ) );
    }

    // append aligned chunk of bottom sequence:
    s.append( alph.rtString( seq + blocks[i].beg2(),
			     seq + blocks[i].end2() ) );
  }

  return s;
}

std::string Alignment::qualityString( const uchar* qualities,
				      unsigned qualsPerBase ) const{
  std::string s;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      // assume we're not doing translated alignment

      // append gaps for unaligned chunk of top sequence:
      s.append( i->beg1() - (i-1)->end1(), '-' );

      //append qualities for unaligned chunk of bottom sequence:
      s.append( qualityBlock( qualities, (i-1)->end2(), i->beg2(),
			      qualsPerBase ) );
    }

    // append qualities for aligned chunk of bottom sequence:
    s.append( qualityBlock( qualities, i->beg2(), i->end2(), qualsPerBase ) );
  }

  return s;
}

std::string Alignment::qualityBlock( const uchar* qualities,
				     indexT beg, indexT end,
				     unsigned qualsPerBase ){
  std::string s;
  for( indexT i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    s += *std::max_element( q, q + qualsPerBase );
  }
  return s;
}
