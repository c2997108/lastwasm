// Copyright 2008, 2009, 2010, 2011, 2012 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstdio>  // sprintf
#include <iterator>  // ostream_iterator

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

  indexT size2 = seq2.finishedSize();
  indexT frameSize2 = isTranslated ? (size2 / 3) : 0;
  indexT alnBeg2 = aaToDna( beg2(), frameSize2 );
  indexT alnEnd2 = aaToDna( end2(), frameSize2 );
  indexT w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  indexT seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  os << score << '\t';

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

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;
      os << ',';
      indexT gapBeg1 = j->end1();
      indexT gapEnd1 = i->beg1();
      writeSignedDifference( gapEnd1, gapBeg1, os );  // allow -1 frameshift
      os << ':';
      indexT gapBeg2 = aaToDna( j->end2(), frameSize2 );
      indexT gapEnd2 = aaToDna( i->beg2(), frameSize2 );
      writeSignedDifference( gapEnd2, gapBeg2, os );  // allow -1 frameshift
      os << ',';
    }
    os << i->size;
  }

  os << '\n';
}

static int mySprintf( char* buffer, unsigned x ){
  return std::sprintf( buffer, "%u", x );
}

void Alignment::writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated, const Alphabet& alph,
			  std::ostream& os ) const{
  indexT alnBeg1 = beg1();
  indexT alnEnd1 = end1();
  indexT w1 = seq1.whichSequence(alnBeg1);
  indexT seqStart1 = seq1.seqBeg(w1);

  indexT size2 = seq2.finishedSize();
  indexT frameSize2 = isTranslated ? (size2 / 3) : 0;
  indexT alnBeg2 = aaToDna( beg2(), frameSize2 );
  indexT alnEnd2 = aaToDna( end2(), frameSize2 );
  indexT w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  indexT seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);

  char b1[32], b2[32], r1[32], r2[32], s1[32], s2[32];
  // sprintf is faster than ostringstream
  const int b1size = mySprintf( b1, alnBeg1 - seqStart1 );
  const int b2size = mySprintf( b2, alnBeg2 - seqStart2 );
  const int r1size = mySprintf( r1, alnEnd1 - alnBeg1 );
  const int r2size = mySprintf( r2, alnEnd2 - alnBeg2 );
  const int s1size = mySprintf( s1, seq1.seqLen(w1) );
  const int s2size = mySprintf( s2, seq2.seqLen(w1) );

  const int nw = std::max( n1.size(), n2.size() );
  const int bw = std::max( b1size, b2size );
  const int rw = std::max( r1size, r2size );
  const int sw = std::max( s1size, s2size );

  std::size_t numCols = numColumns( frameSize2 );
  std::vector<char> bufferVector( numCols + 1 );  // include a final NUL
  char* buffer = &bufferVector[0];

  os << "a";
  os << " score=" << score;
  os << '\n';

  writeTopSeq( seq1.seqReader(), alph, frameSize2, buffer );
  os << "s "
     << std::setw( nw ) << std::left << n1 << std::right << ' '
     << std::setw( bw ) << b1 << ' '
     << std::setw( rw ) << r1 << ' ' << '+' << ' '
     << std::setw( sw ) << s1 << ' '
     << buffer << '\n';

  if( seq1.qualsPerLetter() > 0 ){
    writeTopQual( seq1.qualityReader(), seq1.qualsPerLetter(), buffer );
    os << "q "
       << std::setw( nw ) << std::left << n1 << std::right << ' '
       << std::setw( bw + rw + sw + 5 ) << ""
       << buffer << '\n';
  }

  writeBotSeq( seq2.seqReader(), alph, frameSize2, buffer );
  os << "s "
     << std::setw( nw ) << std::left << n2 << std::right << ' '
     << std::setw( bw ) << b2 << ' '
     << std::setw( rw ) << r2 << ' ' << strand << ' '
     << std::setw( sw ) << s2 << ' '
     << buffer << '\n';

  if( seq2.qualsPerLetter() > 0 ){
    writeBotQual( seq2.qualityReader(), seq2.qualsPerLetter(), buffer );
    os << "q "
       << std::setw( nw ) << std::left << n2 << std::right << ' '
       << std::setw( bw + rw + sw + 5 ) << ""
       << buffer << '\n';
  }

  if( columnAmbiguityCodes.size() > 0 ){
    os << "p "
       << std::setw( nw + bw + rw + sw + 6 ) << "";
    std::copy( columnAmbiguityCodes.begin(), columnAmbiguityCodes.end(),
               std::ostream_iterator<uchar>(os) );
    os << '\n';
  }

  os << '\n';  // blank line afterwards
}

std::size_t Alignment::numColumns( indexT frameSize ) const{
  std::size_t num = 0;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // length of unaligned chunk of top sequence (gaps in bottom sequence):
      num += i->beg1() - j->end1();

      // length of unaligned chunk of bottom sequence (gaps in top sequence):
      indexT gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) ++num;
      num += gap2;
    }

    num += i->size;  // length of aligned chunk
  }

  return num;
}

static char* writeGaps( char* dest, Alignment::indexT num ){
  char* end = dest + num;
  while( dest < end ){
    *dest++ = '-';
  }
  return dest;
}

char* Alignment::writeTopSeq( const uchar* seq, const Alphabet& alph,
			      indexT frameSize, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // append unaligned chunk of top sequence:
      dest = alph.rtCopy( seq + j->end1(), seq + i->beg1(), dest );

      // append gaps for unaligned chunk of bottom sequence:
      indexT gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) *dest++ = '-';
      dest = writeGaps( dest, gap2 );
    }

    // append aligned chunk of top sequence:
    dest = alph.rtCopy( seq + i->beg1(), seq + i->end1(), dest );
  }

  return dest;
}

char* Alignment::writeBotSeq( const uchar* seq, const Alphabet& alph,
			      indexT frameSize, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, i->beg1() - j->end1() );

      //append unaligned chunk of bottom sequence:
      indexT gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 == 1 ) *dest++ = '\\';
      if( frameshift2 == 2 ) *dest++ = '/';
      dest = alph.rtCopy( seq + i->beg2() - gap2, seq + i->beg2(), dest );
    }

    // append aligned chunk of bottom sequence:
    dest = alph.rtCopy( seq + i->beg2(), seq + i->end2(), dest );
  }

  return dest;
}

static char* writeQuals( const uchar* qualities,
			 Alignment::indexT beg, Alignment::indexT end,
			 std::size_t qualsPerBase, char* dest ){
  for( Alignment::indexT i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    *dest++ = *std::max_element( q, q + qualsPerBase );
  }
  return dest;
}

char* Alignment::writeTopQual( const uchar* qualities,
			       std::size_t qualsPerBase, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // assume we're not doing translated alignment

      // append qualities for unaligned chunk of top sequence:
      dest = writeQuals( qualities, j->end1(), i->beg1(), qualsPerBase, dest );

      // append gaps for unaligned chunk of bottom sequence:
      dest = writeGaps( dest, i->beg2() - j->end2() );
    }

    // append qualities for aligned chunk of top sequence:
    dest = writeQuals( qualities, i->beg1(), i->end1(), qualsPerBase, dest );
  }

  return dest;
}

char* Alignment::writeBotQual( const uchar* qualities,
			       std::size_t qualsPerBase, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // assume we're not doing translated alignment

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, i->beg1() - j->end1() );

      // append qualities for unaligned chunk of bottom sequence:
      dest = writeQuals( qualities, j->end2(), i->beg2(), qualsPerBase, dest );
    }

    // append qualities for aligned chunk of bottom sequence:
    dest = writeQuals( qualities, i->beg2(), i->end2(), qualsPerBase, dest );
  }

  return dest;
}
