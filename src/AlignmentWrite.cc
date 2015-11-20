// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "LastEvaluer.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include <algorithm>
#include <cassert>
#include <cstdio>  // sprintf
#include <iostream>

using namespace cbrc;

// This writes a "size_t" integer into a char buffer ending at "end".
// It writes backwards from the end, because that's easier & faster.
static char *writeSize(char *end, size_t x) {
  do {
    --end;
    *end = '0' + x % 10;
    x /= 10;
  } while (x);
  return end;
}

class IntText {  // a text representation of an integer
public:
  IntText() {}
  explicit IntText(size_t x) { set(x); }
  void set(size_t x) { char *e = b + sizeof b; s = e - writeSize(e, x); }
  const char *begin() const { return b + sizeof b - s; }
  size_t size() const { return s; }
private:
  char b[31];
  unsigned char s;
};

class FloatText {  // a text representation of a floating-point number
public:
  FloatText() {}
  FloatText(const char *format, double x) { set(format, x); }
  void set(const char *format, double x) { s = std::sprintf(b, format, x); }
  const char *begin() const { return b; }
  size_t size() const { return s; }
private:
  char b[31];
  unsigned char s;
};

class Writer {  // writes characters to an output buffer
public:
  explicit Writer(char *startOfOutput) : p(startOfOutput) {}
  Writer &operator<<(char c) {
    *p++ = c;
    return *this;
  }
  Writer &operator<<(const IntText &t) {
    std::memcpy(p, t.begin(), t.size());
    p += t.size();
    return *this;
  }
  Writer &operator<<(const FloatText &t) {
    std::memcpy(p, t.begin(), t.size());
    p += t.size();
    return *this;
  }
  Writer &operator<<(const std::string &t) {
    std::memcpy(p, t.c_str(), t.size());
    p += t.size();
    return *this;
  }
  void fill(size_t n, char c) {  // write n copies of character c
    std::memset(p, c, n);
    p += n;
  }
private:
  char *p;
};

// write x - y as a signed integer
static void writeSignedDifference( size_t x, size_t y, std::ostream& os ){
  if( x >= y )  os << x - y;
  else          os << '-' << y - x;
}

void Alignment::write( const MultiSequence& seq1, const MultiSequence& seq2,
		       char strand, const uchar* seqData2,
		       bool isTranslated, const Alphabet& alph,
		       const LastEvaluer& evaluer, int format,
		       std::vector<AlignmentText>& textAlns,
		       const AlignmentExtras& extras ) const{
  assert( !blocks.empty() );
  std::ostream& os = std::cout;
  if( format == 't' )
    writeTab( seq1, seq2, strand,
	      isTranslated, evaluer, os, extras );
  if( format == 'm' )
    writeMaf( seq1, seq2, strand, seqData2,
	      isTranslated, alph, evaluer, os, extras );
  if( format == 'b' )
    writeBlastTab( seq1, seq2, strand, seqData2,
		   isTranslated, alph, evaluer, textAlns );
}

static size_t alignedColumnCount(const std::vector<SegmentPair> &blocks) {
  size_t c = 0;
  for (size_t i = 0; i < blocks.size(); ++i)
    c += blocks[i].size;
  return c;
}

static size_t matchCount(const std::vector<SegmentPair> &blocks,
			 const uchar *seq1, const uchar *seq2,
			 const uchar *numbersToUppercase) {
  // no special treatment of ambiguous bases/residues: same as NCBI BLAST
  size_t matches = 0;
  for (size_t i = 0; i < blocks.size(); ++i) {
    const SegmentPair &b = blocks[i];
    const uchar *x = seq1 + b.beg1();
    const uchar *y = seq2 + b.beg2();
    for (size_t j = 0; j < b.size; ++j)
      if (numbersToUppercase[x[j]] == numbersToUppercase[y[j]])
	++matches;
  }
  return matches;
}

static char* writeTaggedItems( const LastEvaluer& evaluer, double queryLength,
			       int score, double fullScore,
			       char separator, char* out ){
  if( evaluer.isGood() ){
    double epa = evaluer.evaluePerArea( score );
    double area = evaluer.area( score, queryLength );
    *out++ = separator;
    out += std::sprintf( out, "EG2=%.2g", 1e18 * epa );
    *out++ = separator;
    out += std::sprintf( out, "E=%.2g", area * epa );
  }
  if( fullScore > 0 ){
    *out++ = separator;
    out += std::sprintf( out, "fullScore=%.3g", fullScore );
  }
  *out++ = '\n';
  return out;
}

void Alignment::writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated,
			  const LastEvaluer& evaluer, std::ostream& os,
			  const AlignmentExtras& extras ) const{
  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);
  size_t seqLen1 = seq1.seqLen(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);
  size_t seqLen2 = seq2.seqLen(w2);

  IntText sc(score);
  std::string n1 = seq1.seqName(w1);
  std::string n2 = seq2.seqName(w2);
  IntText b1(alnBeg1 - seqStart1);
  IntText b2(alnBeg2 - seqStart2);
  IntText r1(alnEnd1 - alnBeg1);
  IntText r2(alnEnd2 - alnBeg2);
  IntText s1(seqLen1);
  IntText s2(seqLen2);

  size_t s = sc.size() +
    n1.size() + b1.size() + r1.size() + 1 + s1.size() +
    n2.size() + b2.size() + r2.size() + 1 + s2.size() + 11;
  std::vector<char> v(s);
  Writer w(&v[0]);
  const char t = '\t';
  w << sc << t;
  w << n1 << t << b1 << t << r1 << t << '+'    << t << s1 << t;
  w << n2 << t << b2 << t << r2 << t << strand << t << s2 << t;
  os.write(&v[0], s);

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];
      if( x.size ) os << ',';
      size_t gapBeg1 = x.end1();
      size_t gapEnd1 = y.beg1();
      writeSignedDifference( gapEnd1, gapBeg1, os );  // allow -1 frameshift
      os << ':';
      size_t gapBeg2 = aaToDna( x.end2(), frameSize2 );
      size_t gapEnd2 = aaToDna( y.beg2(), frameSize2 );
      writeSignedDifference( gapEnd2, gapBeg2, os );  // allow -1 frameshift
      if( y.size ) os << ',';
    }
    if( y.size ) os << y.size;
  }

  double fullScore = extras.fullScore;
  char line[256];
  char* end = writeTaggedItems( evaluer, seqLen2, score, fullScore, t, line );
  os.write( line, end - line );
}

static void putLeft(Writer &w, const std::string &t, size_t width) {
  w << t;
  w.fill(width - t.size(), ' ');
  w << ' ';
}

static void putRight(Writer &w, const IntText &t, size_t width) {
  w.fill(width - t.size(), ' ');
  w << t;
  w << ' ';
}

static void writeMafData(char *out,
			 const std::string &n, size_t nw,
			 const IntText &b, size_t bw,
			 const IntText &r, size_t rw,
			 char strand,
			 const IntText &s, size_t sw) {
  Writer w(out);
  w << 's' << ' ';
  putLeft(w, n, nw);
  putRight(w, b, bw);
  putRight(w, r, rw);
  w << strand << ' ';
  putRight(w, s, sw);
}

void Alignment::writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, const uchar* seqData2,
			  bool isTranslated, const Alphabet& alph,
			  const LastEvaluer& evaluer, std::ostream& os,
			  const AlignmentExtras& extras ) const{
  double fullScore = extras.fullScore;
  const std::vector<uchar>& columnAmbiguityCodes = extras.columnAmbiguityCodes;
  const std::vector<double>& expectedCounts = extras.expectedCounts;

  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);
  size_t seqLen1 = seq1.seqLen(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);
  size_t seqLen2 = seq2.seqLen(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);
  IntText b1(alnBeg1 - seqStart1);
  IntText b2(alnBeg2 - seqStart2);
  IntText r1(alnEnd1 - alnBeg1);
  IntText r2(alnEnd2 - alnBeg2);
  IntText s1(seqLen1);
  IntText s2(seqLen2);

  const size_t nw = std::max( n1.size(), n2.size() );
  const size_t bw = std::max( b1.size(), b2.size() );
  const size_t rw = std::max( r1.size(), r2.size() );
  const size_t sw = std::max( s1.size(), s2.size() );

  size_t headLen = 2 + nw + 1 + bw + 1 + rw + 3 + sw + 1;
  size_t lineLen = headLen + numColumns( frameSize2 ) + 1;
  std::vector<char> lineVector( lineLen );
  char* line = &lineVector[0];
  char* tail = line + headLen;
  line[ lineLen - 1 ] = '\n';

  char aLine[256];
  char* dest = aLine;
  dest += std::sprintf( dest, "a score=%d", score );
  dest = writeTaggedItems( evaluer, seqLen2, score, fullScore, ' ', dest );
  os.write( aLine, dest - aLine );

  writeMafData( line, n1, nw, b1, bw, r1, rw, '+', s1, sw );
  writeTopSeq( seq1.seqReader(), alph, 0, frameSize2, tail );
  os.write( line, lineLen );

  size_t qualsPerBase1 = seq1.qualsPerLetter();
  if( qualsPerBase1 ){
    *line = 'q';
    std::fill( line + 2 + nw + 1, tail, ' ' );
    writeTopSeq( seq1.qualityReader(), alph, qualsPerBase1, frameSize2, tail );
    os.write( line, lineLen );
  }

  writeMafData( line, n2, nw, b2, bw, r2, rw, strand, s2, sw );
  writeBotSeq( seqData2, alph, 0, frameSize2, tail );
  os.write( line, lineLen );

  size_t qualsPerBase2 = seq2.qualsPerLetter();
  if( qualsPerBase2 && !isTranslated ){
    // for translated alignment: don't write untranslated quality data
    *line = 'q';
    std::fill( line + 2 + nw + 1, tail, ' ' );
    writeBotSeq( seq2.qualityReader(), alph, qualsPerBase2, frameSize2, tail );
    os.write( line, lineLen );
  }

  if( columnAmbiguityCodes.size() > 0 ){
    *line = 'p';
    std::fill( line + 2, tail, ' ' );
    copy( columnAmbiguityCodes.begin(), columnAmbiguityCodes.end(), tail );
    os.write( line, lineLen );
  }

  if( expectedCounts.size() > 0 ){
    os << 'c';
    for( unsigned i = 0; i < expectedCounts.size(); ++i )
      os << ' ' << expectedCounts[i];
    os << '\n';
  }

  os << '\n';  // blank line afterwards
}

void Alignment::writeBlastTab( const MultiSequence& seq1,
			       const MultiSequence& seq2,
			       char strand, const uchar* seqData2,
			       bool isTranslated, const Alphabet& alph,
			       const LastEvaluer& evaluer,
			       std::vector<AlignmentText>& textAlns ) const{
  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  if( strand == '-' ){
    alnBeg2 = size2 - alnBeg2;
    alnEnd2 = size2 - alnEnd2;
  }
  size_t w2 = seq2.whichSequence( alnBeg2 );
  size_t seqStart2 = seq2.seqBeg(w2);

  size_t alnSize = numColumns( frameSize2 );
  size_t matches = matchCount( blocks, seq1.seqReader(), seqData2,
			       alph.numbersToUppercase );
  size_t mismatches = alignedColumnCount(blocks) - matches;
  size_t gapOpens = blocks.size() - 1;
  double matchPercent = 100.0 * matches / alnSize;

  // 1-based coordinates:
  ++alnBeg1;
  ++(strand == '+' ? alnBeg2 : alnEnd2);

  /*
  if( strand == '-' && !isTranslated ){  // xxx this makes it more like BLAST
    std::swap( alnBeg1, alnEnd1 );
    std::swap( alnBeg2, alnEnd2 );
  }
  */

  std::string n1 = seq1.seqName(w1);
  std::string n2 = seq2.seqName(w2);
  FloatText mp("%.2f", matchPercent);
  IntText as(alnSize);
  IntText mm(mismatches);
  IntText go(gapOpens);
  IntText b1(alnBeg1 - seqStart1);
  IntText b2(alnBeg2 - seqStart2);
  IntText e1(alnEnd1 - seqStart1);
  IntText e2(alnEnd2 - seqStart2);
  FloatText ev;
  FloatText bs;
  if( evaluer.isGood() ){
    size_t seqLen2 = seq2.seqLen(w2);
    double area = evaluer.area( score, seqLen2 );
    double epa = evaluer.evaluePerArea( score );
    double bitScore = evaluer.bitScore( score );
    ev.set("%.2g", area * epa);
    bs.set("%.3g", bitScore);
  }

  size_t s =
    n2.size() + n1.size() + mp.size() + as.size() + mm.size() + go.size() +
    b2.size() + e2.size() + b1.size() + e1.size() + 10;
  if (evaluer.isGood()) s += ev.size() + bs.size() + 2;

  char *text = new char[s + 1];
  Writer w(text);
  const char t = '\t';
  w << n2 << t << n1 << t << mp << t << as << t << mm << t << go << t
    << b2 << t << e2 << t << b1 << t << e1;
  if (evaluer.isGood()) w << t << ev << t << bs;
  w << '\n' << '\0';

  AlignmentText at;
  at.queryNum = w2;
  at.score = score;
  at.text = text;
  textAlns.push_back(at);
}

size_t Alignment::numColumns( size_t frameSize ) const{
  size_t num = 0;

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // length of unaligned chunk of top sequence (gaps in bottom sequence):
      num += y.beg1() - x.end1();

      // length of unaligned chunk of bottom sequence (gaps in top sequence):
      size_t gap2, frameshift2;
      sizeAndFrameshift( x.end2(), y.beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) ++num;
      num += gap2;
    }

    num += y.size;  // length of aligned chunk
  }

  return num;
}

static char* writeGaps( char* dest, size_t num ){
  char* end = dest + num;
  while( dest < end ){
    *dest++ = '-';
  }
  return dest;
}

static char* writeQuals( const uchar* qualities, size_t beg, size_t end,
			 size_t qualsPerBase, char* dest ){
  for( size_t i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    *dest++ = *std::max_element( q, q + qualsPerBase );
  }
  return dest;
}

static char* writeSeq( const uchar* seq, size_t beg, size_t end,
		       const Alphabet& alph, size_t qualsPerBase, char* dest ){
  return qualsPerBase ? writeQuals( seq, beg, end, qualsPerBase, dest )
    :                   alph.rtCopy( seq + beg, seq + end, dest );
}

char* Alignment::writeTopSeq( const uchar* seq, const Alphabet& alph,
			      size_t qualsPerBase, size_t frameSize,
			      char* dest ) const{
  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // append unaligned chunk of top sequence:
      dest = writeSeq( seq, x.end1(), y.beg1(), alph, qualsPerBase, dest );

      // append gaps for unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( x.end2(), y.beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) *dest++ = '-';
      dest = writeGaps( dest, gap2 );
    }

    // append aligned chunk of top sequence:
    dest = writeSeq( seq, y.beg1(), y.end1(), alph, qualsPerBase, dest);
  }

  return dest;
}

char* Alignment::writeBotSeq( const uchar* seq, const Alphabet& alph,
			      size_t qualsPerBase, size_t frameSize,
			      char* dest ) const{
  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, y.beg1() - x.end1() );

      //append unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( x.end2(), y.beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 == 1 ) *dest++ = '\\';
      if( frameshift2 == 2 ) *dest++ = '/';
      size_t chunkBeg2 = y.beg2() - gap2;
      dest = writeSeq( seq, chunkBeg2, y.beg2(), alph, qualsPerBase, dest );
    }

    // append aligned chunk of bottom sequence:
    dest = writeSeq( seq, y.beg2(), y.end2(), alph, qualsPerBase, dest );
  }

  return dest;
}
