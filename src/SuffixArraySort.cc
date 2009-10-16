// Copyright 2008, 2009 Martin C. Frith

// Parts of this code are adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#include "SuffixArray.hh"
#include "PeriodicSpacedSeed.hh"

using namespace cbrc;

namespace{
  typedef unsigned indexT;
  typedef unsigned char uchar;
  struct Stack{ indexT* b; indexT* e; indexT d; const uchar* t; };
  Stack stack[1048576];  // big enough???
  Stack* sp = stack;
}

#define PUSH(b1, e1, d1, t1) sp->b = b1, sp->e = e1, sp->d = d1, (sp++)->t = t1
#define  POP(b1, e1, d1, t1) b1 = (--sp)->b, e1 = sp->e, d1 = sp->d, t1 = sp->t

void SuffixArray::sortIndex( const uchar* text,
			     const PeriodicSpacedSeed& seed ){
  PUSH( &index[0], &index[0] + index.size(), 0, text );

  while( sp > stack ){
    indexT* beg;
    indexT* end;
    indexT depth;
    const uchar* textBase;
    POP( beg, end, depth, textBase );

    if( end - beg < 10 ){  // 10 seems good for hg18 chr21
      insertionSort( seed, beg, end, depth, textBase );
      continue;
    }

    // get the next textBase, and increment depth, for sorting within buckets:
    const uchar* textNext = textBase + seed.offsets[ depth % seed.weight ];
    ++depth;

    if( alphSize == 4 )
      radixSortAlph4( textBase, textNext, beg, end, depth );
    else
      radixSort( textBase, textNext, beg, end, depth );
  }
}

void SuffixArray::radixSort( const uchar* textBase, const uchar* textNext,
			     indexT* beg, indexT* end, indexT depth ){
  static indexT bucketSize[256];  // initialized to zero at startup
  /*  */ indexT* bucketEnd[256];  // "static" makes little difference to speed

  // get bucket sizes (i.e. letter counts):
  // The intermediate oracle array makes it faster (see "Engineering
  // Radix Sort for Strings" by J Karkkainen & T Rantala)
  for( indexT* i = beg; i < end; /* noop */ ){
    uchar oracle[256];
    uchar* oracleEnd =
      oracle + std::min( sizeof(oracle), std::size_t(end - i) );
    for( uchar* j = oracle; j < oracleEnd; ++j )
      *j = textBase[ *i++ ];
    for( uchar* j = oracle; j < oracleEnd; ++j )
      ++bucketSize[ *j ];
  }

  // get bucket ends, and put buckets on the stack to sort within them later:
  // (could push biggest bucket first, to ensure logarithmic stack growth)
  indexT* pos = beg;
  for( unsigned i = 0; i < alphSize; ++i ){
    indexT* nextPos = pos + bucketSize[i];
    PUSH( pos, nextPos, depth, textNext );
    pos = nextPos;
    bucketEnd[i] = pos;
  }
  bucketEnd[alphSize] = end;  // don't sort within the delimiter bucket

  // permute items into the correct buckets:
  for( indexT* i = beg; i < end; /* noop */ ) {
    unsigned symbol;  // unsigned is faster than uchar!
    indexT holdOut = *i;
    while( --bucketEnd[ symbol = textBase[holdOut] ] > i ){
      std::swap( *bucketEnd[symbol], holdOut );
    }
    *i = holdOut;
    i += bucketSize[symbol];
    bucketSize[symbol] = 0;  // reset it so we can reuse it
  }
}

// Specialized radix sort for alphSize = 4 (plus one delimiter symbol)
void SuffixArray::radixSortAlph4( const uchar* textBase, const uchar* textNext,
				  indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* end2 = beg;  // end of '2's
  indexT* beg3 = end;  // beginning of '3's

  while( end2 < beg3 ){
    const indexT x = *end2;
    switch( textBase[x] ){
    case 0:
      *end2++ = *end1;
      *end1++ = *end0;
      *end0++ = x;
      break;
    case 1:
      *end2++ = *end1;
      *end1++ = x;
      break;
    case 2:
      end2++;
      break;
    case 3:
      *end2 = *--beg3;
      *beg3 = x;
      break;
    default:  // it's a delimiter
      *end2 = *--beg3;
      *beg3 = *--end;
      *end = x;
      break;
    }
  }

  // put buckets on the stack to sort within them later:
  // (could push biggest bucket first, to ensure logarithmic stack growth)
  PUSH( beg, end0, depth, textNext );   // the '0's
  PUSH( end0, end1, depth, textNext );  // the '1's
  PUSH( end1, beg3, depth, textNext );  // the '2's
  PUSH( beg3, end, depth, textNext );  // the '3's
  // don't sort within the delimiter bucket
}

void SuffixArray::insertionSort( const PeriodicSpacedSeed& seed,
				 indexT* beg, indexT* end,
				 indexT depth, const uchar* textBase ){
  if( seed.weight == 1 && seed.offsets[0] == 1 ){
    insertionSortSimple( beg, end, textBase );  // how much faster?
    return;
  }

  const indexT* maskPtr = &seed.offsets[ depth % seed.weight ];
  const indexT* maskEnd = &seed.offsets[0] + seed.weight;

  for( indexT* i = beg+1; i < end; ++i ){
    for( indexT* j = i; j > beg; --j ){
      const uchar* s = textBase + *(j-1);
      const uchar* t = textBase + *j;
      const indexT* m = maskPtr;

      while( *s == *t && *s < alphSize ){
	indexT off = *m;
	s += off;
	t += off;
	++m;
	if( m == maskEnd ) m = &seed.offsets[0];  // seems faster than modulus
      }

      if( *s <= *t ) break;
      std::iter_swap( j, j-1 );
    }
  }
}

void SuffixArray::insertionSortSimple( indexT* beg, indexT* end,
				       const uchar* textBase ){
  for( indexT* i = beg+1; i < end; ++i ){
    for( indexT* j = i; j > beg; --j ){
      const uchar* s = textBase + *(j-1);
      const uchar* t = textBase + *j;
      while( *s == *t && *s < alphSize ){
	++s;
	++t;
      }
      if( *s <= *t ) break;
      std::iter_swap( j, j-1 );
    }
  }
}
