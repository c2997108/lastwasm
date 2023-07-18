// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

#include "SubsetSuffixArray.hh"
//#include <iostream>  // for debugging

using namespace cbrc;

static size_t posGetAt(const PosPart *p, size_t i) {
  return posGet(p + i * posParts);
}

static size_t offGet(const OffPart *p) {
  size_t x = 0;
  for (int i = 0; i < offParts; ++i) {
    size_t y = p[i];  // must convert to size_t before shifting!
    x += y << (i * sizeof(OffPart) * CHAR_BIT);
  }
  return x;
}

static size_t lowerBound(const PosPart *sufArray, size_t beg, size_t end,
			 const BigPtr &textBase, const uchar *subsetMap,
			 uchar subset) {
  while (beg < end) {
    size_t mid = beg + (end - beg) / 2;
    if (subsetMap[textBase[posGetAt(sufArray, mid)]] < subset) {
      beg = mid + 1;
    } else {
      end = mid;
    }
  }
  return beg;
}

static size_t upperBound(const PosPart *sufArray, size_t beg, size_t end,
			 const BigPtr &textBase, const uchar *subsetMap,
			 uchar subset) {
  while (beg < end) {
    size_t mid = beg + (end - beg) / 2;
    if (subsetMap[textBase[posGetAt(sufArray, mid)]] <= subset) {
      beg = mid + 1;
    } else {
      end = mid;
    }
  }
  return end;
}

// Find the suffix array range of one letter, whose subset is
// "subset", within the suffix array range [beg, end)
static void equalRange(const PosPart *sufArray, size_t &beg, size_t &end,
		       const BigPtr &textBase, const uchar *subsetMap,
		       uchar subset) {
  while (beg < end) {
    size_t mid = beg + (end - beg) / 2;
    uchar s = subsetMap[textBase[posGetAt(sufArray, mid)]];
    if (s < subset) {
      beg = mid + 1;
    } else if (s > subset) {
      end = mid;
    } else {
      beg = lowerBound(sufArray, beg, mid, textBase, subsetMap, subset);
      end = upperBound(sufArray, mid + 1, end, textBase, subsetMap, subset);
      return;
    }
  }
}

// Same as the 1st equalRange, but uses more info and may be faster
static void equalRange(const PosPart *sufArray, size_t &beg, size_t &end,
		       const BigPtr &textBase, const uchar *subsetMap,
		       uchar subset, uchar begSubset, uchar endSubset,
		       size_t begOffset, size_t endOffset) {
  size_t b = beg + begOffset;
  size_t e = end - endOffset;
  if (subset == begSubset) {
    end = upperBound(sufArray, b, e, textBase, subsetMap, subset);
  } else if (subset == endSubset) {
    beg = lowerBound(sufArray, b, e, textBase, subsetMap, subset);
  } else {
    beg = b;
    end = e;
    equalRange(sufArray, beg, end, textBase, subsetMap, subset);
  }
}

// Same as the 1st equalRange, but tries to be faster by checking endpoints
static void fastEqualRange(const PosPart *sufArray, size_t &beg, size_t &end,
			   const BigPtr &textBase, const uchar *subsetMap,
			   uchar subset) {
  uchar b = subsetMap[textBase[posGetAt(sufArray, beg)]];
  if (subset < b) { end = beg; return; }
  uchar e = subsetMap[textBase[posGetAt(sufArray, end - 1)]];
  if (subset > e) { beg = end; return; }
  if (b == e) return;
  equalRange(sufArray, beg, end, textBase, subsetMap, subset, b, e, 1, 1);
}

static size_t lowerBound2(const PosPart *sufArray, size_t beg, size_t end,
			  BigSeq text, size_t depth, const uchar *subsetMap,
			  const uchar *queryBeg, const uchar *queryEnd,
			  const CyclicSubsetSeed &seed) {
  while (beg < end) {
    size_t mid = beg + (end - beg) / 2;
    size_t offset = posGetAt(sufArray, mid);
    size_t t = depth + offset;
    const uchar *q = queryBeg;
    const uchar *s = subsetMap;
    for (;;) {  // loop over consecutive letters
      const uchar *textSubsetMap = seed.originalSubsetMap(s);
      if (textSubsetMap[text[t]] < s[*q]) {
	beg = mid + 1;
	// the next 3 lines are unnecessary, but make it faster:
	queryBeg = q;
	depth = t - offset;
	subsetMap = s;
	break;
      }
      ++q;  // next query letter
      if (q == queryEnd) {  // we found a full match to [queryBeg, queryEnd)
	end = mid;
	break;
      }
      ++t;  // next text letter
      s = seed.nextMap(s);  // next mapping from letters to subsets
    }
  }
  return beg;
}

static size_t upperBound2(const PosPart *sufArray, size_t beg, size_t end,
			  BigSeq text, size_t depth, const uchar *subsetMap,
			  const uchar *queryBeg, const uchar *queryEnd,
			  const CyclicSubsetSeed &seed) {
  while (beg < end) {
    size_t mid = beg + (end - beg) / 2;
    size_t offset = posGetAt(sufArray, mid);
    size_t t = depth + offset;
    const uchar *q = queryBeg;
    const uchar *s = subsetMap;
    for (;;) {  // loop over consecutive letters
      const uchar *textSubsetMap = seed.originalSubsetMap(s);
      if (textSubsetMap[text[t]] > s[*q]) {
	end = mid;
	// the next 3 lines are unnecessary, but make it faster:
	queryBeg = q;
	depth = t - offset;
	subsetMap = s;
        break;
      }
      ++q;  // next query letter
      if (q == queryEnd) {  // we found a full match to [queryBeg, queryEnd)
	beg = mid + 1;
	break;
      }
      ++t;  // next text letter
      s = seed.nextMap(s);  // next mapping from letters to subsets
    }
  }
  return end;
}

// Find the suffix array range of string [queryBeg, queryEnd) within
// the suffix array range [beg, end)
static void equalRange2(const PosPart *sufArray, size_t &beg, size_t &end,
			BigSeq text, size_t depth, const uchar *subsetMap,
			const uchar *queryBeg, const uchar *queryEnd,
			const CyclicSubsetSeed &seed) {
  const uchar *qBeg = queryBeg;
  const uchar *qEnd = qBeg;
  size_t tBeg = depth;
  size_t tEnd = tBeg;
  const uchar *sBeg = subsetMap;
  const uchar *sEnd = sBeg;

  while (beg < end) {
    size_t mid = beg + (end - beg) / 2;
    size_t offset = posGetAt(sufArray, mid);
    const uchar *q;
    size_t t;
    const uchar *s;
    if (qBeg < qEnd) {
      q = qBeg;
      t = tBeg + offset;
      s = sBeg;
    } else {
      q = qEnd;
      t = tEnd + offset;
      s = sEnd;
    }
    uchar x, y;
    for (;;) {  // loop over consecutive letters
      const uchar *textSubsetMap = seed.originalSubsetMap(s);
      x = textSubsetMap[text[t]];  // this text letter's subset
      y = s[*q];  // this query letter's subset
      if (x != y) break;
      ++q;  // next query letter
      if (q == queryEnd) {  // we found a full match to [queryBeg, queryEnd)
	beg = lowerBound2(sufArray, beg, mid, text, tBeg, sBeg,
			  qBeg, queryEnd, seed);
	end = upperBound2(sufArray, mid + 1, end, text, tEnd, sEnd,
			  qEnd, queryEnd, seed);
	return;
      }
      ++t;  // next text letter
      s = seed.nextMap(s);  // next mapping from letters to subsets
    }
    if (x < y) {
      beg = mid + 1;
      // the next 3 lines are unnecessary, but make it faster:
      qBeg = q;
      tBeg = t - offset;
      sBeg = s;
    } else {
      end = mid;
      // the next 3 lines are unnecessary, but make it faster:
      qEnd = q;
      tEnd = t - offset;
      sEnd = s;
    }
  }
}

// Find the suffix array range of string [queryBeg, queryBeg+d) within
// the suffix array range [beg, end), where d is not known in advance.
// The routine tries to find a large-as-possible d: it guarantees that
// the suffix array range for d-1 is longer than maxHits.  d is
// returned.  Actually, this routine may find *part* of the SA range
// of [queryBeg, queryBeg+d), which is guaranteed to include the whole
// range for the smallest d whose range is no longer than maxHits.
static size_t equalRange3(const PosPart *sufArray, size_t &beg, size_t &end,
			  const uchar *&subsetMap, BigSeq text, size_t depth,
			  const uchar *queryBeg, const CyclicSubsetSeed &seed,
			  size_t maxHits) {
  if (subsetMap[*queryBeg] == CyclicSubsetSeed::DELIMITER) return 0;
  assert(end - beg > maxHits);
  const uchar *qBeg = queryBeg;
  size_t tBeg = depth;
  const uchar *sBeg = subsetMap;
  const uchar *qBegOld = qBeg;
  size_t tBegOld = tBeg;
  const uchar *sBegOld = sBeg;
  const uchar *qEnd = queryBeg;
  size_t tEnd = depth;
  const uchar *sEnd = subsetMap;
  const uchar *qEndOld = qEnd;
  size_t tEndOld = tEnd;
  const uchar *sEndOld = sEnd;
  const uchar *qMid = queryBeg;
  size_t tMid = depth;
  const uchar *sMid = subsetMap;

  while (end - beg > maxHits * 2) {
    size_t mid = beg + (end - beg) / 2;
    size_t offset = posGetAt(sufArray, mid);
    tMid += offset;
    int iterations = 1023;  // xxx ???
    uchar tChar, qChar;
    for (;;) {  // loop over consecutive letters
      const uchar *textSubsetMap = seed.originalSubsetMap(sMid);
      tChar = textSubsetMap[text[tMid]];  // this text letter's subset
      qChar = sMid[*qMid];  // this query letter's subset
      if (tChar != qChar || qChar == CyclicSubsetSeed::DELIMITER) break;
      if (--iterations == 0) {  // avoid huge self-comparisons
	subsetMap = qBeg < qEnd ? sBeg : sEnd;
	return std::min(qBeg, qEnd) - queryBeg;
      }
      ++qMid;  // next query letter
      ++tMid;  // next text letter
      sMid = seed.nextMap(sMid);  // next mapping from letters to subsets
    }
    tMid -= offset;
    if (tChar <= qChar) {
      beg = mid + 1;
      qBegOld = qBeg;
      tBegOld = tBeg;
      sBegOld = sBeg;
      qBeg = qMid;
      tBeg = tMid;
      sBeg = sMid;
      if (qEnd < qMid) {
	if (qChar == CyclicSubsetSeed::DELIMITER && qBeg == qBegOld) {
	  beg = end;
	  return 0;
	}
	qMid = qEnd;
	tMid = tEnd;
	sMid = sEnd;
      }
    } else {
      end = mid;
      qEndOld = qEnd;
      tEndOld = tEnd;
      sEndOld = sEnd;
      qEnd = qMid;
      tEnd = tMid;
      sEnd = sMid;
      if (qBeg < qMid) {
	qMid = qBeg;
	tMid = tBeg;
	sMid = sBeg;
      }
    }
  }

  if (qBeg < qEnd) {
    qMid = std::max(qBeg, qEndOld) + 1;
    subsetMap = seed.nextMap(qBeg > qEndOld ? sBeg : sEndOld);
    if (qMid > qEnd) {
      equalRange2(sufArray, beg, end, text, tBeg, sBeg, qBeg, qMid, seed);
    } else {
      beg = lowerBound2(sufArray, beg, end,
			text, tBeg, sBeg, qBeg, qMid, seed);
      end = upperBound2(sufArray, end + 1, end + maxHits + 1,
			text, tEndOld, sEndOld, qEndOld, qMid, seed);
    }
  } else {
    qMid = std::max(qEnd, qBegOld) + 1;
    subsetMap = seed.nextMap(qEnd > qBegOld ? sEnd : sBegOld);
    if (qMid > qBeg) {
      equalRange2(sufArray, beg, end, text, tEnd, sEnd, qEnd, qMid, seed);
    } else {
      beg = lowerBound2(sufArray, beg - maxHits - 1, beg - 1,
			text, tBegOld, sBegOld, qBegOld, qMid, seed);
      end = upperBound2(sufArray, beg, end,
			text, tEnd, sEnd, qEnd, qMid, seed);
    }
  }

  return qMid - queryBeg;
}

// use past results to speed up long matches?
// could & probably should return the match depth
void SubsetSuffixArray::match(const PosPart *&begPtr, const PosPart *&endPtr,
			      const uchar *queryPtr, BigSeq text,
			      unsigned seedNum, size_t maxHits,
			      size_t minDepth, size_t maxDepth) const {
  // the next line is unnecessary, but makes it faster in some cases:
  if( maxHits == 0 && minDepth < maxDepth ) minDepth = maxDepth;

  size_t depth = 0;
  const CyclicSubsetSeed &seed = seeds[seedNum];
  const uchar* subsetMap = seed.firstMap();

  // match using buckets:
  size_t bucketDepth = maxBucketPrefix(seedNum);
  size_t startDepth = std::min( bucketDepth, maxDepth );
  const OffPart *bucketPtr = bucketEnds[seedNum];
  const size_t *myBucketSteps = bucketStepEnds[seedNum];

  while( depth < startDepth ){
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) break;
    ++depth;
    bucketPtr += subset * myBucketSteps[depth];
    subsetMap = seed.nextMap( subsetMap );
  }

  size_t beg = offGet(bucketPtr);
  size_t end = offGet(bucketPtr + myBucketSteps[depth]);

  while( depth > minDepth && end - beg < maxHits ){
    // maybe we lengthened the match too far: try shortening it again
    const uchar* oldMap = seed.prevMap( subsetMap );
    uchar subset = oldMap[ queryPtr[depth-1] ];
    bucketPtr -= subset * myBucketSteps[depth];
    size_t oldBeg = offGet(bucketPtr);
    size_t oldEnd = offGet(bucketPtr + myBucketSteps[depth-1]);
    if( oldEnd - oldBeg > maxHits ) break;
    subsetMap = oldMap;
    beg = oldBeg;
    end = oldEnd;
    --depth;
  }

  const PosPart *sufArray = suffixArray.begin();

  // match using binary search:

  if( depth < minDepth ){
    size_t d = depth;
    const uchar* s = subsetMap;
    while( depth < minDepth ){
      uchar subset = subsetMap[ queryPtr[depth] ];
      if( subset == CyclicSubsetSeed::DELIMITER ){
	beg = end;
	break;
      }
      ++depth;
      subsetMap = seed.nextMap( subsetMap );
    }
    equalRange2(sufArray, beg, end, text, d, s,
		queryPtr + d, queryPtr + depth, seed);
  }

  if (end - beg > maxHits * 2 && maxDepth + 1 == 0 &&
      childTable.empty() && kiddyTable.empty() && chibiTable.empty()) {
    depth += equalRange3(sufArray, beg, end, subsetMap,
			 text, depth, queryPtr + depth, seed, maxHits);
  }

  ChildDirection childDirection = UNKNOWN;

  while( end - beg > maxHits && depth < maxDepth ){
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ){
      beg = end;
      break;
    }
    const uchar *textSubsetMap = seed.originalSubsetMap(subsetMap);
    childRange(beg, end, childDirection, text + depth, textSubsetMap, subset);
    ++depth;
    subsetMap = seed.nextMap( subsetMap );
  }

  begPtr = sufArray + beg * posParts;
  endPtr = sufArray + end * posParts;
}

void SubsetSuffixArray::countMatches(std::vector<unsigned long long> &counts,
				     const uchar *queryPtr, BigSeq text,
				     unsigned seedNum, size_t maxDepth) const {
  size_t depth = 0;
  const CyclicSubsetSeed &seed = seeds[seedNum];
  const uchar* subsetMap = seed.firstMap();

  // match using buckets:
  size_t bucketDepth = maxBucketPrefix(seedNum);
  const OffPart *bucketPtr = bucketEnds[seedNum];
  const size_t *myBucketSteps = bucketStepEnds[seedNum];
  size_t beg = offGet(bucketPtr);
  size_t end = offGet(bucketPtr + myBucketSteps[depth]);

  while( depth < bucketDepth ){
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    if( depth >= maxDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    ++depth;
    size_t step = myBucketSteps[depth];
    bucketPtr += subset * step;
    beg = offGet(bucketPtr);
    end = offGet(bucketPtr + step);
    subsetMap = seed.nextMap( subsetMap );
  }

  // match using binary search:
  ChildDirection childDirection = UNKNOWN;

  while( beg < end ){
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    if( depth >= maxDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    const uchar *textSubsetMap = seed.originalSubsetMap(subsetMap);
    childRange(beg, end, childDirection, text + depth, textSubsetMap, subset);
    ++depth;
    subsetMap = seed.nextMap( subsetMap );
  }
}

void SubsetSuffixArray::childRange(size_t &beg, size_t &end,
				   ChildDirection &childDirection,
				   BigPtr textBase, const uchar *subsetMap,
				   uchar subset) const {
  const PosPart *sufArray = suffixArray.begin();

  if( childDirection == UNKNOWN ){
    size_t mid = getChildForward( beg );
    if( mid == beg ){  // failure: never happens with the full childTable
      mid = getChildReverse( end );
      if( mid == end ){  // failure: never happens with the full childTable
	fastEqualRange(sufArray, beg, end, textBase, subsetMap, subset);
	return;
      }
      childDirection = REVERSE;
    }else{
      childDirection = (mid < end) ? FORWARD : REVERSE;
    }
  }

  if( childDirection == FORWARD ){
    uchar e = subsetMap[textBase[posGetAt(sufArray, end - 1)]];
    if( subset > e ){ beg = end; return; }
    if( subset < e ) childDirection = REVERSE;  // flip it for next time
    while( 1 ){
      uchar b = subsetMap[textBase[posGetAt(sufArray, beg)]];
      if( subset < b ) { end = beg; return; }
      if( b == e ) return;
      size_t mid = getChildForward( beg );
      if( mid == beg ){  // failure: never happens with the full childTable
	size_t offset = kiddyTable.empty() ? UCHAR_MAX : USHRT_MAX;
	equalRange(sufArray, beg, end, textBase, subsetMap, subset,
		   b, e, offset, 1);
	return;
      }
      if( subset == b ) { end = mid; return; }
      beg = mid;
      if( b + 1 == e ) return;  // unnecessary, but may be faster
    }
  }else{
    uchar b = subsetMap[textBase[posGetAt(sufArray, beg)]];
    if( subset < b ) { end = beg; return; }
    if( subset > b ) childDirection = FORWARD;  // flip it for next time
    while( 1 ){
      uchar e = subsetMap[textBase[posGetAt(sufArray, end - 1)]];
      if( subset > e ){ beg = end; return; }
      if( b == e ) return;
      size_t mid = getChildReverse( end );
      if( mid == end ){  // failure: never happens with the full childTable
	size_t offset = kiddyTable.empty() ? UCHAR_MAX : USHRT_MAX;
	equalRange(sufArray, beg, end, textBase, subsetMap, subset,
		   b, e, 1, offset);
	return;
      }
      if( subset == e ) { beg = mid; return; }
      end = mid;
      if( b + 1 == e ) return;  // unnecessary, but may be faster
    }
  }
}
