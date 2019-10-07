// Copyright 2010, 2011 Martin C. Frith

// These definitions are for manipulating a score matrix, or a PSSM,
// as a 2-dimensional array.  The hope is that 2D arrays are extremely
// fast.  But arrays are "evil", so this may not be the best way...

#ifndef SCORE_MATRIX_ROW_HH
#define SCORE_MATRIX_ROW_HH

#include <climits>  // INT_MAX

namespace cbrc{

enum { scoreMatrixRowSize = ALPHABET_CAPACITY };

typedef int ScoreMatrixRow[scoreMatrixRowSize];

// An "infinite" score.  Delimiters at the ends of sequences get a
// score of -INF.  We want it high enough to terminate alignments
// immediately, but not so high that it causes overflow errors.
enum { INF = INT_MAX / 2 };

}

#endif
