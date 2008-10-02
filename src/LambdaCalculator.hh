// Copyright 2008 Michiaki Hamada

#ifndef LAMBDACALCULATOR_HH
#define LAMBDACALCULATOR_HH

namespace cbrc{
  class LambdaCalculator{
  public:
    enum { MAT = 64 };

    // compute $\lambda$ value for a given score matrix.
    // return -1.0 if the score matrix is invalid.
    static double calculate( const int matrix[MAT][MAT], int alphSize );
  };
}
#endif // LAMBDACALCULATOR_HH
