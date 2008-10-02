// Copyright 2008 Michiaki Hamada

#include "LambdaCalculator.hh"
#include <vector>
#include <cassert>

namespace lambda
{
extern "C" {
#include "CA_code/lambda_calculator.h"
}
}

namespace cbrc{
  double LambdaCalculator::calculate( const int matrix[MAT][MAT],
				      int alphSize ){
    assert( alphSize < MAT );

    std::vector< double > x( alphSize * alphSize + 1 );
    std::vector< const double* > y( alphSize + 1 );

    for( int i = 0; i < alphSize; ++i ){
      y[ i+1 ] = &x[ i * alphSize ];
      for( int j = 0; j < alphSize; ++j ){
	x[ i * alphSize + j + 1 ] = matrix[i][j];
      }
    }

    return lambda::calculate_lambda( &y[0], alphSize );
  }
}
