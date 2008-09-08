// Copyright 2008 Martin C. Frith

#include "GeneralizedAffineGapCosts.hh"
#include <algorithm>

namespace cbrc{

int GeneralizedAffineGapCosts::cost( int gapSize1, int gapSize2 ) const{
  int sizeMin = std::min( gapSize1, gapSize2 );
  int sizeMax = std::max( gapSize1, gapSize2 );
  int sizeDif = sizeMax - sizeMin;

  int cost1 = exist + extend * sizeDif + extendPair * sizeMin;
  int cost2 = exist * 2 + extend * ( gapSize1 + gapSize2 );

  return std::min( cost1, cost2 );
}

}  // end namespace cbrc
