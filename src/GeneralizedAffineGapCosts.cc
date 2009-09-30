// Copyright 2008, 2009 Martin C. Frith

#include "GeneralizedAffineGapCosts.hh"
#include <algorithm>
#include <cassert>

namespace cbrc{

int GeneralizedAffineGapCosts::cost( int gapSize1, int gapSize2 ) const{
  int sizeMin = std::min( gapSize1, gapSize2 );
  int sizeMax = std::max( gapSize1, gapSize2 );
  int sizeDif = sizeMax - sizeMin;

  if( sizeMax == 0 ) return 0;

  int cost1 = exist + extend * sizeDif + extendPair * sizeMin;
  int cost2 = exist * 2 + extend * ( gapSize1 + gapSize2 );

  assert( cost1 >= exist );  // try to catch overflow errors

  return std::min( cost1, cost2 );
}

}  // end namespace cbrc
