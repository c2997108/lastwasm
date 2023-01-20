// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SegmentPair.hh"

namespace cbrc {

void SegmentPair::maxIdenticalRun(mcf::BigSeq seq1, const uchar *seq2,
				  const uchar *map1, const uchar *map2) {
  indexT runSize = 0;
  indexT bestSize = 0;
  indexT bestEnd = 0;

  for (indexT i = 0; i < size; ++i) {
    if (map1[seq1[start1 + i]] == map2[seq2[start2 + i]]) {
      ++runSize;
      if (runSize > bestSize) {
	bestSize = runSize;
	bestEnd = i + 1;
      }
    } else {
      runSize = 0;
    }
  }

  start1 += bestEnd - bestSize;
  start2 += bestEnd - bestSize;
  size = bestSize;
}

}
