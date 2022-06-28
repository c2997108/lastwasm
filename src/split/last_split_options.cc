// Author: Martin C. Frith 2022
// SPDX-License-Identifier: GPL-3.0-or-later

#include "last_split_options.hh"

#include <ctype.h>
#include <math.h>

#include <iostream>

#define OPT_d 1
#define OPT_c 0.004
#define OPT_t 1e-05
#define OPT_M 7.0
#define OPT_S 1.7
#define OPT_m 1.0

LastSplitOptions::LastSplitOptions()
  : format(0),
    isTopSeqQuery(false),
    direction(OPT_d),
    cis(OPT_c),
    trans(OPT_t),
    mean(OPT_M),
    sdev(OPT_S),
    mismap(OPT_m),
    score(-1),
    no_split(false),
    bytes(0),
    verbose(false),
    isSplicedAlignment(false) {}

static size_t defaultBytes(bool isSplicedAlignment) {
  size_t b = isSplicedAlignment ? 8 : 8 * 1024;
  for (int i = 0; i < 3; ++i) {
    size_t n = b * 1024;
    if (n / 1024 != b) return -1;
    b = n;
  }
  return b;
}

void LastSplitOptions::setUnspecifiedValues(int lastalMinScore, double scale) {
  if (!bytes) bytes = defaultBytes(isSplicedAlignment);

  if (score < 0) {
    score = lastalMinScore;
    if (isSplicedAlignment) score += floor(scale * log(100.0) + 0.5);
  }
}

void LastSplitOptions::print() const {
  std::streamsize p = std::cout.precision(12);
  std::cout << '#'
	    << " m=" << mismap
	    << " s=" << score;
  if (isSplicedAlignment) {
    std::cout << " d=" << direction
	      << " c=" << cis
	      << " t=" << trans
	      << " M=" << mean
	      << " S=" << sdev;
  }
  std::cout << '\n';
  std::cout.precision(p);
}

char LastSplitOptions::parseOutputFormat(const char *text) {
  std::string s = text;
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = tolower(s[i]);
  }
  if (s == "maf")  return 'm';
  if (s == "maf+") return 'M';
  return 0;
}
