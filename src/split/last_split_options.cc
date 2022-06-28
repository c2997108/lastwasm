// Author: Martin C. Frith 2022
// SPDX-License-Identifier: GPL-3.0-or-later

#include "last_split_options.hh"

#include <ctype.h>

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

char LastSplitOptions::parseOutputFormat(const char *text) {
  std::string s = text;
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = tolower(s[i]);
  }
  if (s == "maf")  return 'm';
  if (s == "maf+") return 'M';
  return 0;
}
