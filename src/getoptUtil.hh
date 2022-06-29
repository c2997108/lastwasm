// Author: Martin C. Frith 2017
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef GETOPT_UTIL_HH
#define GETOPT_UTIL_HH

#include <getopt.h>

inline void resetGetopt() {
  optind = 1;  // xxx ???
}

#endif
