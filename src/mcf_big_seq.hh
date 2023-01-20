// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MCF_BIG_SEQ
#define MCF_BIG_SEQ

#include <stddef.h>

namespace mcf {

// A pointer to the start of a sequence, where each element is stored
// in either 4 bits or 1 byte
struct BigSeq {
  const unsigned char *beg;
  bool is4bit;

  static int from4bit(const unsigned char *b, size_t i) {
    return (b[i >> 1] >> ((i & 1) << 2)) & 15;
  }

  int operator[](size_t i) const {
    return is4bit ? from4bit(beg, i) : beg[i];
  }
};

// A pointer to some position in a sequence, where each element is
// stored in either 4 bits or 1 byte
struct BigPtr {
  const unsigned char *beg;
  size_t pos;
  bool is4bit;

  int operator[](size_t i) const {
    return is4bit ? BigSeq::from4bit(beg, pos + i) : beg[i];
  }

  int operator*() const {
    return is4bit ? BigSeq::from4bit(beg, pos) : beg[pos];
  }

  BigPtr &operator+=(int i) {
    pos += i;
    return *this;
  }
};

inline int getNext(BigPtr &x) {
  return x.is4bit ? BigSeq::from4bit(x.beg, x.pos++) : *x.beg++;
}

inline int getPrev(BigPtr &x) {
  return x.is4bit ? BigSeq::from4bit(x.beg, --x.pos) : *--x.beg;
}

inline BigPtr operator+(BigSeq s, size_t i) {
  BigPtr p = {s.beg + i * !s.is4bit, i * s.is4bit, s.is4bit};
  return p;
}

}

#endif
