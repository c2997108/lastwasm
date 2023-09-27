// Author: Martin C. Frith 2023
// SPDX-License-Identifier: GPL-3.0-or-later

// Get and set items, using a given number of bits per item.  The
// items are packed next to each other without wasted bits (apart from
// some wasted bits at the end).

#ifndef MCF_PACKED_ARRAY_HH
#define MCF_PACKED_ARRAY_HH

#include <limits.h>
#include <stddef.h>

namespace mcf {

// Number of bits needed to represent the value n
inline int numOfBitsNeededFor(size_t n) {
  int b = 1;
  while (n >>= 1) ++b;
  return b;
}

// Number of size_t "words" needed to hold the items.  May include one
// extra unused word, because getBits/setBits always use 2 words.
inline size_t numOfWordsNeededFor(int bitsPerItem, size_t numOfItems) {
  const int w = sizeof(size_t) * CHAR_BIT;
  unsigned long long bpi = bitsPerItem;
  return (numOfItems * bpi + w - bpi + w) / w;
}

// Get the i-th item
inline size_t getBits(int bitsPerItem, const size_t *items, size_t i) {
  const int w = sizeof(size_t) * CHAR_BIT;
  const size_t one = 1;
  size_t ones = (one << bitsPerItem) - 1;
  unsigned long long bpi = bitsPerItem;
  size_t q = (i * bpi) / w;
  int    r = (i * bpi) % w;
  int    s = w - r;
  return (items[q+0] >> r | items[q+1] << 1 << (s-1)) & ones;
  // the stupid << 1 << (s-1) avoids undefined behavior when s == w
}

// Set the i-th item
inline void setBits(int bitsPerItem, size_t *items, size_t i, size_t value) {
  const int w = sizeof(size_t) * CHAR_BIT;
  const size_t one = 1;
  size_t ones = (one << bitsPerItem) - 1;
  unsigned long long bpi = bitsPerItem;
  size_t q = (i * bpi) / w;
  int    r = (i * bpi) % w;
  int    s = w - r;
  items[q+0] = (items[q+0] & ~(ones << r         )) | (value << r         );
  items[q+1] = (items[q+1] & ~(ones >> 1 >> (s-1))) | (value >> 1 >> (s-1));
}

// Unpack the items from "packed" into "unpacked" (xxx could be much faster)
inline void unpackBits(int bitsPerItem, const size_t *packed, size_t *unpacked,
		       size_t beg, size_t end) {
  const int w = sizeof(size_t) * CHAR_BIT;
  const size_t one = 1;
  size_t ones = (one << bitsPerItem) - 1;
  unsigned long long bpi = bitsPerItem;
  unsigned long long b = beg * bpi;
  unsigned long long e = end * bpi;
  for (unsigned long long i = b; i < e; i += bpi) {
    size_t q = i / w;
    int    r = i % w;
    int    s = w - r;
    *unpacked++ = (packed[q+0] >> r | packed[q+1] << 1 << (s-1)) & ones;
  }
}

// Pack the items from "unpacked" into "packed" (xxx could be much faster)
inline void packBits(int bitsPerItem, size_t *packed, const size_t *unpacked,
		     size_t beg, size_t end) {
  const int w = sizeof(size_t) * CHAR_BIT;
  const size_t one = 1;
  size_t ones = (one << bitsPerItem) - 1;
  unsigned long long bpi = bitsPerItem;
  unsigned long long b = beg * bpi;
  unsigned long long e = end * bpi;
  for (unsigned long long i = b; i < e; i += bpi) {
    size_t q = i / w;
    int    r = i % w;
    int    s = w - r;
    size_t val = *unpacked++;
    packed[q+0] = (packed[q+0] & ~(ones << r         )) | (val << r         );
    packed[q+1] = (packed[q+1] & ~(ones >> 1 >> (s-1))) | (val >> 1 >> (s-1));
  }
}

struct PackedArray {
  size_t *items;
  unsigned char bitsPerItem;

  size_t operator[](size_t i) const { return getBits(bitsPerItem, items, i); }
  void set(size_t i, size_t x) { setBits(bitsPerItem, items, i, x); }
};

struct ConstPackedArray {
  const size_t *items;
  unsigned char bitsPerItem;

  size_t operator[](size_t i) const { return getBits(bitsPerItem, items, i); }
};

}

#endif
