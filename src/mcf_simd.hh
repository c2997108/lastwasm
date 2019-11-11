// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MCF_SIMD_HH
#define MCF_SIMD_HH

#include <immintrin.h>

namespace mcf {

//#if defined __AVX2__
#if defined AVX2_SEEMS_SLOW

typedef __m256i SimdInt;

const int simdLen = 8;

static inline SimdInt simdSet(int i7, int i6, int i5, int i4,
			      int i3, int i2, int i1, int i0) {
  return _mm256_set_epi32(i7, i6, i5, i4, i3, i2, i1, i0);
}

static inline SimdInt simdSet1(int x) {
  return _mm256_set1_epi32(x);
}

static inline SimdInt simdLoad(const int *p) {
  return _mm256_loadu_si256((const SimdInt *)p);
}

static inline void simdStore(int *p, SimdInt x) {
  _mm256_storeu_si256((SimdInt *)p, x);
}

static inline SimdInt simdGt(SimdInt x, SimdInt y) {
  return _mm256_cmpgt_epi32(x, y);
}

static inline SimdInt simdAdd(SimdInt x, SimdInt y) {
  return _mm256_add_epi32(x, y);
}

static inline SimdInt simdSub(SimdInt x, SimdInt y) {
  return _mm256_sub_epi32(x, y);
}

static inline SimdInt simdMax(SimdInt x, SimdInt y) {
  return _mm256_max_epi32(x, y);
}

static inline SimdInt simdBlend(SimdInt x, SimdInt y, SimdInt mask) {
  return _mm256_blendv_epi8(x, y, mask);
}

static inline int simdHorizontalMax(SimdInt x) {
  x = simdMax(x, _mm256_permute4x64_epi64(x, 0x4E));
  x = simdMax(x, _mm256_shuffle_epi32(x, 0x4E));
  x = simdMax(x, _mm256_shuffle_epi32(x, 0xB1));
  //return _mm256_cvtsi256_si32(x);  // doesn't work with gcc 4.8.5?
  return _mm256_extract_epi32(x, 0);
}

#elif defined __SSE4_1__

typedef __m128i SimdInt;

const int simdLen = 4;

static inline SimdInt simdSet(int i3, int i2, int i1, int i0) {
  return _mm_set_epi32(i3, i2, i1, i0);
}

static inline SimdInt simdSet1(int x) {
  return _mm_set1_epi32(x);
}

static inline SimdInt simdLoad(const int *p) {
  return _mm_loadu_si128((const SimdInt *)p);
}

static inline void simdStore(int *p, SimdInt x) {
  _mm_storeu_si128((SimdInt *)p, x);
}

static inline SimdInt simdGt(SimdInt x, SimdInt y) {
  return _mm_cmpgt_epi32(x, y);
}

static inline SimdInt simdAdd(SimdInt x, SimdInt y) {
  return _mm_add_epi32(x, y);
}

static inline SimdInt simdSub(SimdInt x, SimdInt y) {
  return _mm_sub_epi32(x, y);
}

static inline SimdInt simdMax(SimdInt x, SimdInt y) {
  return _mm_max_epi32(x, y);  // SSE4.1
}

static inline SimdInt simdBlend(SimdInt x, SimdInt y, SimdInt mask) {
  return _mm_blendv_epi8(x, y, mask);  // SSE4.1
}

static inline int simdHorizontalMax(SimdInt x) {
  x = simdMax(x, _mm_shuffle_epi32(x, 0x4E));
  x = simdMax(x, _mm_shuffle_epi32(x, 0xB1));
  return _mm_cvtsi128_si32(x);
}

#else

typedef int SimdInt;
const int simdLen = 1;
static inline int simdSet(int x) { return x; }
static inline int simdSet1(int x) { return x; }
static inline int simdLoad(const int *p) { return *p; }
static inline void simdStore(int *p, int x) { *p = x; }
static inline int simdGt(int x, int y) { return x > y; }
static inline int simdAdd(int x, int y) { return x + y; }
static inline int simdSub(int x, int y) { return x - y; }
static inline int simdMax(int x, int y) { return x > y ? x : y; }
static inline int simdBlend(int x, int y, int mask) { return mask ? y : x; }
static inline int simdHorizontalMax(int a) { return a; }

#endif

}

#endif
