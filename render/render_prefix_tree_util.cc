/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <new>
#include <iostream>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <math.h> // M_PI
#else
#include <cmath>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <limits>
#include <algorithm>

#include "render_common.h"
#include "render_timerutil.h"
#include "render_prefix_tree_util.h"
#include "render_simd_util.h"

using namespace lsgl::render;

#if defined(__SSE2__) || (_M_IX86_FP >= 2)
namespace {

inline __m128i vpartby2(const int *nn) {
  const __m128i nn_ld_1 = _mm_loadu_si128((__m128i const *)(nn + 0));
  const __m128i t_ivec1_1 = nn_ld_1;
  __m128i n_1 = t_ivec1_1;

  const __m128i t_ivec2_1 = n_1;
  const __m128i t_ivec3_1 = n_1;
  const int t_int5_1 = 16U;
  const int t_int4_1 = t_int5_1;
  const __m128i t_ivec7_1 = _mm_slli_epi32((__m128i)t_ivec3_1, t_int4_1);
  const __m128i t_ivec8_1 = _mm_xor_si128(t_ivec2_1, t_ivec7_1);
  const int t_int9_1 = 4278190335U;
  const int t_int10_1 = t_int9_1;
  const __m128i t_ivec11_1 = _mm_set1_epi32(t_int10_1);
  const __m128i t_ivec12_1 = (__m128i)(t_ivec11_1);
  const __m128i t_ivec13_1 = _mm_and_si128(t_ivec8_1, t_ivec12_1);
  n_1 = t_ivec13_1;

  const __m128i t_ivec14_1 = n_1;
  const __m128i t_ivec15_1 = n_1;
  const int t_int16_1 = 8U;
  const int t_int17_1 = t_int16_1;
  const __m128i t_ivec19_1 = _mm_slli_epi32((__m128i)t_ivec15_1, t_int17_1);
  const __m128i t_ivec20_1 = _mm_xor_si128(t_ivec14_1, t_ivec19_1);
  const int t_int21_1 = 50393103U;
  const int t_int22_1 = t_int21_1;
  const __m128i t_ivec23_1 = _mm_set1_epi32(t_int22_1);
  const __m128i t_ivec24_1 = (__m128i)(t_ivec23_1);
  const __m128i t_ivec25_1 = _mm_and_si128(t_ivec20_1, t_ivec24_1);
  n_1 = t_ivec25_1;

  const __m128i t_ivec26_1 = n_1;
  const __m128i t_ivec27_1 = n_1;
  const int t_int28_1 = 4U;
  const int t_int29_1 = t_int28_1;
  const __m128i t_ivec31_1 = _mm_slli_epi32((__m128i)t_ivec27_1, t_int29_1);
  const __m128i t_ivec32_1 = _mm_xor_si128(t_ivec26_1, t_ivec31_1);
  const int t_int33_1 = 51130563U;
  const int t_int34_1 = t_int33_1;
  const __m128i t_ivec35_1 = _mm_set1_epi32(t_int34_1);
  const __m128i t_ivec36_1 = (__m128i)(t_ivec35_1);
  const __m128i t_ivec37_1 = _mm_and_si128(t_ivec32_1, t_ivec36_1);
  n_1 = t_ivec37_1;

  const __m128i t_ivec38_1 = n_1;
  const __m128i t_ivec39_1 = n_1;
  const int t_int40_1 = 2U;
  const int t_int41_1 = t_int40_1;
  const __m128i t_ivec43_1 = _mm_slli_epi32((__m128i)t_ivec39_1, t_int41_1);
  const __m128i t_ivec44_1 = _mm_xor_si128(t_ivec38_1, t_ivec43_1);
  const int t_int45_1 = 153391689U;
  const int t_int46_1 = t_int45_1;
  const __m128i t_ivec47_1 = _mm_set1_epi32(t_int46_1);
  const __m128i t_ivec48_1 = (__m128i)(t_ivec47_1);
  const __m128i t_ivec49_1 = _mm_and_si128(t_ivec44_1, t_ivec48_1);
  n_1 = t_ivec49_1;

  const __m128i t_ivec50_1 = n_1;
  return t_ivec50_1;
}

inline void v_get_morton_code(unsigned int *ret, const int4 &ix, const int4 &iy,
                              const int4 &iz) {
  const __m128i x_ld_1 = ix;
  const __m128i y_ld_1 = iy;
  const __m128i z_ld_1 = iz;
  const __m128i t_ivec52_1 = x_ld_1;
  const __m128i t_ivec51_1 = (__m128i)vpartby2((int *)&(t_ivec52_1));
  const int t_int54_1 = 2U;
  const int t_int53_1 = t_int54_1;
  const __m128i t_ivec56_1 = _mm_slli_epi32((__m128i)t_ivec51_1, t_int53_1);
  const __m128i t_ivec57_1 = y_ld_1;
  const __m128i t_ivec58_1 = (__m128i)vpartby2((int *)&(t_ivec57_1));
  const int t_int59_1 = 1U;
  const int t_int60_1 = t_int59_1;
  const __m128i t_ivec62_1 = _mm_slli_epi32((__m128i)t_ivec58_1, t_int60_1);
  const __m128i t_ivec63_1 = _mm_or_si128(t_ivec56_1, t_ivec62_1);
  const __m128i t_ivec64_1 = z_ld_1;
  const __m128i t_ivec65_1 = (__m128i)vpartby2((int *)&(t_ivec64_1));
  const __m128i t_ivec66_1 = _mm_or_si128(t_ivec63_1, t_ivec65_1);
  __m128i result_1 = t_ivec66_1;

  const __m128i t_ivec67_1 = result_1;
  (*((__m128i *)(ret + 0))) = t_ivec67_1;

  return;
}

static inline void get_morton_code4(unsigned int *dst, const float4 &px,
                                    const float4 &py, const float4 &pz,
                                    const float4 &minx, const float4 &miny,
                                    const float4 &minz, const float4 &invx,
                                    const float4 &invy, const float4 &invz) {
  const float4 dx = vsub_f4(px, minx);
  const float4 dy = vsub_f4(py, miny);
  const float4 dz = vsub_f4(pz, minz);
  const float4 mdx = vmul_f4(dx, invx);
  const float4 mdy = vmul_f4(dy, invy);
  const float4 mdz = vmul_f4(dz, invz);
  int4 ix = (int4)_mm_cvttps_epi32(mdx);
  int4 iy = (int4)_mm_cvttps_epi32(mdy);
  int4 iz = (int4)_mm_cvttps_epi32(mdz);

  v_get_morton_code(dst, ix, iy, iz);
}
}
#endif

void lsgl::render::CalculateMortonCodes30(uint32_t *codes, const float *points,
                                          const real3 &bmin, const real3 &bmax,
                                          int64_t startIdx, int64_t endIdx) {
  assert(startIdx <= endIdx);

  int kDIV = (1 << 10);
  real invx = kDIV / (bmax[0] - bmin[0]);
  real invy = kDIV / (bmax[1] - bmin[1]);
  real invz = kDIV / (bmax[2] - bmin[2]);

  // printf("inv: %f, %f, %f\n", invx, invy, invz);

  int64_t n = endIdx - startIdx;

#ifdef _OPENMP
//#pragma omp parallel for if (n > 4096)
//#pragma omp parallel for schedule(dynamic, 1) if (n > 4096)
#pragma omp parallel for if (n > 4096)
#endif
  for (int64_t i = startIdx; i < endIdx; i++) {
    real3 p_i(points[3 * i + 0], points[3 * i + 1], points[3 * i + 2]);
    codes[i] = MortionCode30(p_i, bmin, invx, invy, invz);
  }
}

void
lsgl::render::CalculateMortonCodes30SIMD(uint32_t *codes, const float *points,
                                         const real3 &bmin, const real3 &bmax,
                                         int64_t startIdx, int64_t endIdx) {
#if defined(__SSE2__) || (_M_IX86_FP >= 2)
  assert(startIdx <= endIdx);

  int kDIV = (1 << 10);
  real invx = kDIV / (bmax[0] - bmin[0]);
  real invy = kDIV / (bmax[1] - bmin[1]);
  real invz = kDIV / (bmax[2] - bmin[2]);

  // printf("inv: %f, %f, %f\n", invx, invy, invz);

  int64_t n = endIdx - startIdx;
  int64_t n4 = n / 4;

  float4 vminx = vset1_f4(bmin[0]);
  float4 vminy = vset1_f4(bmin[1]);
  float4 vminz = vset1_f4(bmin[2]);

  float4 vinvx = vset1_f4(invx);
  float4 vinvy = vset1_f4(invy);
  float4 vinvz = vset1_f4(invz);

#ifdef _OPENMP
#pragma omp parallel for if (n > 1024 * 128)
#endif
  for (int64_t i = startIdx; i < startIdx + 4 * n4; i += 4) {

    // @todo { optimize load. }
    float4 vpx = vsetr4_f4(points[3 * (i + 0) + 0], points[3 * (i + 1) + 0],
                           points[3 * (i + 2) + 0], points[3 * (i + 3) + 0]);
    float4 vpy = vsetr4_f4(points[3 * (i + 0) + 1], points[3 * (i + 1) + 1],
                           points[3 * (i + 2) + 1], points[3 * (i + 3) + 1]);
    float4 vpz = vsetr4_f4(points[3 * (i + 0) + 2], points[3 * (i + 1) + 2],
                           points[3 * (i + 2) + 2], points[3 * (i + 3) + 2]);

    get_morton_code4(&codes[i], vpx, vpy, vpz, vminx, vminy, vminz, vinvx,
                     vinvy, vinvz);
  }

  // Remainder
  for (int64_t i = startIdx + 4 * n4; i < endIdx; i++) {
    real3 p_i(points[3 * i + 0], points[3 * i + 1], points[3 * i + 2]);
    codes[i] = MortionCode30(p_i, bmin, invx, invy, invz);
  }

#else
  // Fallback
  CalculateMortonCodes30(codes, points, bmin, bmax, startIdx, endIdx);
#endif
}

void lsgl::render::CalculateMortonCodesTriangleFloat30(
    uint32_t *codes, const float *points, const uint32_t *faces,
    const real3 &bmin, const real3 &bmax, int64_t startIdx, int64_t endIdx) {
  assert(startIdx <= endIdx);

  int kDIV = (1 << 10);
  real invx = kDIV / (bmax[0] - bmin[0]);
  real invy = kDIV / (bmax[1] - bmin[1]);
  real invz = kDIV / (bmax[2] - bmin[2]);

  // printf("inv: %f, %f, %f\n", invx, invy, invz);

  int64_t n = endIdx - startIdx;

  float one_third = 1.0f / 3.0f;

#ifdef _OPENMP
//#pragma omp parallel for if (n > 4096)
//#pragma omp parallel for schedule(dynamic, 1) if (n > 4096)
#pragma omp parallel for if (n > 4096)
#endif
  for (int64_t i = startIdx; i < endIdx; i++) {
    uint32_t f0 = faces[3 * i + 0];
    uint32_t f1 = faces[3 * i + 1];
    uint32_t f2 = faces[3 * i + 2];
    real3 p0(points[3 * f0 + 0], points[3 * f0 + 1], points[3 * f0 + 2]);
    real3 p1(points[3 * f1 + 0], points[3 * f1 + 1], points[3 * f1 + 2]);
    real3 p2(points[3 * f2 + 0], points[3 * f2 + 1], points[3 * f2 + 2]);
    real3 p_i;

    p_i[0] = one_third * (p0[0] + p1[0] + p2[0]);
    p_i[1] = one_third * (p0[1] + p1[1] + p2[1]);
    p_i[2] = one_third * (p0[2] + p1[2] + p2[2]);
    codes[i] = MortionCode30(p_i, bmin, invx, invy, invz);
  }
}

void lsgl::render::CalculateMortonCodesTetraFloat30(
    uint32_t *codes, const float *points, const uint32_t *faces,
    const real3 &bmin, const real3 &bmax, int64_t startIdx, int64_t endIdx) {
  assert(startIdx <= endIdx);

  int kDIV = (1 << 10);
  real invx = kDIV / (bmax[0] - bmin[0]);
  real invy = kDIV / (bmax[1] - bmin[1]);
  real invz = kDIV / (bmax[2] - bmin[2]);

  // printf("inv: %f, %f, %f\n", invx, invy, invz);

  int64_t n = endIdx - startIdx;

  float one_fourth = 1.0f / 4.0f;

#ifdef _OPENMP
//#pragma omp parallel for if (n > 4096)
//#pragma omp parallel for schedule(dynamic, 1) if (n > 4096)
#pragma omp parallel for if (n > 4096)
#endif
  for (int64_t i = startIdx; i < endIdx; i++) {
    uint32_t f0 = faces[4 * i + 0];
    uint32_t f1 = faces[4 * i + 1];
    uint32_t f2 = faces[4 * i + 2];
    uint32_t f3 = faces[4 * i + 3];
    real3 p0(points[3 * f0 + 0], points[3 * f0 + 1], points[3 * f0 + 2]);
    real3 p1(points[3 * f1 + 0], points[3 * f1 + 1], points[3 * f1 + 2]);
    real3 p2(points[3 * f2 + 0], points[3 * f2 + 1], points[3 * f2 + 2]);
    real3 p3(points[3 * f3 + 0], points[3 * f3 + 1], points[3 * f3 + 2]);
    real3 p_i;

    p_i[0] = one_fourth * (p0[0] + p1[0] + p2[0] + p3[0]);
    p_i[1] = one_fourth * (p0[1] + p1[1] + p2[1] + p3[1]);
    p_i[2] = one_fourth * (p0[2] + p1[2] + p2[2] + p3[2]);
    codes[i] = MortionCode30(p_i, bmin, invx, invy, invz);
  }
}

void lsgl::render::CalculateMortonCodesTetraDouble30(
    uint32_t *codes, const double *points, const uint32_t *faces,
    const real3 &bmin, const real3 &bmax, int64_t startIdx, int64_t endIdx) {
  assert(startIdx <= endIdx);

  int kDIV = (1 << 10);
  real invx = kDIV / (bmax[0] - bmin[0]);
  real invy = kDIV / (bmax[1] - bmin[1]);
  real invz = kDIV / (bmax[2] - bmin[2]);

  // printf("inv: %f, %f, %f\n", invx, invy, invz);

  int64_t n = endIdx - startIdx;

  float one_fourth = 1.0f / 4.0f;

#ifdef _OPENMP
//#pragma omp parallel for if (n > 4096)
//#pragma omp parallel for schedule(dynamic, 1) if (n > 4096)
#pragma omp parallel for if (n > 4096)
#endif
  for (int64_t i = startIdx; i < endIdx; i++) {
    uint32_t f0 = faces[4 * i + 0];
    uint32_t f1 = faces[4 * i + 1];
    uint32_t f2 = faces[4 * i + 2];
    uint32_t f3 = faces[4 * i + 3];

    // may truncate precision.
    real3 p0(points[3 * f0 + 0], points[3 * f0 + 1], points[3 * f0 + 2]);
    real3 p1(points[3 * f1 + 0], points[3 * f1 + 1], points[3 * f1 + 2]);
    real3 p2(points[3 * f2 + 0], points[3 * f2 + 1], points[3 * f2 + 2]);
    real3 p3(points[3 * f3 + 0], points[3 * f3 + 1], points[3 * f3 + 2]);
    real3 p_i;

    p_i[0] = one_fourth * (p0[0] + p1[0] + p2[0] + p3[0]);
    p_i[1] = one_fourth * (p0[1] + p1[1] + p2[1] + p3[1]);
    p_i[2] = one_fourth * (p0[2] + p1[2] + p2[2] + p3[2]);
    codes[i] = MortionCode30(p_i, bmin, invx, invy, invz);
  }
}

void lsgl::render::CalculateMortonCodes60(uint64_t *codes, const float *points,
                                          const real3 &bmin, const real3 &bmax,
                                          int64_t startIdx, int64_t endIdx) {
  assert(startIdx <= endIdx);

  int kDIV = (1 << 10);
  real invx = kDIV / (bmax[0] - bmin[0]);
  real invy = kDIV / (bmax[1] - bmin[1]);
  real invz = kDIV / (bmax[2] - bmin[2]);

  // printf("inv: %f, %f, %f\n", invx, invy, invz);

  int64_t n = endIdx - startIdx;

#ifdef _OPENMP
#pragma omp parallel for if (n > 4096)
#endif
  for (int64_t i = 0; i < n; i++) {
    real3 p_i(points[3 * i + 0], points[3 * i + 1], points[3 * i + 2]);
    codes[i] = MortionCode60(p_i, bmin, invx, invy, invz);
  }
}

static inline int CountLeadingZeros32(uint32_t x) {
  // Note: we can use clz(), count leading zeros instruction, or
  //       convert integert to float then caculate 31 - floor(log2(x))
  //       for better performance
  //
  // If you have popc instruction, clz() can be computed as,
  //
  // function clz(x):
  //       for each y in {1, 2, 4, 8, 16}: x ← x | (x >> y)
  //       return 32 − popc(x)
  //

  if (x == 0)
    return 32;

  int n = 0;
  if ((x & 0xFFFF0000) == 0) {
    n = n + 16;
    x = x << 16;
  }
  if ((x & 0xFF000000) == 0) {
    n = n + 8;
    x = x << 8;
  }
  if ((x & 0xF0000000) == 0) {
    n = n + 4;
    x = x << 4;
  }
  if ((x & 0xC0000000) == 0) {
    n = n + 2;
    x = x << 2;
  }
  if ((x & 0x80000000) == 0) {
    n = n + 1;
  }
  return n;
}

static inline int CountLeadingZeros64(uint64_t x) {
  uint32_t x0 = static_cast<uint32_t>(x);
  uint32_t y = x >> 32;

  return ((y > 0) ? CountLeadingZeros32(y) : 32) + CountLeadingZeros32(x0);
}

// The length of logest common prefix between keys k_i and k_j
// Property: delta(i', j') >= delta(i, j) for any i', j' in range [i, j]
// Note: j could be negative.
//
// The prefix tree algorithm doesn't allow multiple key, so we assign
// k_i' = k_i |+| i and k_j' = k_j |+| j, when k_i == k_j.
// |+| is bit concatenation.
// Thus we simply fall back to compare i and j if k_i == k_j when evaluating
// delta(i, j)
static inline int Delta30(const IndexKey30 *keys, int32_t i, int32_t j,
                          int32_t n) {
  assert((i >= 0) && (i <= (n - 1)));

  // -1 if j is not in range [0, n-1]
  if (j < 0 || j >= n)
    return -1;

  uint32_t k_i = keys[i].code;
  uint32_t k_j = keys[j].code;
  int offset = 0;

  if (k_i == k_j) {
    // fallback
    k_i = i;
    k_j = j;
    offset = 32;
  }

  // Compututing logical XOR between two keys and counting the leading zero
  // bits in the resulting value.
  uint32_t x = k_i ^ k_j;
  int c = CountLeadingZeros32(x);

  return c + offset;
}

static inline int Delta60(const IndexKey60 *keys, int64_t i, int64_t j,
                          int64_t n) {
  assert((i >= 0) && (i <= (n - 1)));

  // -1 if j is not in range [0, n-1]
  if (j < 0 || j >= n)
    return -1;

  uint64_t k_i = keys[i].code;
  uint64_t k_j = keys[j].code;
  int offset = 0;

  if (k_i == k_j) {
    // fallback
    k_i = i;
    k_j = j;
    offset = 64;
  }

  // Compututing logical XOR between two keys and counting the leading zero
  // bits in the resulting value.
  uint64_t x = k_i ^ k_j;
  int c = CountLeadingZeros64(x);

  return c + offset;
}

static inline int Sign(int x) {
  if (x >= 0)
    return 1;
  return -1;
}

// Construct binary radix tree by 30bit Morton code.
NodeInfo32
lsgl::render::ConstructBinaryRadixTree30(const IndexKey30 *keys,
                                         int i, // i in range [0, n-2]
                                         uint32_t n) {
  assert((i >= 0) && (i <= (n - 2)));

  // Determine direction of the randge(+1 or -1)
  int d = Sign(Delta30(keys, i, i + 1, n) - Delta30(keys, i, i - 1, n));

  // Compute upper bound for the length of the range.
  int delta_min = Delta30(keys, i, i - d, n);
  int l_max = 2;

  while (Delta30(keys, i, i + l_max * d, n) > delta_min) {
    assert(i != (i + l_max * d));
    l_max = l_max * 2;
  }

  // Find the other end using binary search.
  int l = 0;
  int t = l_max / 2;

  for (; t >= 1; t /= 2) {
    assert(i != (i + (l + t) * d));
    if (Delta30(keys, i, i + (l + t) * d, n) > delta_min) {
      l = l + t;
    }
  }

  int j = i + l * d;

  // Find the split position using binary search.
  int delta_node = Delta30(keys, i, j, n);

  int s = 0;

  int k = 2;
  t = (int)ceil((double)l / (double)k);

  if (t > 0) {
    // for t <- {ceil(l/2), ceil(l/4), ..., 1} do
    while (1) {
      assert(t >= 1);
      // printf("ceil(%d / %d) = %d\n", l, k, (int)ceil((double)l/(double)k));
      // printf("t = %d\n", t);
      // printf("s = %d\n", s);
      // printf(" d = %d, dnode = %d\n", Delta30(keys, i, i+(s+t)*d,n),
      // delta_node);

      assert(i != (i + (s + t) * d));
      if (Delta30(keys, i, i + (s + t) * d, n) > delta_node) {
        s = s + t;
      }

      if (t == 1)
        break;

      k = k * 2;
      t = (int)ceil((double)l / (double)k);
    }
  }

  // if (i == 14618) {
  //  printf("s = %d, t = %d\n", l, t);
  //}

  // printf("s = %d, d = %d\n", s, d);
  int r = i + s * d + (std::min)(d, 0);

  NodeInfo32 Ii;

  // Output child index
  Ii.childIndex = r;

  if ((std::min)(i, j) == r) {
    Ii.leftType = NODE_TYPE_LEAF;
  } else {
    Ii.leftType = NODE_TYPE_INTERMEDIATE;
  }

  if ((std::max)(i, j) == (r + 1)) {
    Ii.rightType = NODE_TYPE_LEAF;
  } else {
    Ii.rightType = NODE_TYPE_INTERMEDIATE;
  }

  return Ii;
}

NodeInfo64
lsgl::render::ConstructBinaryRadixTree60(const IndexKey60 *keys,
                                         int i, // i in range [0, n-2]
                                         uint64_t n) {
  assert((i >= 0) && (i <= (n - 2)));

  // Determine direction of the randge(+1 or -1)
  int d = Sign(Delta60(keys, i, i + 1, n) - Delta60(keys, i, i - 1, n));

  // Compute upper bound for the length of the range.
  int delta_min = Delta60(keys, i, i - d, n);
  int l_max = 2;

  while (Delta60(keys, i, i + l_max * d, n) > delta_min) {
    assert(i != (i + l_max * d));
    l_max = l_max * 2;
  }

  // Find the other end using binary search.
  int l = 0;
  int t = l_max / 2;

  for (; t >= 1; t /= 2) {
    assert(i != (i + (l + t) * d));
    if (Delta60(keys, i, i + (l + t) * d, n) > delta_min) {
      l = l + t;
    }
  }

  int j = i + l * d;

  // Find the split position using binary search.
  int delta_node = Delta60(keys, i, j, n);

  int s = 0;

  int k = 2;
  t = (int)ceil((double)l / (double)k);

  if (t > 0) {
    // for t <- {ceil(l/2), ceil(l/4), ..., 1} do
    while (1) {
      assert(t >= 1);
      // printf("ceil(%d / %d) = %d\n", l, k, (int)ceil((double)l/(double)k));
      // printf("t = %d\n", t);
      // printf("s = %d\n", s);
      // printf(" d = %d, dnode = %d\n", Delta30(keys, i, i+(s+t)*d,n),
      // delta_node);

      assert(i != (i + (s + t) * d));
      if (Delta60(keys, i, i + (s + t) * d, n) > delta_node) {
        s = s + t;
      }

      if (t == 1)
        break;

      k = k * 2;
      t = (int)ceil((double)l / (double)k);
    }
  }

  // if (i == 14618) {
  //  printf("s = %d, t = %d\n", l, t);
  //}

  // printf("s = %d, d = %d\n", s, d);
  int r = i + s * d + (std::min)(d, 0);

  NodeInfo64 Ii;

  // Output child index
  Ii.index = r;

  if ((std::min)(i, j) == r) {
    Ii.leftType = NODE_TYPE_LEAF;
  } else {
    Ii.leftType = NODE_TYPE_INTERMEDIATE;
  }

  if ((std::max)(i, j) == (r + 1)) {
    Ii.rightType = NODE_TYPE_LEAF;
  } else {
    Ii.rightType = NODE_TYPE_INTERMEDIATE;
  }

  return Ii;
}
