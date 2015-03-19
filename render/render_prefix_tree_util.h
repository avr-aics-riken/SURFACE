/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

//
// References
// - Maximizing Parallelism in the Construction of BVHs, Octrees, and k-d Trees
//   Tero Karras, HPG 2012
// - LUT-based morton code calculation
//   https://github.com/aavenel/FastMortonKeys
// - Produce interleaving bit patterns (morton keys) for 32 bit , 64 bit and
// 128bit
//   http://stackoverflow.com/questions/18529057/produce-interleaving-bit-patterns-morton-keys-for-32-bit-64-bit-and-128bit
#ifndef __PREFIX_TREE_UTIL_H__
#define __PREFIX_TREE_UTIL_H__

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

namespace lsgl {
namespace render {

typedef enum {
  NODE_TYPE_INTERMEDIATE = 0,
  NODE_TYPE_LEAF = 1,
} NodeType;

// @todo { save memory }
struct NodeInfo32 {
  uint32_t childIndex; // left = childIndex, right = childIndex + 1
  NodeType leftType;
  NodeType rightType;
};

struct NodeInfo64 {
  uint64_t index; // left = index, right = index + 1
  NodeType leftType;
  NodeType rightType;
};

// 30bit = up to 2G primitives.
// 60bit = For large scale dataset.

struct IndexKey30 {
  uint32_t index; // Primitive index
  uint32_t code;  // 30bit Morton code
};

struct IndexKey60 {
  uint64_t index; // Primitive index
  uint64_t code;  // 60bit Morton code
};

// Comparator for 30bit uint Radix-sort
class RadixComparatorUInt30 {
  const int bit; // bit position [0..29] to examine
public:
  RadixComparatorUInt30(int offset) : bit(offset) {}

  bool operator()(IndexKey30 value) const // functor
  {
    // Assume bit is less than 30
    return !(value.code & (1 << bit));
  }
};

// Comparator for 60bit uint Radix-sort
class RadixComparatorUInt60 {
  const int bit; // bit position [0..59] to examine
public:
  RadixComparatorUInt60(int offset) : bit(offset) {}

  bool operator()(IndexKey60 value) const // functor
  {
    // Assume bit is less than 60
    return !(value.code & (1ULL << bit));
  }
};

//
// --
//

// Morton code utilities
static inline uint32_t PartBy2_30(uint32_t n) {
  n &= 0x000003ff;
  n = (n ^ (n << 16)) & 0xff0000ff;
  n = (n ^ (n << 8)) & 0x0300f00f;
  n = (n ^ (n << 4)) & 0x030c30c3;
  n = (n ^ (n << 2)) & 0x09249249;
  return n;
}

static inline uint64_t PartBy2_60(uint32_t n) {
  uint64_t x = n & 0x1fffff; // 21 bits
  x = (x | x << 32) & 0x1f00000000ffff;
  x = (x | x << 16) & 0x1f0000ff0000ff;
  x = (x | x << 8) & 0x100f00f00f00f00f;
  x = (x | x << 4) & 0x10c30c30c30c30c3;
  x = (x | x << 2) & 0x1249249249249249;
  return x;
}

static inline uint32_t ConstructMortionBits30(uint32_t x, uint32_t y,
                                              uint32_t z) {
  return (PartBy2_30(x) << 2) | (PartBy2_30(y) << 1) | (PartBy2_30(z));
}

static inline uint64_t ConstructMortionBits60(uint32_t x, uint32_t y,
                                              uint32_t z) {
  return (PartBy2_60(x) << 2) | (PartBy2_60(y) << 1) | (PartBy2_60(z));
}

// Get 30bit mortion code.
static inline uint32_t MortionCode30(const real3 &p, const real3 &bmin,
                                     real invx, real invy, real invz) {

  uint32_t ix = (uint32_t)((p[0] - bmin[0]) * invx);
  uint32_t iy = (uint32_t)((p[1] - bmin[1]) * invy);
  uint32_t iz = (uint32_t)((p[2] - bmin[2]) * invz);

  uint32_t code = ConstructMortionBits30(ix, iy, iz);
  // printf("p = %f, %f, %f\n", p[0], p[1], p[2]);
  // printf("bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  // printf("inv = %f, %f, %f\n", invx, invy, invz);
  // printf("ijk = %d, %d, %d, code = %d\n", ix, iy, iz, code);

  return code;
}

// Get 60bit mortion code.
static inline uint64_t MortionCode60(const real3 &p, const real3 &bmin,
                                     real invx, real invy, real invz) {

  uint32_t ix = (uint32_t)((p[0] - bmin[0]) * invx);
  uint32_t iy = (uint32_t)((p[1] - bmin[1]) * invy);
  uint32_t iz = (uint32_t)((p[2] - bmin[2]) * invz);

  uint64_t code = ConstructMortionBits60(ix, iy, iz);
  // printf("p = %f, %f, %f\n", p[0], p[1], p[2]);
  // printf("bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  // printf("inv = %f, %f, %f\n", invx, invy, invz);
  // printf("ijk = %d, %d, %d, code = %d\n", ix, iy, iz, code);

  return code;
}

static std::string BitString32(int x) {
  char buf[32 + 1];
  buf[32] = '\0';

  for (int i = 0; i < 32; i++) {
    buf[i] = '0';
  }

  for (int i = 0; i < 32; i++) {
    if (i != 0 && x == 0)
      break;
    if (x % 2 == 0)
      buf[(32 - 1) - i] = '0';
    else
      buf[(32 - 1) - i] = '1';
    x = x / 2;
  }

  return std::string(buf);
}

static std::string BitString64(uint64_t x) {
  char buf[64 + 1];
  buf[64] = '\0';

  for (int i = 0; i < 64; i++) {
    buf[i] = '0';
  }

  for (int i = 0; i < 64; i++) {
    if (i != 0 && x == 0)
      break;
    if (x % 2 == 0)
      buf[(64 - 1) - i] = '0';
    else
      buf[(64 - 1) - i] = '1';
    x = x / 2;
  }

  return std::string(buf);
}

static void QuickSortKey30(IndexKey30 *v, int firstIdx, int lastIdx) {
  int startIdx[2], endIdx[2];
  IndexKey30 pivot, temp;

  if (firstIdx < lastIdx) {
    startIdx[1] = firstIdx;
    endIdx[0] = lastIdx;
    pivot = v[(firstIdx + lastIdx) / 2];
    while (startIdx[1] <= endIdx[0]) {
      while (v[startIdx[1]].code < pivot.code)
        startIdx[1]++;
      while (pivot.code < v[endIdx[0]].code)
        endIdx[0]--;
      if (startIdx[1] <= endIdx[0]) {
        temp = v[startIdx[1]];
        v[startIdx[1]] = v[endIdx[0]];
        v[endIdx[0]] = temp;
        startIdx[1]++;
        endIdx[0]--;
      }
    }
    startIdx[0] = firstIdx;
    endIdx[1] = lastIdx;

    {
      for (size_t i = 0; i <= 1; i++) {
        QuickSortKey30(v, startIdx[i], endIdx[i]);
      }
    }
  }
}

// Simple in-place OpenMP parallel quick sorter
static void QuickSortKey30OMP(IndexKey30 *v, int firstIdx, int lastIdx) {
  int startIdx[2], endIdx[2];
  IndexKey30 pivot, temp;

  if (firstIdx < lastIdx) {
    startIdx[1] = firstIdx;
    endIdx[0] = lastIdx;
    pivot = v[(firstIdx + lastIdx) / 2];
    while (startIdx[1] <= endIdx[0]) {
      while (v[startIdx[1]].code < pivot.code)
        startIdx[1]++;
      while (pivot.code < v[endIdx[0]].code)
        endIdx[0]--;
      if (startIdx[1] <= endIdx[0]) {
        temp = v[startIdx[1]];
        v[startIdx[1]] = v[endIdx[0]];
        v[endIdx[0]] = temp;
        startIdx[1]++;
        endIdx[0]--;
      }
    }
    startIdx[0] = firstIdx;
    endIdx[1] = lastIdx;

#pragma omp parallel
    {
#pragma omp for nowait
      for (size_t i = 0; i <= 1; i++) {
        QuickSortKey30(v, startIdx[i], endIdx[i]);
      }
    }
  }
}

#define MAX_RADIX_SORT_THREADS (8)

// Assume less than 2G items
// Reference:
// - Introduction to GPU Radix Sort
//   http://www.heterogeneouscompute.org/wordpress/wp-content/uploads/2011/06/RadixSort.pdf
#ifdef _OPENMP
static void RadixSort30OMP(IndexKey30 *begin, IndexKey30 *end) {

  unsigned int n = end - begin;

  if (n < 1024) {
    QuickSortKey30(begin, 0, n - 1);
    return;
  }

  IndexKey30 *begin1 = new IndexKey30[end - begin];
  IndexKey30 *end1 = begin1 + (end - begin);

  // Process 8bits(256 counts) each.
  for (unsigned shift = 0; shift < 32; shift += 8) {
    // unsigned int  count[0x100] = {0};

    unsigned int local_count[MAX_RADIX_SORT_THREADS][0x100];
    unsigned int local_offset[MAX_RADIX_SORT_THREADS][0x100];

    IndexKey30 *keys = begin;

// 1. Count
#pragma omp parallel num_threads(MAX_RADIX_SORT_THREADS)
    {

      int tid = omp_get_thread_num();
      memset(local_count[tid], 0, sizeof(unsigned int) * 0x100);

      unsigned int startIdx = (tid * n) / (MAX_RADIX_SORT_THREADS);
      unsigned int endIdx =
          std::min(n, ((tid + 1) * n) / (MAX_RADIX_SORT_THREADS));

      // printf("range (%d, %d)\n", startIdx, endIdx);

      for (size_t i = startIdx; i < endIdx; i++) {
        local_count[tid][((keys[i].code) >> shift) & 0xFF]++;
      }
    }

    // 2. Scan
    {
      for (int j = 0; j < MAX_RADIX_SORT_THREADS; j++) {
        memset(local_offset[j], 0, sizeof(unsigned int) * 0x100);
      }

      unsigned int local_sum[0x100];

      for (int i = 0; i < 0x100; i++) {
        for (int j = 1; j < MAX_RADIX_SORT_THREADS; j++) {
          // printf("local_count[%d] = %d\n", i, local_count[j-1][i]);
          local_offset[j][i] = local_offset[j - 1][i] + local_count[j - 1][i];
        }

        local_sum[i] = local_offset[MAX_RADIX_SORT_THREADS - 1][i] +
                       local_count[MAX_RADIX_SORT_THREADS - 1][i];
        // printf("local_sum[%d] = %d\n", i, local_sum[i]);
      }

      for (int i = 1; i < 0x100; i++) {
        local_sum[i] += local_sum[i - 1];
      }

      // printf("local_sum = %d\n", local_sum[0xff]);
      assert(local_sum[0xff] <= n);

      for (int i = 1; i < 0x100; i++) {
        for (int j = 0; j < MAX_RADIX_SORT_THREADS; j++) {
          local_offset[j][i] += local_sum[i - 1];
        }
      }
    }

// IndexKey30 *bucket[0x100], *q = begin1;

// Store offset
// for (int i = 0; i < 0x100; q += count[i++]) {
//  bucket[i] = q;
//}

// 3. Scatter
#pragma omp parallel num_threads(MAX_RADIX_SORT_THREADS)
    {

      int tid = omp_get_thread_num();

      // Compute per-thread offset pointer
      IndexKey30 *bucket[0x100];
      for (int i = 0; i < 0x100; i++) {
        bucket[i] = begin1 + local_offset[tid][i];
      }

      unsigned int startIdx = (tid * n) / (MAX_RADIX_SORT_THREADS);
      unsigned int endIdx =
          std::min(n, ((tid + 1) * n) / (MAX_RADIX_SORT_THREADS));

      for (size_t i = startIdx; i < endIdx; i++) {
        IndexKey30 *p = begin + i;
        *bucket[((p->code) >> shift) & 0xFF]++ = *p;
      }
    }

    std::swap(begin, begin1);
    std::swap(end, end1);
  }

  delete[] begin1;
}
#endif

static void RadixSort30(IndexKey30 *begin, IndexKey30 *end) {

  unsigned int n = end - begin;

  if (n < 1024) {
    QuickSortKey30(begin, 0, n - 1);
    return;
  }

  IndexKey30 *begin1 = new IndexKey30[end - begin];
  IndexKey30 *end1 = begin1 + (end - begin);

  // Process 8bits(256 counts) each.
  for (unsigned shift = 0; shift < 32; shift += 8) {

    unsigned int count[0x100] = {0};

    IndexKey30 *keys = begin;

    // 1. Count
    {

      unsigned int startIdx = 0;
      unsigned int endIdx = end1 - begin1;

      for (size_t i = startIdx; i < endIdx; i++) {
        count[((keys[i].code) >> shift) & 0xFF]++;
      }
    }

    IndexKey30 *bucket[0x100], *q = begin1;
    for (int i = 0; i < 0x100; q += count[i++]) {
      bucket[i] = q;
    }

    // 2. Scatter
    {

      unsigned int startIdx = 0;
      unsigned int endIdx = end1 - begin1;

      for (size_t i = startIdx; i < endIdx; i++) {
        IndexKey30 *p = begin + i;
        *bucket[((p->code) >> shift) & 0xFF]++ = *p;
      }
    }

    std::swap(begin, begin1);
    std::swap(end, end1);
  }

  delete[] begin1;
}

static void Merge30(IndexKey30 *a, IndexKey30 *b, unsigned int low,
                    unsigned int pivot, unsigned int high) {
  unsigned int h, i, j, k;
  h = low;
  i = low;
  j = pivot + 1;

  while ((h <= pivot) && (j <= high)) {
    if (a[h].code <= a[j].code) {
      b[i] = a[h];
      h++;
    } else {
      b[i] = a[j];
      j++;
    }
    i++;
  }
  if (h > pivot) {
    for (k = j; k <= high; k++) {
      b[i] = a[k];
      i++;
    }
  } else {
    for (k = h; k <= pivot; k++) {
      b[i] = a[k];
      i++;
    }
  }
  for (k = low; k <= high; k++)
    a[k] = b[k];
}

// Sort keys in range [low, high]
static void MergeSort30(IndexKey30 *a, IndexKey30 *b, unsigned int low,
                        unsigned int high) {
  unsigned int pivot;

  if (low < high) {
    pivot = (low + high) / 2;

    if (high - low < 1024 * 10) {
      // MergeSort30(a, b, low, pivot);
      // MergeSort30(a, b, pivot + 1, high);
      QuickSortKey30(a, low, high);
    } else {

#pragma omp task
      MergeSort30(a, b, low, pivot);

#pragma omp task
      MergeSort30(a, b, pivot + 1, high);

#pragma omp taskwait
    }

    Merge30(a, b, low, pivot, high);
  }
}

// Sort Morton code with radix-sort by LSB
static void RadixSortByMortionCode30LSB(IndexKey30 *firstIdx,
                                        IndexKey30 *lastIdx) {
#if 1
  for (int lsb = 0; lsb < 30; ++lsb) {
    std::stable_partition(firstIdx, lastIdx, RadixComparatorUInt30(lsb));
  }
#else
  QuickSortKey30(firstIdx, 0, lastIdx - firstIdx);
#endif
}

#ifdef _OPENMP
static void RadixSortByMortionCode30MSBTask(IndexKey30 *firstIdx,
                                            IndexKey30 *lastIdx, int msb = 29) {
  if ((firstIdx != lastIdx) && (msb >= 0)) {

    // Sequential
    IndexKey30 *midIdx =
        std::partition(firstIdx, lastIdx, RadixComparatorUInt30(msb));
    msb--; // decrement most-significant-bit

#pragma omp task firstprivate(midIdx) if (msb > 20) // Do not spawn many
    // tasks for better performance.
    RadixSortByMortionCode30MSBTask(firstIdx, midIdx,
                                    msb); // sort left partition

#pragma omp task firstprivate(midIdx) if (msb > 20) // Do not spawn many
    // tasks for better performance.
    RadixSortByMortionCode30MSBTask(midIdx, lastIdx,
                                    msb); // sort right partition

#pragma omp taskwait
  }
}
#endif

// Sort Morton code with radix-sort by MSB(recursive)
static void RadixSortByMortionCode30MSB(IndexKey30 *firstIdx,
                                        IndexKey30 *lastIdx, int msb = 29) {

#ifdef _OPENMP

#pragma omp parallel
  {
#pragma omp single
    { RadixSortByMortionCode30MSBTask(firstIdx, lastIdx, msb); }
  }

#else
  if ((firstIdx != lastIdx) && (msb >= 0)) {
    IndexKey30 *midIdx =
        std::partition(firstIdx, lastIdx, RadixComparatorUInt30(msb));
    msb--; // decrement most-significant-bit
    RadixSortByMortionCode30MSB(firstIdx, midIdx, msb); // sort left partition
    RadixSortByMortionCode30MSB(midIdx, lastIdx, msb);  // sort right partition
  }
#endif
}

static void RadixSortByMortionCode60LSB(IndexKey60 *firstIdx,
                                        IndexKey60 *lastIdx) {
  for (int lsb = 0; lsb < 60; ++lsb) {
    std::stable_partition(firstIdx, lastIdx, RadixComparatorUInt60(lsb));
  }
}

static void RadixSortByMortionCode60MSB(IndexKey60 *firstIdx,
                                        IndexKey60 *lastIdx, int msb = 59) {
  if ((firstIdx != lastIdx) && (msb >= 0)) {
    IndexKey60 *midIdx =
        std::partition(firstIdx, lastIdx, RadixComparatorUInt60(msb));
    msb--; // decrement most-significant-bit
    RadixSortByMortionCode60MSB(firstIdx, midIdx, msb); // sort left partition
    RadixSortByMortionCode60MSB(midIdx, lastIdx, msb);  // sort right partition
  }
}

void CalculateMortonCodes30(uint32_t *codes, const float *points,
                            const real3 &bmin, const real3 &bmax,
                            int64_t startIdx, int64_t endIdx);

void CalculateMortonCodes30SIMD(uint32_t *codes, const float *points,
                                const real3 &bmin, const real3 &bmax,
                                int64_t startIdx, int64_t endIdx);

void CalculateMortonCodesTriangleFloat30(uint32_t *codes, const float *vertices,
                                         const uint32_t *faces,
                                         const real3 &bmin, const real3 &bmax,
                                         int64_t startIdx, int64_t endIdx);

void CalculateMortonCodesTetraFloat30(uint32_t *codes, const float *vertices,
                                         const uint32_t *faces,
                                         const real3 &bmin, const real3 &bmax,
                                         int64_t startIdx, int64_t endIdx);

void CalculateMortonCodesTetraDouble30(uint32_t *codes, const double *vertices,
                                       const uint32_t *faces,
                                       const real3 &bmin, const real3 &bmax,
                                       int64_t startIdx, int64_t endIdx);

void CalculateMortonCodes60(uint64_t *codes, const float *points,
                            const real3 &bmin, const real3 &bmax,
                            int64_t startIdx, int64_t endIdx);

// Construct binary radix tree by 30bit Morton code.
NodeInfo32 ConstructBinaryRadixTree30(const IndexKey30 *keys,
                                      int i, // i in range [0, n-2]
                                      uint32_t n);

NodeInfo64 ConstructBinaryRadixTree60(const IndexKey60 *keys,
                                      int i, // i in range [0, n-2]
                                      uint64_t n);

} // namespace
} // namespace

#endif // __PREFIX_TREE_UTIL_H__
