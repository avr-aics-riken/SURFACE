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
//

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

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10
#include <emmintrin.h>
#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2
#include <emmintrin.h>
#endif

#include "render_accel_particle.h"
#include "render_timerutil.h"
#include "render_prefix_tree_util.h"

using namespace lsgl::render;

#define MAX_LEAF_ELEMENTS (16)
#define MAX_TREE_DEPTH_32BIT                                                   \
  (22) // FYI, log2(1G/16) ~= 25.897, log2(1G/32) ~= 21

#define USE_BLOB (0)
#define ENABLE_SIMD_ISECTOR (1)

#define ENABLE_TRACE_PRINT (0)
#define ENABLE_DEBUG_PRINT (0)

#define trace(f, ...)                                                          \
  {                                                                            \
    if (ENABLE_TRACE_PRINT)                                                    \
      printf(f, __VA_ARGS__);                                                  \
  }

#if ENABLE_DEBUG_PRINT
#define debug(message, ...)                                                    \
  fprintf(stderr, "[LSGL] %s(%s:%d) " message "\n", __FUNCTION__, __FILE__,    \
          __LINE__, ##__VA_ARGS__)
#else
#define debug(message, ...) (void(0))
#endif

namespace {

#if ENABLE_SIMD_ISECTOR

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10

typedef __m128d double2;

#define FORCEINLINE __attribute__((always_inline))

#define vseld2(a, b, cond) _fjsp_selmov_v2r8((a), (b), (cond))
#define vabsd2(a) _fjsp_abs_v2r8((a))
#define vzerod2() _mm_setzero_pd()
#define voned2() _mm_set_pd(1.0, 1.0)

static FORCEINLINE double2 vmaxd2(double2 x, double2 y) {
  return _fjsp_max_v2r8(x, y);
}

static FORCEINLINE double2 vmind2(double2 x, double2 y) {
  return _fjsp_min_v2r8(x, y);
}

static FORCEINLINE double2 vinvd2(double2 x) {
  // Reciprocal + Newton Raphson improvement
  // (2 * Rcp(x)) - (x * Rcp(x) * Rcp(x))]
  // = (2 - (x * Rcp(x)) * Rcp(x)
  const __m128d two = _mm_set_pd(2.0, 2.0);
  const __m128d rcp = _fjsp_rcpa_v2r8(x);
  const __m128d rcp_nr0 = _fjsp_nmsub_v2r8(x, rcp, two);
  const __m128d rcp_nr1 = _mm_mul_pd(rcp_nr0, rcp);

  return rcp_nr1;
}

static FORCEINLINE double2 vsqrtd2(double2 arg) {
  const __m128d zero = _mm_setzero_pd();
  const __m128d three = _mm_set_pd(3.0, 3.0);
  const __m128d half = _mm_set_pd(0.5, 0.5);
  const __m128d approx = _fjsp_rsqrta_v2r8(arg);
  const __m128d cond = _mm_cmpeq_pd(arg, zero);
  const __m128d aa = _mm_mul_pd(approx, arg);
  const __m128d aaa3 = _fjsp_nmsub_v2r8(approx, aa, three);
  const __m128d aaaa3 = _mm_mul_pd(aaa3, approx);
  const __m128d nr_rsqrt = _mm_mul_pd(aaaa3, half);
  const __m128d sqr = _mm_mul_pd(arg, nr_rsqrt);

  return _fjsp_selmov_v2r8(zero, sqr, cond);
}

inline double2 vdotd2(double2 ax, double2 ay, double2 az, double2 bx,
                      double2 by, double2 bz) {
  return _mm_add_pd(_mm_mul_pd(ax, bx),
                    _mm_add_pd(_mm_mul_pd(ay, by), _mm_mul_pd(az, bz)));
}
#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

typedef __m128d double2;

#ifndef FORCEINLINE
#define FORCEINLINE __attribute__((always_inline))
#endif

#define vseld2(a, b, m) _mm_or_pd(_mm_and_pd((m), (a)), _mm_andnot_pd((m), (b)))
#define vabsd2(a)                                                              \
  _mm_andnot_pd(_mm_castsi128_pd(_mm_set1_epi64x(0x8000000000000000ULL)), a)
#define vzerod2() _mm_setzero_pd()
#define voned2() _mm_set_pd(1.0, 1.0)

static FORCEINLINE double2 vmaxd2(double2 x, double2 y) {
  return _mm_max_pd(x, y);
}

static FORCEINLINE double2 vmind2(double2 x, double2 y) {
  return _mm_min_pd(x, y);
}

static FORCEINLINE double2 vinvd2(double2 x) {
  __m128d b = _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x)));
  b = _mm_sub_pd(_mm_add_pd(b, b), _mm_mul_pd(_mm_mul_pd(x, b), b));
  b = _mm_sub_pd(_mm_add_pd(b, b), _mm_mul_pd(_mm_mul_pd(x, b), b));
  return b;
}

static FORCEINLINE double2 vsqrtd2(double2 a) {
  double2 v0_5 = _mm_set1_pd(0.5);
  double2 v3_0 = _mm_set1_pd(3.0);
  double2 b = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(a)));
  b = _mm_mul_pd(_mm_mul_pd(v0_5, b),
                 _mm_sub_pd(v3_0, _mm_mul_pd(_mm_mul_pd(a, b), b)));
  b = _mm_mul_pd(_mm_mul_pd(v0_5, b),
                 _mm_sub_pd(v3_0, _mm_mul_pd(_mm_mul_pd(a, b), b)));
  return _mm_mul_pd(a, b); // sqrt(a) = (1/sqrt(a))*a
}

inline void vcrossd2(double2 &cx, double2 &cy, double2 &cz, double2 ax,
                     double2 ay, double2 az, double2 bx, double2 by,
                     double2 bz) {
  cx = _mm_sub_pd(_mm_mul_pd(ay, bz), _mm_mul_pd(az, by));
  cy = _mm_sub_pd(_mm_mul_pd(az, bx), _mm_mul_pd(ax, bz));
  cz = _mm_sub_pd(_mm_mul_pd(ax, by), _mm_mul_pd(ay, bx));
}

inline double2 vdotd2(double2 ax, double2 ay, double2 az, double2 bx,
                      double2 by, double2 bz) {
  return _mm_add_pd(_mm_mul_pd(ax, bx),
                    _mm_add_pd(_mm_mul_pd(ay, by), _mm_mul_pd(az, bz)));
}

inline void
IsectD2(double2 &tInOut, double2 &tidInOut, // Use as integer(53bit)
        double2 &uOut, double2 &vOut, const double2 &v0x, const double2 &v0y,
        const double2 &v0z, const double2 &v1x, const double2 &v1y,
        const double2 &v1z, const double2 &v2x, const double2 &v2y,
        const double2 &v2z, const double2 &rox, const double2 &roy,
        const double2 &roz, const double2 &rdx, const double2 &rdy,
        const double2 &rdz,
        const double2 &tid,             // Triangle ID as interger(53bit)
        const double2 &doubleSidedMask) // Assume all bits are 1 or 0.
{
  const double kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  // e1 = v1 - v0;
  const double2 e1x = _mm_sub_pd(v1x, v0x);
  const double2 e1y = _mm_sub_pd(v1y, v0y);
  const double2 e1z = _mm_sub_pd(v1z, v0z);

  // e2 = v2 - v0;
  const double2 e2x = _mm_sub_pd(v2x, v0x);
  const double2 e2y = _mm_sub_pd(v2y, v0y);
  const double2 e2z = _mm_sub_pd(v2z, v0z);

  double2 px, py, pz;

  // p = cross(rayDir, e2)
  vcrossd2(px, py, pz, rdx, rdy, rdz, e2x, e2y, e2z);

  // det = vdot(e1, p)
  const double2 det = vdotd2(e1x, e1y, e1z, px, py, pz);
  const double2 ddet = vseld2(vabsd2(det), det, doubleSidedMask);

  // 1.0 / det
  const double2 invDet = vinvd2(ddet);

  const double2 detMask = _mm_cmpgt_pd(ddet, _mm_set_pd(kEPS, kEPS));

  // s = rayOrg - p0;
  const double2 sx = _mm_sub_pd(rox, v0x);
  const double2 sy = _mm_sub_pd(roy, v0y);
  const double2 sz = _mm_sub_pd(roz, v0z);

  // q = vcross(s, e1);
  double2 qx, qy, qz;
  vcrossd2(qx, qy, qz, sx, sy, sz, e1x, e1y, e1z);

  // double u = vdotd(s, p) * invDet;
  // double v = vdotd(q, rayDir) * invDet;
  // double t = vdotd(e2, q) * invDet;
  double2 u = _mm_mul_pd(vdotd2(sx, sy, sz, px, py, pz), invDet);
  double2 v = _mm_mul_pd(vdotd2(qx, qy, qz, rdx, rdy, rdz), invDet);
  double2 t = _mm_mul_pd(vdotd2(e2x, e2y, e2z, qx, qy, qz), invDet);

  double2 vzero = vzerod2();
  double2 vone = voned2();

  // if (u < 0.0 || u > 1.0)
  //  return false;
  double2 uMask = _mm_and_pd(_mm_cmpge_pd(u, vzero), _mm_cmple_pd(u, vone));

  // if (v < 0.0 || u + v > 1.0)
  //  return false;
  double2 uvMask =
      _mm_and_pd(_mm_cmpge_pd(v, vzero), _mm_cmple_pd(_mm_add_pd(u, v), vone));

  // if (t < 0.0 || t > tInOut)
  //  return false;
  double2 tMask = _mm_and_pd(_mm_cmpge_pd(t, vzero), _mm_cmple_pd(t, tInOut));

  double2 finalMask =
      _mm_and_pd(_mm_and_pd(detMask, tMask), _mm_and_pd(uMask, uvMask));

  tInOut = vseld2(t, tInOut, finalMask);
  tidInOut = vseld2(tid, tidInOut, finalMask);
  uOut = vseld2(u, uOut, finalMask);
  vOut = vseld2(v, vOut, finalMask);

  // if (1) { // dbg
  //  double buf[2];
  //  _mm_store_pd(buf, tInOut);
  //  printf("msk: %f, %f\n", buf[0], buf[1]);
  //}
}
#else
#error SIMD is not supported on this architecture. Disable ENABLE_SIMD_ISECTOR to solve this error.
#endif

#endif // ENABLE_SIMD_ISECTOR

inline real GetRadius(const Particles *particles, unsigned int index) {
  if (particles->radius) {
    return particles->radius[index];
  } else {
    return particles->constantRadius;
  }
}

class PartPred : public std::unary_function<unsigned int, bool> {
public:
  PartPred(int axis, real pos, const Particles *particles)
      : axis_(axis), pos_(pos), particles_(particles) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;

    real3 p(particles_->positions[3 * i + 0], particles_->positions[3 * i + 1],
            particles_->positions[3 * i + 2]);

    real center = p[axis];

    return (center < pos_);
  }

private:
  int axis_;
  real pos_;
  const Particles *particles_;
};

#ifdef _OPENMP
void ComputeBoundingBoxOMP(real3 &bmin, real3 &bmax, const Particles *particles,
                           unsigned int *indices, unsigned int leftIndex,
                           unsigned int rightIndex) {

  // assert(leftIndex < rightIndex);
  // assert(rightIndex - leftIndex > 0);

  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  {
    size_t idx = indices[leftIndex];
    real radius = GetRadius(particles, idx);
    bmin[0] = particles->positions[3 * idx + 0] - radius - kEPS;
    bmin[1] = particles->positions[3 * idx + 1] - radius - kEPS;
    bmin[2] = particles->positions[3 * idx + 2] - radius - kEPS;
    bmax[0] = particles->positions[3 * idx + 0] + radius + kEPS;
    bmax[1] = particles->positions[3 * idx + 1] + radius + kEPS;
    bmax[2] = particles->positions[3 * idx + 2] + radius + kEPS;
  }

  // We could use min and max reduction if OpenMP 3.1 ready compiler was
  // available.

  real local_bmin[3] = {bmin[0], bmin[1], bmin[2]};
  real local_bmax[3] = {bmax[0], bmax[1], bmax[2]};

  size_t n = rightIndex - leftIndex;

#pragma omp parallel firstprivate(local_bmin, local_bmax) if (n > (1024 * 1024))
  {

#pragma omp for
    for (size_t i = leftIndex; i < rightIndex; i++) {

      size_t idx = indices[i];
      real radius = GetRadius(particles, idx);

      real minval_x = particles->positions[3 * idx + 0] - radius - kEPS;
      real minval_y = particles->positions[3 * idx + 1] - radius - kEPS;
      real minval_z = particles->positions[3 * idx + 2] - radius - kEPS;

      real maxval_x = particles->positions[3 * idx + 0] + radius + kEPS;
      real maxval_y = particles->positions[3 * idx + 1] + radius + kEPS;
      real maxval_z = particles->positions[3 * idx + 2] + radius + kEPS;

      local_bmin[0] = std::min(local_bmin[0], minval_x);
      local_bmin[1] = std::min(local_bmin[1], minval_y);
      local_bmin[2] = std::min(local_bmin[2], minval_z);

      local_bmax[0] = std::max(local_bmax[0], maxval_x);
      local_bmax[1] = std::max(local_bmax[1], maxval_y);
      local_bmax[2] = std::max(local_bmax[2], maxval_z);
    }

#pragma omp critical
    {
      for (int k = 0; k < 3; k++) {

        if (local_bmin[k] < bmin[k]) {
          {
            if (local_bmin[k] < bmin[k])
              bmin[k] = local_bmin[k];
          }
        }

        if (local_bmax[k] > bmax[k]) {
          {
            if (local_bmax[k] > bmax[k])
              bmax[k] = local_bmax[k];
          }
        }
      }
    }
  }
}
#endif

void ComputeBoundingBox(real3 &bmin, real3 &bmax, const Particles *particles,
                        unsigned int *indices, unsigned int leftIndex,
                        unsigned int rightIndex) {

  // assert(leftIndex < rightIndex);
  // assert(rightIndex - leftIndex > 0);

  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  {
    size_t idx = indices[leftIndex];
    real radius = GetRadius(particles, idx);
    bmin[0] = particles->positions[3 * idx + 0] - radius - kEPS;
    bmin[1] = particles->positions[3 * idx + 1] - radius - kEPS;
    bmin[2] = particles->positions[3 * idx + 2] - radius - kEPS;
    bmax[0] = particles->positions[3 * idx + 0] + radius + kEPS;
    bmax[1] = particles->positions[3 * idx + 1] + radius + kEPS;
    bmax[2] = particles->positions[3 * idx + 2] + radius + kEPS;
  }

  for (size_t i = leftIndex; i < rightIndex; i++) {

    size_t idx = indices[i];
    real radius = GetRadius(particles, idx);

    for (int k = 0; k < 3; k++) { // xyz
      real minval = particles->positions[3 * idx + k] - radius - kEPS;
      real maxval = particles->positions[3 * idx + k] + radius + kEPS;
      if (bmin[k] > minval)
        bmin[k] = minval;
      if (bmax[k] < maxval)
        bmax[k] = maxval;
    }
  }
}

void ComputeBoundingBox30(real3 &bmin, real3 &bmax, const Particles *particles,
                          const std::vector<IndexKey30> &keys,
                          uint32_t leftIndex, uint32_t rightIndex) {

  bmin[0] = std::numeric_limits<real>::max();
  bmin[1] = std::numeric_limits<real>::max();
  bmin[2] = std::numeric_limits<real>::max();
  bmax[0] = -std::numeric_limits<real>::max();
  bmax[1] = -std::numeric_limits<real>::max();
  bmax[2] = -std::numeric_limits<real>::max();

  if ((rightIndex - leftIndex) == 0) {
    // empty.
    return;
  }

  if (leftIndex >= rightIndex) {
    printf("left = %d, right = %d\n", leftIndex, rightIndex);
    assert(leftIndex < rightIndex);
  }

  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  size_t i = leftIndex;
  size_t idx = keys[i].index;
  real radius = GetRadius(particles, idx);
  bmin[0] = particles->positions[3 * idx + 0] - radius - kEPS;
  bmin[1] = particles->positions[3 * idx + 1] - radius - kEPS;
  bmin[2] = particles->positions[3 * idx + 2] - radius - kEPS;
  bmax[0] = particles->positions[3 * idx + 0] + radius + kEPS;
  bmax[1] = particles->positions[3 * idx + 1] + radius + kEPS;
  bmax[2] = particles->positions[3 * idx + 2] + radius + kEPS;

  for (i = leftIndex; i < rightIndex; i++) {

    idx = keys[i].index;
    radius = GetRadius(particles, idx);

    for (int k = 0; k < 3; k++) { // xyz
      real minval = particles->positions[3 * idx + k] - radius - kEPS;
      real maxval = particles->positions[3 * idx + k] + radius + kEPS;
      if (bmin[k] > minval)
        bmin[k] = minval;
      if (bmax[k] < maxval)
        bmax[k] = maxval;
    }
  }
}

inline void InvalidateBoundingBox(real3 &bmin, real3 &bmax) {
  bmin[0] = bmin[1] = bmin[2] = std::numeric_limits<real>::max();
  bmax[0] = bmax[1] = bmax[2] = -std::numeric_limits<real>::max();
}

inline void MergeBoundingBox(real3 &bmin, real3 &bmax, const real3 &leftBMin,
                             const real3 &leftBMax, const real3 &rightBMin,
                             const real3 &rightBMax) {
  bmin = leftBMin;
  bmax = leftBMax;

  for (int k = 0; k < 3; k++) {
    bmin[k] = std::min(bmin[k], rightBMin[k]);
    bmax[k] = std::max(bmax[k], rightBMax[k]);
  }
}

inline int CountLeadingZeros32(uint32_t x) {
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

int GetSplitAxis(uint32_t key) {
  int clz = CountLeadingZeros32(key);

  int n = clz - 2;
  assert(n >= 0);
  assert(n < 30);

  // 0 -> x, 1 -> y, 2 -> z, 3 -> x, ...
  return n % 3;
}

void MakeLeaf32(ParticleNode &leaf, const Particles *particles, real3 &bmin,
                real3 &bmax, const std::vector<IndexKey30> &keys,
                uint32_t leftIndex, uint32_t rightIndex) {
  // 1. Compute leaf AABB
  ComputeBoundingBox30(bmin, bmax, particles, keys, leftIndex, rightIndex);

  // 2. Create leaf node. `n' may be null(create empty leaf node for that case.)
  int n = rightIndex - leftIndex;

  leaf.bmin[0][0] = bmin[0];
  leaf.bmin[0][1] = bmin[1];
  leaf.bmin[0][2] = bmin[2];

  leaf.bmax[0][0] = bmax[0];
  leaf.bmax[0][1] = bmax[1];
  leaf.bmax[0][2] = bmax[2];

  leaf.flag = 1; // leaf
  leaf.data[0] = n;
  leaf.data[1] = (uint32_t)leftIndex;

  return;
}

//
// Build BVH tree from bottom to up manner
//
size_t BuildTreeRecursive32(std::vector<ParticleNode> &nodes, real3 &bmin,
                            real3 &bmax, const std::vector<IndexKey30> &keys,
                            const std::vector<NodeInfo32> &nodeInfos,
                            const Particles *particles, uint32_t rootIndex,
                            uint32_t leftIndex, uint32_t rightIndex,
                            bool isLeaf, int depth) {
  InvalidateBoundingBox(bmin, bmax);

  uint32_t n = rightIndex - leftIndex;
  // printf("single: [%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n",
  // depth, rootIndex, leftIndex, rightIndex, isLeaf, n);

  // if (depth > MAX_TREE_DEPTH_32BIT) {
  //  return nodes.size();
  //}

  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= MAX_LEAF_ELEMENTS) || (depth > MAX_TREE_DEPTH_32BIT)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;
    // if (leftIndex == rightIndex) { // this would be OK. 1 tri in 1 leaf case.
    //  endIndex++;
    //}

    ParticleNode leaf;
    MakeLeaf32(leaf, particles, bmin, bmax, keys, leftIndex, endIndex);

    size_t offset = nodes.size();
    nodes.push_back(leaf); // need atomic update.

    return offset;
  }

  //
  // Intermediate node. We already know split position.
  //
  uint32_t midIndex = nodeInfos[rootIndex].childIndex;

  bool isLeftLeaf =
      (nodeInfos[rootIndex].leftType == NODE_TYPE_LEAF) ? true : false;
  bool isRightLeaf =
      (nodeInfos[rootIndex].rightType == NODE_TYPE_LEAF) ? true : false;

  ParticleNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  size_t offset = nodes.size();
  nodes.push_back(node);

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  if (midIndex < leftIndex) {
    printf("rootIndex = %d, mid = %d, left = %d, leftLeaf = %d\n", rootIndex,
           midIndex, leftIndex, isLeftLeaf);
    assert(leftIndex <= midIndex);
  }

  if (midIndex > rightIndex) {
    printf("rootIndex = %d, mid = %d, right = %d, rightLeaf = %d\n", rootIndex,
           midIndex, rightIndex, isRightLeaf);
    assert(midIndex <= rightIndex);
  }

  size_t leftChildIndex = BuildTreeRecursive32(
      nodes, leftBMin, leftBMax, keys, nodeInfos, particles, midIndex,
      leftIndex, midIndex, isLeftLeaf, depth + 1);
  size_t rightChildIndex = BuildTreeRecursive32(
      nodes, rightBMin, rightBMax, keys, nodeInfos, particles, midIndex + 1,
      midIndex + 1, rightIndex, isRightLeaf, depth + 1);

  MergeBoundingBox(bmin, bmax, leftBMin, leftBMax, rightBMin, rightBMax);

  // printf("[%d] -> (%d, %d)\n", offset, leftChildIndex, rightChildIndex);

  node.data[0] = leftChildIndex;
  node.data[1] = rightChildIndex;

  node.bmin[0][0] = leftBMin[0];
  node.bmin[0][1] = leftBMin[1];
  node.bmin[0][2] = leftBMin[2];

  node.bmax[0][0] = leftBMax[0];
  node.bmax[0][1] = leftBMax[1];
  node.bmax[0][2] = leftBMax[2];

  node.bmin[1][0] = rightBMin[0];
  node.bmin[1][1] = rightBMin[1];
  node.bmin[1][2] = rightBMin[2];

  node.bmax[1][0] = rightBMax[0];
  node.bmax[1][1] = rightBMax[1];
  node.bmax[1][2] = rightBMax[2];

  nodes[offset] = node;

  return offset;
}

#if 0
void BottomUpBuildTreeRecursive(std::vector<ParticleNode> &nodes, real3 &bmin,
                            real3 &bmax, const std::vector<IndexKey30> &keys,
                            const std::vector<NodeInfo32> &nodeInfos,
                            const Particles *particles, uint32_t rootIndex,
                            uint32_t leftIndex, uint32_t rightIndex,
                            int depth) {
  InvalidateBoundingBox(bmin, bmax);

  if (rootIndex == 0) { // Top of tree.
    // nothing to do.
    return;
  }

  NodeInfo32 node = nodeInfos[rootIndex];

  if (node.isLeaf) {

  }


  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= 0)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex;
    if (leftIndex == rightIndex) { // this would be OK. 1 tri in 1 leaf case.
      endIndex++;
    }

    ParticleNode leaf;
    MakeLeaf32(leaf, particles, bmin, bmax, keys, leftIndex, endIndex);

    // @{ critical section }
    size_t offset;
#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(leaf); // need atomic update.
    }

    return offset;
  }

  //
  // Intermediate node. We already know split position.
  //
  uint32_t midIndex = nodeInfos[rootIndex].childIndex;
  bool isLeftLeaf =
      (nodeInfos[rootIndex].leftType == NODE_TYPE_LEAF) ? true : false;
  bool isRightLeaf =
      (nodeInfos[rootIndex].rightType == NODE_TYPE_LEAF) ? true : false;

  ParticleNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  // @{ critical section }
  size_t offset;
#pragma omp critical
  { 
    offset = nodes.size();
    nodes.push_back(node);
  }

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  size_t leftChildIndex = -1, rightChildIndex = -1;

  if (depth > 3) {

    // Enough number of tasks was launched. Switch to sequential code.
    
    std::vector<ParticleNode> sub_nodes;

    leftChildIndex = BuildTreeRecursive32(
        sub_nodes, leftBMin, leftBMax, keys, nodeInfos, particles, midIndex,
        leftIndex, midIndex, isLeftLeaf, depth + 1);

    rightChildIndex = BuildTreeRecursive32(
        sub_nodes, rightBMin, rightBMax, keys, nodeInfos, particles, midIndex + 1,
        midIndex + 1, rightIndex, isRightLeaf, depth + 1);

  } else {

#pragma omp task shared(leftChildIndex, nodes, leftBMin, leftBMax, keys,       \
                        nodeInfos, particles) firstprivate(midIndex)
    leftChildIndex = BuildTreeRecursive32OMP(
        nodes, leftBMin, leftBMax, keys, nodeInfos, particles, midIndex,
        leftIndex, midIndex, isLeftLeaf, depth + 1);

#pragma omp task shared(rightChildIndex, nodes, rightBMin, rightBMax, keys,    \
                        nodeInfos) firstprivate(midIndex)
    rightChildIndex = BuildTreeRecursive32OMP(
        nodes, rightBMin, rightBMax, keys, nodeInfos, particles, midIndex + 1,
        midIndex + 1, rightIndex, isRightLeaf, depth + 1);

#pragma omp taskwait
  }

  MergeBoundingBox(bmin, bmax, leftBMin, leftBMax, rightBMin, rightBMax);

  node.data[0] = leftChildIndex;
  node.data[1] = rightChildIndex;

  node.bmin[0] = bmin[0];
  node.bmin[1] = bmin[1];
  node.bmin[2] = bmin[2];

  node.bmax[0] = bmax[0];
  node.bmax[1] = bmax[1];
  node.bmax[2] = bmax[2];

  // @{ critical section }
#pragma omp critical
  {
    nodes[offset] = node;
  }

  return offset;
}
#endif

#ifdef _OPENMP
size_t BuildTreeRecursive32OMP(std::vector<ParticleNode> &nodes, real3 &bmin,
                               real3 &bmax, const std::vector<IndexKey30> &keys,
                               const std::vector<NodeInfo32> &nodeInfos,
                               const Particles *particles, uint32_t rootIndex,
                               uint32_t leftIndex, uint32_t rightIndex,
                               bool isLeaf, int depth) {
  InvalidateBoundingBox(bmin, bmax);

  uint32_t n = rightIndex - leftIndex;

  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= MAX_LEAF_ELEMENTS) || (depth > MAX_TREE_DEPTH_32BIT)) {
    // if (isLeaf || (n <= 0)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;

    ParticleNode leaf;
    MakeLeaf32(leaf, particles, bmin, bmax, keys, leftIndex, endIndex);

    // @{ critical section }
    size_t offset;
#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(leaf); // need atomic update.
    }

    return offset;
  }

  //
  // Intermediate node. We already know split position.
  //
  uint32_t midIndex = nodeInfos[rootIndex].childIndex;
  bool isLeftLeaf =
      (nodeInfos[rootIndex].leftType == NODE_TYPE_LEAF) ? true : false;
  bool isRightLeaf =
      (nodeInfos[rootIndex].rightType == NODE_TYPE_LEAF) ? true : false;

  ParticleNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  size_t offset;

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

  if (depth > 6) {

    // Enough number of tasks was launched. Switch to sequential code.
    std::vector<ParticleNode> sub_nodes;

    leftChildIndex = BuildTreeRecursive32(
        sub_nodes, leftBMin, leftBMax, keys, nodeInfos, particles, midIndex,
        leftIndex, midIndex, isLeftLeaf, depth + 1);

    rightChildIndex = BuildTreeRecursive32(
        sub_nodes, rightBMin, rightBMax, keys, nodeInfos, particles,
        midIndex + 1, midIndex + 1, rightIndex, isRightLeaf, depth + 1);

#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(node);

      // printf("offset = %d\n", offset);

      // add sub nodes
      for (size_t i = 0; i < sub_nodes.size(); i++) {
        if (sub_nodes[i].flag == 0) { // branch node
          // printf("sub[%d] = %d, %d\n", i, sub_nodes[i].data[0],
          // sub_nodes[i].data[1]);
          sub_nodes[i].data[0] += offset + 1;
          sub_nodes[i].data[1] += offset + 1;
          assert(sub_nodes[i].data[0] < nodes.size() + sub_nodes.size());
          assert(sub_nodes[i].data[1] < nodes.size() + sub_nodes.size());
          // printf("-> offt = %d, sublen = %d, sub[%d] = %d, %d\n", offset,
          // sub_nodes.size(), i, sub_nodes[i].data[0], sub_nodes[i].data[1]);
        }
      }

      nodes.insert(nodes.end(), sub_nodes.begin(), sub_nodes.end());
      leftChildIndex += offset + 1;
      rightChildIndex += offset + 1;
      assert(leftChildIndex < nodes.size() + sub_nodes.size());
      assert(rightChildIndex < nodes.size() + sub_nodes.size());
    }

  } else {

#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(node);
    }

#pragma omp task shared(leftChildIndex, rightChildIndex, nodes, leftBMin,      \
                        leftBMax, keys, nodeInfos,                             \
                        particles) firstprivate(midIndex) if (depth < 10)
    leftChildIndex = BuildTreeRecursive32OMP(
        nodes, leftBMin, leftBMax, keys, nodeInfos, particles, midIndex,
        leftIndex, midIndex, isLeftLeaf, depth + 1);

#pragma omp task shared(leftIndex, rightChildIndex, nodes, rightBMin,          \
                        rightBMax, keys,                                       \
                        nodeInfos) firstprivate(midIndex) if (depth < 10)
    rightChildIndex = BuildTreeRecursive32OMP(
        nodes, rightBMin, rightBMax, keys, nodeInfos, particles, midIndex + 1,
        midIndex + 1, rightIndex, isRightLeaf, depth + 1);

#pragma omp taskwait
  }

  assert(leftChildIndex != (size_t)(-1));
  assert(rightChildIndex != (size_t)(-1));

  MergeBoundingBox(bmin, bmax, leftBMin, leftBMax, rightBMin, rightBMax);

  // printf("[%d] -> (%d, %d)\n", offset, leftChildIndex, rightChildIndex);

  node.data[0] = leftChildIndex;
  node.data[1] = rightChildIndex;

  node.bmin[0][0] = leftBMin[0];
  node.bmin[0][1] = leftBMin[1];
  node.bmin[0][2] = leftBMin[2];

  node.bmax[0][0] = leftBMax[0];
  node.bmax[0][1] = leftBMax[1];
  node.bmax[0][2] = leftBMax[2];

  node.bmin[1][0] = rightBMin[0];
  node.bmin[1][1] = rightBMin[1];
  node.bmin[1][2] = rightBMin[2];

  node.bmax[1][0] = rightBMax[0];
  node.bmax[1][1] = rightBMax[1];
  node.bmax[1][2] = rightBMax[2];

#pragma omp critical
  { nodes[offset] = node; }

  return offset;
}
#endif

} // namespace

//
// --
//

//#ifdef _OPENMP
#if 0
size_t ParticleAccel::BuildTree(const Particles *particles,
                                unsigned int leftIdx, unsigned int rightIdx,
                                int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  size_t offset;

#pragma omp critical
  offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
#pragma omp critical
    stats_.maxTreeDepth = depth;
  }

  real3 bmin, bmax;
  ComputeBoundingBox(bmin, bmax, particles, &indices_.at(0), leftIdx, rightIdx);

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {

    // Create leaf node. `n' may be null(create empty leaf node for that case.)
    ParticleNode leaf;

    leaf.bmin[0] = bmin[0];
    leaf.bmin[1] = bmin[1];
    leaf.bmin[2] = bmin[2];

    leaf.bmax[0] = bmax[0];
    leaf.bmax[1] = bmax[1];
    leaf.bmax[2] = bmax[2];

    assert(leftIdx < std::numeric_limits<unsigned int>::max());

    leaf.flag = 1; // leaf
    leaf.data[0] = n;
    leaf.data[1] = (unsigned int)leftIdx;

#pragma omp critical
    {
      nodes_.push_back(leaf);
      stats_.numLeafNodes++;
    }

    return offset;
  }

  //
  // Create branch node.
  //

  //
  // Simple spatial median.
  //

  // Begin with longest axis
  real3 bsize;
  bsize[0] = bmax[0] - bmin[0];
  bsize[1] = bmax[1] - bmin[1];
  bsize[2] = bmax[2] - bmin[2];

  real maxsize = bsize[0];
  int longestAxis = 0;

  if (maxsize < bsize[1]) {
    maxsize = bsize[1];
    longestAxis = 1;
  }

  if (maxsize < bsize[2]) {
    maxsize = bsize[2];
    longestAxis = 2;
  }

  // Try all 3 axis until good cut position avaiable.
  unsigned int midIdx = 0;
  int cutAxis = 0;
  for (int axisTry = 0; axisTry < 3; axisTry++) {

    unsigned int *begin = &indices_[leftIdx];
    // unsigned int* end   = &indices_[rightIdx];
    unsigned int *end = begin + (rightIdx - leftIdx);
    unsigned int *mid = 0;

    // Try longest axis first.
    cutAxis = (longestAxis + axisTry) % 3;
    real cutPos = bmin[cutAxis] + 0.5 * bsize[cutAxis]; // spatial median

    //
    // Split at (cutAxis, cutPos)
    // indices_ will be modified.
    //
    mid = std::partition(begin, end, PartPred(cutAxis, cutPos, particles));

    midIdx = leftIdx + (mid - begin);
    if ((midIdx == leftIdx) || (midIdx == rightIdx)) {

      // Can't split well.
      // Switch to object median(which may create unoptimized tree, but stable)
      midIdx = leftIdx + (n >> 1);

      // Try another axis if there's axis to try.

    } else {

      // Found good cut. exit loop.
      break;
    }
  }

  ParticleNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch

#pragma omp critical
  nodes_.push_back(node);

  // Recurively split tree.
  unsigned int leftChildIndex, rightChildIndex;
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

#pragma omp task firstprivate(midIdx, depth) if (depth < 16)
  leftChildIndex = BuildTree(particles, leftBMin, leftBMax, leftIdx, midIdx, depth + 1);

#pragma omp task firstprivate(midIdx, depth) if (depth < 16)
  rightChildIndex = BuildTree(particles, rightBMin, rightBMax, midIdx, rightIdx, depth + 1);

#pragma omp taskwait
 
  node.data[0] = leftChildIndex;
  node.data[1] = rightChildIndex;

  node.bmin[0][0] = leftBMin[0];
  node.bmin[0][1] = leftBMin[1];
  node.bmin[0][2] = leftBMin[2];

  node.bmax[0][0] = leftBMax[0];
  node.bmax[0][1] = leftBMax[1];
  node.bmax[0][2] = leftBMax[2];

  node.bmin[1][0] = rightBMin[0];
  node.bmin[1][1] = rightBMin[1];
  node.bmin[1][2] = rightBMin[2];

  node.bmax[1][0] = rightBMax[0];
  node.bmax[1][1] = rightBMax[1];
  node.bmax[1][2] = rightBMax[2];

#pragma omp critical
  {
    nodes_[offset] = node;
    stats_.numBranchNodes++;
  }

  return offset;
}
#else
size_t ParticleAccel::BuildTree(const Particles *particles, real3 &bmin,
                                real3 &bmax, unsigned int leftIdx,
                                unsigned int rightIdx, int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  size_t offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
    stats_.maxTreeDepth = depth;
  }

  ComputeBoundingBox(bmin, bmax, particles, &indices_.at(0), leftIdx, rightIdx);

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {

    // Create leaf node. `n' may be null(create empty leaf node for that case.)
    ParticleNode leaf;

    leaf.bmin[0][0] = bmin[0];
    leaf.bmin[0][1] = bmin[1];
    leaf.bmin[0][2] = bmin[2];

    leaf.bmax[0][0] = bmax[0];
    leaf.bmax[0][1] = bmax[1];
    leaf.bmax[0][2] = bmax[2];

    assert(leftIdx < std::numeric_limits<unsigned int>::max());

    leaf.flag = 1; // leaf
    leaf.data[0] = n;
    leaf.data[1] = (unsigned int)leftIdx;

    nodes_.push_back(leaf);

    stats_.numLeafNodes++;

    return offset;
  }

  //
  // Create branch node.
  //

  //
  // Simple spatial median.
  //

  // Begin with longest axis
  real3 bsize;
  bsize[0] = bmax[0] - bmin[0];
  bsize[1] = bmax[1] - bmin[1];
  bsize[2] = bmax[2] - bmin[2];

  real maxsize = bsize[0];
  int longestAxis = 0;

  if (maxsize < bsize[1]) {
    maxsize = bsize[1];
    longestAxis = 1;
  }

  if (maxsize < bsize[2]) {
    maxsize = bsize[2];
    longestAxis = 2;
  }

  // Try all 3 axis until good cut position avaiable.
  unsigned int midIdx = 0;
  int cutAxis = 0;
  for (int axisTry = 0; axisTry < 1; axisTry++) {

    unsigned int *begin = &indices_[leftIdx];
    // unsigned int* end   = &indices_[rightIdx];
    unsigned int *end = begin + (rightIdx - leftIdx);
    unsigned int *mid = 0;

    // Try longest axis first.
    cutAxis = (longestAxis + axisTry) % 3;
    real cutPos = bmin[cutAxis] + 0.5 * bsize[cutAxis]; // spatial median

    //
    // Split at (cutAxis, cutPos)
    // indices_ will be modified.
    //
    mid = std::partition(begin, end, PartPred(cutAxis, cutPos, particles));

    midIdx = leftIdx + (mid - begin);
    if ((midIdx == leftIdx) || (midIdx == rightIdx)) {

      // Can't split well.
      // Switch to object median(which may create unoptimized tree, but stable)
      midIdx = leftIdx + (n >> 1);

      // Try another axis if there's axis to try.

    } else {

      // Found good cut. exit loop.
      break;
    }
  }

  ParticleNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch
  nodes_.push_back(node);

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;
  unsigned int leftChildIndex =
      BuildTree(particles, leftBMin, leftBMax, leftIdx, midIdx, depth + 1);
  unsigned int rightChildIndex =
      BuildTree(particles, rightBMin, rightBMax, midIdx, rightIdx, depth + 1);

  // refit.
  MergeBoundingBox(bmin, bmax, leftBMin, leftBMax, rightBMin, rightBMax);

  nodes_[offset].data[0] = leftChildIndex;
  nodes_[offset].data[1] = rightChildIndex;

  nodes_[offset].bmin[0][0] = leftBMin[0];
  nodes_[offset].bmin[0][1] = leftBMin[1];
  nodes_[offset].bmin[0][2] = leftBMin[2];

  nodes_[offset].bmax[0][0] = leftBMax[0];
  nodes_[offset].bmax[0][1] = leftBMax[1];
  nodes_[offset].bmax[0][2] = leftBMax[2];

  nodes_[offset].bmin[1][0] = rightBMin[0];
  nodes_[offset].bmin[1][1] = rightBMin[1];
  nodes_[offset].bmin[1][2] = rightBMin[2];

  nodes_[offset].bmax[1][0] = rightBMax[0];
  nodes_[offset].bmax[1][1] = rightBMax[1];
  nodes_[offset].bmax[1][2] = rightBMax[2];

  stats_.numBranchNodes++;

  return offset;
}
#endif // !_OPENMP

#if 0
size_t ParticleAccel::BuildTreeBlob(const Particles *particles,
                                    std::vector<unsigned int> &indices,
                                    unsigned int prevNum, int depth) {
  debug("d: %d, n: %d\n", depth, (int)indices.size());

  size_t offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
    stats_.maxTreeDepth = depth;
  }

  size_t n = indices.size();

  real3 bmin, bmax;
  ComputeBoundingBox(bmin, bmax, particles, &indices.at(0), 0, n);

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {
    // Create leaf node.
    ParticleNode leaf;

    leaf.bmin[0][0] = bmin[0];
    leaf.bmin[0][1] = bmin[1];
    leaf.bmin[0][2] = bmin[2];

    leaf.bmax[0][0] = bmax[0];
    leaf.bmax[0][1] = bmax[1];
    leaf.bmax[0][2] = bmax[2];

    assert(n < std::numeric_limits<unsigned int>::max());

    leaf.flag = 1; // leaf
    leaf.data[0] = n;
    leaf.data[1] = indices_.size();

    nodes_.push_back(leaf);

    // Add indices.
    indices_.insert(indices_.end(), indices.begin(), indices.end());

    stats_.numLeafNodes++;

    return offset;
  }

  //
  // Create branch node.
  //

  //
  // Simple spatial median.
  //

  // Begin with longest axis
  real3 bsize;
  bsize[0] = bmax[0] - bmin[0];
  bsize[1] = bmax[1] - bmin[1];
  bsize[2] = bmax[2] - bmin[2];

  real maxsize = bsize[0];
  int longestAxis = 0;

  if (maxsize < bsize[1]) {
    maxsize = bsize[1];
    longestAxis = 1;
  }

  if (maxsize < bsize[2]) {
    maxsize = bsize[2];
    longestAxis = 2;
  }

  // Try all 3 axis until good cut position avaiable.
  unsigned int midIdx = 0;
  int cutAxis = 0;
  for (int axisTry = 0; axisTry < 1; axisTry++) {

    // Try longest axis first.
    cutAxis = (longestAxis + axisTry) % 3;
    real cutPos = bmin[cutAxis] + 0.5 * bsize[cutAxis]; // spatial median

    //
    // Split at (cutAxis, cutPos)
    // indices will be modified.
    //
    std::vector<unsigned int>::iterator mid = std::partition(
        indices.begin(), indices.end(), PartPred(cutAxis, cutPos, particles));

    midIdx = (mid - indices.begin());
    if ((midIdx == 0) || (midIdx >= (n - 1))) {

      // Can't split well.
      // Switch to object median(which may create unoptimized tree, but stable)
      midIdx = (n >> 1);

      // Try another axis if there's axis to try.

    } else {

      // Found good cut. exit loop.
      break;
    }
  }

  // Generate new indices.
  // For blob, some blob overlappes both bbox.
  // Handle it by allowing multiple registration of index to both nodes.
  //
  std::vector<unsigned int> leftIndices;
  std::vector<unsigned int> rightIndices;
  leftIndices.insert(leftIndices.begin(), indices.begin(),
                     indices.begin() + midIdx);
  rightIndices.insert(rightIndices.begin(), indices.begin() + midIdx,
                      indices.end());

  debug("left Sz: %ld\n", leftIndices.size());
  debug("rightSz: %ld\n", rightIndices.size());
  debug("center = %f\n", bmin[cutAxis] + 0.5 * bsize[cutAxis]);

  {
    real3 leftBMin, leftBMax;
    real3 rightBMin, rightBMax;
    ComputeBoundingBox(leftBMin, leftBMax, particles, &leftIndices.at(0), 0,
                       leftIndices.size());
    ComputeBoundingBox(rightBMin, rightBMax, particles, &rightIndices.at(0), 0,
                       rightIndices.size());

    std::vector<unsigned int> rightOverlapped;
    std::vector<unsigned int> leftOverlapped;

    // Find blob which overlappes both node.
    for (size_t i = 0; i < leftIndices.size(); i++) {
      unsigned int idx = leftIndices[i];
      real p = particles->positions[3 * idx + cutAxis];
      real radius = GetRadius(particles, idx);

      debug("p+rad = %f, rightBmin[%d] = %f\n", p + radius, cutAxis,
            rightBMin[cutAxis]);
      if (p + radius >= rightBMin[cutAxis]) {
        debug("right got\n");
        rightOverlapped.push_back(idx);
      }
    }

    for (size_t i = 0; i < rightIndices.size(); i++) {
      unsigned int idx = rightIndices[i];
      real p = particles->positions[3 * idx + cutAxis];
      real radius = GetRadius(particles, idx);

      debug("p-rad = %f, leftBMax[%d] = %f\n", p - radius, cutAxis,
            leftBMax[cutAxis]);
      if (p - radius <= leftBMax[cutAxis]) {
        debug("left got\n");
        leftOverlapped.push_back(idx);
      }
    }

    debug("leftBMin: (%f, %f, %f) - (%f, %f, %f)\n", leftBMin[0], leftBMin[1],
          leftBMin[2], leftBMax[0], leftBMax[1], leftBMax[2]);
    debug("rightBMin: (%f, %f, %f) - (%f, %f, %f)\n", rightBMin[0],
          rightBMin[1], rightBMin[2], rightBMax[0], rightBMax[1], rightBMax[2]);
    debug("axis: %d\n", cutAxis);
    debug("indices.size = %ld\n", indices.size());
    debug("rightOverlap: %ld\n", rightOverlapped.size());
    debug("leftOverlap: %ld\n", leftOverlapped.size());

    leftIndices.insert(leftIndices.end(), leftOverlapped.begin(),
                       leftOverlapped.end());
    rightIndices.insert(rightIndices.end(), rightOverlapped.begin(),
                        rightOverlapped.end());
  }

  ParticleNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch
  nodes_.push_back(node);

  // Recurively split tree.
  unsigned int leftChildIndex =
      BuildTreeBlob(particles, leftIndices, indices.size(), depth + 1);
  unsigned int rightChildIndex =
      BuildTreeBlob(particles, rightIndices, indices.size(), depth + 1);

  nodes_[offset].data[0] = leftChildIndex;
  nodes_[offset].data[1] = rightChildIndex;

  nodes_[offset].bmin[0] = bmin[0];
  nodes_[offset].bmin[1] = bmin[1];
  nodes_[offset].bmin[2] = bmin[2];

  nodes_[offset].bmax[0] = bmax[0];
  nodes_[offset].bmax[1] = bmax[1];
  nodes_[offset].bmax[2] = bmax[2];

  stats_.numBranchNodes++;

  return offset;
}
#endif

bool ParticleAccel::Build(const Particles *particles,
                          const ParticleBuildOptions &options) {
  options_ = options;
  stats_ = ParticleBuildStatistics();

  assert(particles);

  // @todo { double precision }
  assert(particles->isDoublePrecisionPos == false);

  size_t n = particles->numParticles;
  trace("[ParticleAccel] Input # of particles = %lu\n", n);

#if USE_BLOB

  //
  // 1. Create initial temporal indices
  //
  std::vector<unsigned int> indices(n);
  for (size_t i = 0; i < n; i++) {
    indices[i] = i;
  }

  //
  // 2. Build tree
  //
  indices_.clear(); // indices_ will be generated in BuildTreeBlob
  BuildTreeBlob(particles, indices, indices.size(), 0);

#else

  //
  // 1. Create indices(this will be permutated in BuildTree)
  //
  timerutil t;
  t.start();
  indices_.resize(n);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < n; i++) {
    indices_[i] = i;
  }
  t.end();
  printf("[1:indexing] %d msec\n", (int)t.msec());

//
// 2. Build tree
//
//#ifdef _OPENMP
#if 0
#pragma omp parallel
  {
    // Task parallel.
#pragma omp single
    BuildTree(particles, 0, n, 0);
  }
#else
  real3 bmin, bmax;
  BuildTree(particles, bmin, bmax, 0, n, 0);
#endif

#endif

  // Store bbox
  bmin_[0] = bmin[0];
  bmin_[1] = bmin[1];
  bmin_[2] = bmin[2];

  bmax_[0] = bmax[0];
  bmax_[1] = bmax[1];
  bmax_[2] = bmax[2];

  // Tree will be null if input primitive count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    // real3 bmin(&nodes_[0].bmin[0]);
    // real3 bmax(&nodes_[0].bmax[0]);
    trace("[ParticleAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1],
          bmin[2]);
    trace("[ParticleAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1],
          bmax[2]);
  }

  trace("[ParticleAccel] # of nodes = %lu\n", nodes_.size());

  // Store pointer to the primitive
  particles_ = particles;

  return true;
}

bool ParticleAccel::Build32(const Particles *particles,
                            const ParticleBuildOptions &options) {
  options_ = options;
  stats_ = ParticleBuildStatistics();

  assert(particles);

  // @todo { double precision }
  assert(particles->isDoublePrecisionPos == false);

  size_t n = particles->numParticles;

  if (n < 1024) {
    // Use non-optimized accel builder
    return Build(particles, options);
  }

  trace("[ParticleAccel] Input # of particles = %lu\n", n);

#if USE_BLOB
  assert(0);
#endif

  //
  // 1. Create indices(this will be permutated in BuildTree)
  //
  timerutil t;
  t.start();
  indices_.resize(n);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < n; i++) {
    indices_[i] = i;
  }
  t.end();
  trace("[1:indexing] %d msec\n", (int)t.msec());

  {
    timerutil t;
    t.start();

    real3 bmin, bmax;

#ifdef _OPENMP
    ComputeBoundingBoxOMP(bmin, bmax, particles, &indices_.at(0), 0, n);
#else
    ComputeBoundingBox(bmin, bmax, particles, &indices_.at(0), 0, n);
#endif

    t.end();
    trace("[2:scene bbox calculation] %d msec\n", (int)t.msec());

    trace(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    trace(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

    std::vector<uint32_t> codes(n);
    assert(sizeof(real) == sizeof(float));

    {
      timerutil t;
      t.start();

      CalculateMortonCodes30SIMD(&codes.at(0), particles->positions, bmin, bmax,
                                 0, n);
      t.end();
      trace("[3:morton calculation] %d msec\n", (int)t.msec());

      // for (size_t i = 0; i < n; i++) {
      //  printf("code[%d] = %d\n", i, codes[i]);
      //}
    }

    std::vector<IndexKey30> keys(n);

#ifdef _OPENMP
#pragma omp parallel for if (n > (1024 * 1024))
#endif
    for (size_t i = 0; i < n; i++) {
      keys[i].index = indices_[i];
      keys[i].code = codes[i];
    }

    // Free memory
    std::vector<uint32_t>().swap(indices_);
    std::vector<uint32_t>().swap(codes);

    { // sort
      timerutil t;
      t.start();
      std::vector<IndexKey30> temp(n);

#ifdef _OPENMP
      RadixSort30OMP(&keys.at(0), &keys.at(0) + n);
#pragma omp barrier
      t.end();
      trace("[4:radix sort] %d msec\n", (int)t.msec());

#else

      RadixSortByMortionCode30LSB(&keys.at(0), &keys.at(0) + n);

      t.end();
      trace("[4:sort by Morton codes LSB] %ld msec\n", t.msec());
#endif

      //// CHECK
      // for (size_t i = 0; i < n-1; i++) {
      //  //std::string b = BitString32(keys[i].code);
      //  //printf("[%08d] i = %010d, c = %010d(%s)\n", i, keys[i].index,
      // keys[i].code, b.c_str());
      //  assert(keys[i].code <= keys[i+1].code);
      //}
      // printf("check ok\n");
    }

    std::vector<NodeInfo32> nodeInfos(n - 1);
    {
      timerutil t;
      t.start();

#ifdef _OPENMP
#pragma omp parallel for if (n > (1024 * 1024))
#endif
      for (size_t i = 0; i < n - 1; i++) {
        nodeInfos[i] = ConstructBinaryRadixTree30(&keys.at(0), i, n);
        // printf("I[%d].index   = %d\n", i, nodeInfos[i].index);
        // printf("I[%d].leftTy  = %d\n", i, nodeInfos[i].leftType);
        // printf("I[%d].rightTy = %d\n", i, nodeInfos[i].rightType);
      }

      t.end();
      trace("[5:Construct binary radix tree: %d msec\n", (int)t.msec());
    }

    {
      timerutil t;
      t.start();

      nodes_.clear();

      // Explicitly create root node and reserve storage here.
      ParticleNode rootNode;
      nodes_.push_back(rootNode);

      bool isLeftLeaf =
          (nodeInfos[0].leftType == NODE_TYPE_LEAF) ? true : false;
      bool isRightLeaf =
          (nodeInfos[0].rightType == NODE_TYPE_LEAF) ? true : false;
      uint32_t midIndex = nodeInfos[0].childIndex;

      real3 leftBMin, leftBMax;
      real3 rightBMin, rightBMax;

// printf("root: midIndex = %d, range (%d, %d), flag = %d/%d\n", midIndex,
// 0, n-1, isLeftLeaf, isRightLeaf);

//#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10
//  // OMP version sometimes doesn't work. Use sequential code for a while.
// size_t leftChildIndex =
//    BuildTreeRecursive32(nodes_, leftBMin, leftBMax, keys, nodeInfos,
//                         particles, midIndex, 0, midIndex, isLeftLeaf, 0);
// size_t rightChildIndex = BuildTreeRecursive32(
//    nodes_, rightBMin, rightBMax, keys, nodeInfos, particles,
//    midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);

//#elif defined(_OPENMP)
#if defined(_OPENMP)

      size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

#pragma omp parallel shared(leftChildIndex, rightChildIndex, keys, nodeInfos,  \
                            particles)
      {
#pragma omp single
        {
          leftChildIndex = BuildTreeRecursive32OMP(
              nodes_, leftBMin, leftBMax, keys, nodeInfos, particles, midIndex,
              0, midIndex, isLeftLeaf, 0);
          // leftChildIndex = BuildTreeRecursive32(nodes_, leftBMin, leftBMax,
          // keys, nodeInfos,
          //                         particles, midIndex, 0, midIndex,
          // isLeftLeaf, 0);
        }

#pragma omp single
        {
          rightChildIndex = BuildTreeRecursive32OMP(
              nodes_, rightBMin, rightBMax, keys, nodeInfos, particles,
              midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);
          // rightChildIndex = BuildTreeRecursive32(nodes_, rightBMin,
          // rightBMax, keys, nodeInfos, particles,
          //    midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);
          // printf("right = %d, nodes = %d\n", rightChildIndex, nodes_.size());
        }
      }
#pragma omp barrier

      assert(leftChildIndex != (size_t)(-1));
      assert(rightChildIndex != (size_t)(-1));

#else

      size_t leftChildIndex =
          BuildTreeRecursive32(nodes_, leftBMin, leftBMax, keys, nodeInfos,
                               particles, midIndex, 0, midIndex, isLeftLeaf, 0);
      size_t rightChildIndex = BuildTreeRecursive32(
          nodes_, rightBMin, rightBMax, keys, nodeInfos, particles,
          midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);

#endif

      // printf("  leftbmin = %f, %f, %f\n", leftBMin[0], leftBMin[1],
      // leftBMin[2]);
      // printf("  leftbmax = %f, %f, %f\n", leftBMax[0], leftBMax[1],
      // leftBMax[2]);
      // printf("  rightbmin = %f, %f, %f\n", rightBMin[0], rightBMin[1],
      // rightBMin[2]);
      // printf("  rightbmax = %f, %f, %f\n", rightBMax[0], rightBMax[1],
      // rightBMax[2]);
      // printf("  leftIndex = %d, rightIndex = %d\n", leftChildIndex,
      // rightChildIndex);

      real3 rootBMin, rootBMax;

      MergeBoundingBox(rootBMin, rootBMax, leftBMin, leftBMax, rightBMin,
                       rightBMax);

      rootNode.axis = GetSplitAxis(keys[0].code);
      rootNode.flag = 0; // = branch

      rootNode.data[0] = leftChildIndex;
      rootNode.data[1] = rightChildIndex;

      rootNode.bmin[0][0] = leftBMin[0];
      rootNode.bmin[0][1] = leftBMin[1];
      rootNode.bmin[0][2] = leftBMin[2];

      rootNode.bmax[0][0] = leftBMax[0];
      rootNode.bmax[0][1] = leftBMax[1];
      rootNode.bmax[0][2] = leftBMax[2];

      rootNode.bmin[1][0] = rightBMin[0];
      rootNode.bmin[1][1] = rightBMin[1];
      rootNode.bmin[1][2] = rightBMin[2];

      rootNode.bmax[1][0] = rightBMax[0];
      rootNode.bmax[1][1] = rightBMax[1];
      rootNode.bmax[1][2] = rightBMax[2];

      // @atomic update
      nodes_[0] = rootNode;

      t.end();
      // printf("root: midIndex = %d, range (%d, %d), flag = %d/%d\n", midIndex,
      // 0, n-1, isLeftLeaf, isRightLeaf);
      // printf("node.total = %d", nodes_.size());
      // printf("[%d] -> (%d, %d)\n", 0, leftChildIndex, rightChildIndex);
      trace("[6:Construct final AABB tree: %d msec\n", (int)t.msec());

      printf("  bmin = %f, %f, %f\n", rootBMin[0], rootBMin[1], rootBMin[2]);
      printf("  bmax = %f, %f, %f\n", rootBMax[0], rootBMax[1], rootBMax[2]);

      bmin_[0] = rootBMin[0];
      bmin_[1] = rootBMin[1];
      bmin_[2] = rootBMin[2];

      bmax_[0] = rootBMax[0];
      bmax_[1] = rootBMax[1];
      bmax_[2] = rootBMax[2];
    }

    std::vector<NodeInfo32>().swap(nodeInfos); // clear memory

    {
      indices_.resize(keys.size());

      // Store sorted indices.
      assert(indices_.size() == keys.size());
      for (size_t i = 0; i < keys.size(); i++) {
        indices_[i] = keys[i].index;
      }
    }
  }

  // Tree will be null if input primitive count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    // real3 bmin(&nodes_[0].bmin[0]);
    // real3 bmax(&nodes_[0].bmax[0]);
    trace("[ParticleAccel] bound min = (%f, %f, %f)\n", bmin_[0], bmin_[1],
          bmin_[2]);
    trace("[ParticleAccel] bound max = (%f, %f, %f)\n", bmax_[0], bmax_[1],
          bmax_[2]);
  }

  trace("[ParticleAccel] # of nodes = %lu\n", nodes_.size());

  trace("  ParticleAccel : %llu MB\n",
        sizeof(ParticleNode) * nodes_.size() / (1024ULL * 1014ULL));

  // Store pointer
  particles_ = particles;

  return true;
}

bool ParticleAccel::Dump(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "[ParticleAccel] Cannot write a file: %s\n", filename);
    return false;
  }

  unsigned long long numNodes = nodes_.size();
  assert(nodes_.size() > 0);

  unsigned long long numIndices = indices_.size();

  int r = 0;
  r = fwrite(&numNodes, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&nodes_.at(0), sizeof(ParticleNode), numNodes, fp);
  assert(r == numNodes);

  r = fwrite(&numIndices, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&indices_.at(0), sizeof(unsigned int), numIndices, fp);
  assert(r == numIndices);

  fclose(fp);

  return true;
}

namespace {

const int kMaxStackDepth = 512;

#if ENABLE_SIMD_ISECTOR

#if defined(__sparc__) && defined(__HPC_ACE__)

// 1 ray - 2 AABB Intersection test.
inline void IntersectRayAABBD2(
    uint64_t retMask[2],         // output
    const double2 bboxes[2][3],  // 2boxes : min-max[2] of xyz[3] of boxes[4]
    const double2 rox,           // ray origin
    const double2 roy,           // ray origin
    const double2 roz,           // ray origin
    const double2 ridx,          // ray inverse direction
    const double2 ridy,          // ray inverse direction
    const double2 ridz,          // ray inverse direction
    const double2 roidx,         // org * invdir
    const double2 roidy,         //
    const double2 roidz,         //
    const int sign[3],           // ray xyz direction -> +:0,-:1
    double2 vtmin, double2 vtmax // ray range tmin-tmax
    ) {
  // x coordinate
  double2 tmin =
      vmaxd2(vtmin, _mm_sub_pd(_mm_mul_pd(bboxes[sign[0]][0], ridx), roidx));
  double2 tmax = vmind2(
      vtmax, _mm_sub_pd(_mm_mul_pd(bboxes[1 - sign[0]][0], ridx), roidx));

  // y coordinate
  tmin = vmaxd2(tmin, _mm_sub_pd(_mm_mul_pd(bboxes[sign[1]][1], ridy), roidy));
  tmax =
      vmind2(tmax, _mm_sub_pd(_mm_mul_pd(bboxes[1 - sign[1]][1], ridy), roidy));

  // z coordinate
  tmin = vmaxd2(tmin, _mm_sub_pd(_mm_mul_pd(bboxes[sign[2]][2], ridz), roidz));
  tmax =
      vmind2(tmax, _mm_sub_pd(_mm_mul_pd(bboxes[1 - sign[2]][2], ridz), roidz));

  double2 vret = _mm_cmpge_pd(tmax, tmin);

  // double ret[2];
  _mm_store_pd(reinterpret_cast<double *>(retMask), vret);

  // mimics _mm_movemask_pd
  // return (ret[0] ? 1 : 0) | ((ret[1] ? 1 : 0) << 1);
}

#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

// 1 ray - 2 AABB Intersection test.
inline int IntersectRayAABBD2(
    const double2 bboxes[2][3],  // 2boxes : min-max[2] of xyz[3] of boxes[2]
    const double2 rox,           // ray origin
    const double2 roy,           // ray origin
    const double2 roz,           // ray origin
    const double2 ridx,          // ray inverse direction
    const double2 ridy,          // ray inverse direction
    const double2 ridz,          // ray inverse direction
    const double2 roidx,         // org * invdir
    const double2 roidy,         //
    const double2 roidz,         //
    const int sign[3],           // ray xyz direction -> +:0,-:1
    double2 vtmin, double2 vtmax // ray range tmin-tmax
    ) {
  // x coordinate
  double2 tmin =
      vmaxd2(vtmin, _mm_sub_pd(_mm_mul_pd(bboxes[sign[0]][0], ridx), roidx));
  double2 tmax = vmind2(
      vtmax, _mm_sub_pd(_mm_mul_pd(bboxes[1 - sign[0]][0], ridx), roidx));

  // y coordinate
  tmin = vmaxd2(tmin, _mm_sub_pd(_mm_mul_pd(bboxes[sign[1]][1], ridy), roidy));
  tmax =
      vmind2(tmax, _mm_sub_pd(_mm_mul_pd(bboxes[1 - sign[1]][1], ridy), roidy));

  // z coordinate
  tmin = vmaxd2(tmin, _mm_sub_pd(_mm_mul_pd(bboxes[sign[2]][2], ridz), roidz));
  tmax =
      vmind2(tmax, _mm_sub_pd(_mm_mul_pd(bboxes[1 - sign[2]][2], ridz), roidz));

  double2 vret =
      _mm_and_pd(_mm_cmpge_pd(tmax, vzerod2()), _mm_cmpge_pd(tmax, tmin));

  return _mm_movemask_pd(vret);
}
#else
#error SIMD is not supported on this architecture. Disable ENABLE_SIMD_ISECTOR to solve this error.
#endif
#endif // ENABLE_SIMD_ISECTOR

inline bool IntersectRayAABB(real &tminOut, // [out]
                             real &tmaxOut, // [out]
                             real maxT, real const bmin[3], real const bmax[3],
                             real3 rayOrg, real3 rayInvDir, int rayDirSign[3]) {
  real tmin, tmax;

  const real min_x = rayDirSign[0] ? bmax[0] : bmin[0];
  const real min_y = rayDirSign[1] ? bmax[1] : bmin[1];
  const real min_z = rayDirSign[2] ? bmax[2] : bmin[2];
  const real max_x = rayDirSign[0] ? bmin[0] : bmax[0];
  const real max_y = rayDirSign[1] ? bmin[1] : bmax[1];
  const real max_z = rayDirSign[2] ? bmin[2] : bmax[2];

  // X
  const double tmin_x = (min_x - rayOrg[0]) * rayInvDir[0];
  const double tmax_x = (max_x - rayOrg[0]) * rayInvDir[0];

  // Y
  const double tmin_y = (min_y - rayOrg[1]) * rayInvDir[1];
  const double tmax_y = (max_y - rayOrg[1]) * rayInvDir[1];

  tmin = (tmin_x > tmin_y) ? tmin_x : tmin_y;
  tmax = (tmax_x < tmax_y) ? tmax_x : tmax_y;

  // Z
  const double tmin_z = (min_z - rayOrg[2]) * rayInvDir[2];
  const double tmax_z = (max_z - rayOrg[2]) * rayInvDir[2];

  tmin = (tmin > tmin_z) ? tmin : tmin_z;
  tmax = (tmax < tmax_z) ? tmax : tmax_z;

  tmax = std::min(tmax, maxT);

  //
  // Hit include (tmin == tmax) edge case(hit 2D plane).
  //
  if ((tmax > 0.0) && (tmin <= tmax)) {

    tminOut = tmin;
    tmaxOut = tmax;

    return true;
  }

  return false; // no hit
}

inline real vdot(real3 a, real3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline real3 vnormalize(real3 a) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  real3 v;
  real len = sqrt(vdot(a, a));
  v = a;

  if (len > kEPS) {
    real invlen = 1.0 / len;
    v[0] *= invlen;
    v[1] *= invlen;
    v[2] *= invlen;
  }

  return v;
}

#if ENABLE_SIMD_ISECTOR
#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10

//
// Do 1 ray - 2 spheres SIMD intersection test.
//
inline void SphereIsectD2(double2 &vtInOut,
                          double2 &vtidInOut, // store integer(53bit)
                          const double2 &vtid, const double2 &vcx, // xx
                          const double2 &vcy,                      // yy
                          const double2 &vcz,                      // zz
                          const double2 &vr,                       // rr
                          const double2 &rox,                      // rayOrg.xx
                          const double2 &roy,                      // rayOrg.yy
                          const double2 &roz,                      // rayOrg.zz
                          const double2 &rdx,                      // rayDir.xx
                          const double2 &rdy,                      // rayDir.yy
                          const double2 &rdz,                      // rayDir.zz
                          const double2 &va) // vdot(rayDir, rayDir).xx;
{
  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;
  const double2 veps = _mm_set_pd(kEPS, kEPS);

  // real3 center(v[0], v[1], v[2]);
  // real3 oc = rayOrg - center;
  double2 ocx = _mm_sub_pd(rox, vcx);
  double2 ocy = _mm_sub_pd(roy, vcy);
  double2 ocz = _mm_sub_pd(roz, vcz);

  double2 vzero = vzerod2();
  double2 vtwo = _mm_set_pd(2.0, 2.0);
  double2 vfour = _mm_set_pd(4.0, 4.0);

  // real b = 2.0f * vdot(rayDir, oc);
  double2 vb = _mm_mul_pd(vtwo, vdotd2(rdx, rdy, rdz, ocx, ocy, ocz));

  // real c = vdot(oc, oc) - radius * radius;
  double2 vc =
      _mm_sub_pd(vdotd2(ocx, ocy, ocz, ocx, ocy, ocz), _mm_mul_pd(vr, vr));

  // real disc = b * b - 4.0 * a * c;
  double2 vdisc =
      _mm_sub_pd(_mm_mul_pd(vb, vb), _mm_mul_pd(vfour, _mm_mul_pd(va, vc)));

  // if (disc > 0)
  double2 vrootMask = _mm_cmpge_pd(vdisc, vzerod2());

  //{ // early exit test
  //  double tRet[2];

  //  _mm_store_pd(tRet, vrootMask);

  //  const uint64_t *iRet = reinterpret_cast<uint64_t*>(tRet);

  //  if (iRet[0] && iRet[1]) {
  //    return;
  //  }
  //}

  //// compute q as described above
  // real distSqrt = sqrt(disc);
  // real q;
  // if (b < 0)
  //  q = (-b - distSqrt) / 2.0;
  // else
  //  q = (-b + distSqrt) / 2.0;
  double2 vdistSqrt = vsqrtd2(vdisc);
  double2 vcond_b = _mm_cmplt_pd(vb, vzero);
  double2 vhalf = _mm_set_pd(0.5, 0.5);

  // (b <0) ? -distSqrt : distSqrt
  double2 vdistSqrt_s =
      vseld2(_mm_sub_pd(vzero, vdistSqrt), vdistSqrt, vcond_b);
  double2 vq =
      _mm_mul_pd(_mm_add_pd(_mm_sub_pd(vzero, vb), vdistSqrt_s), vhalf);

  //// compute t0 and t1
  // real t0 = q / a;
  // real t1 = c / q;
  double2 vt0_0 = _mm_mul_pd(vq, vinvd2(va));
  double2 vt1_0 = _mm_mul_pd(vc, vinvd2(vq));

  //// make sure t0 is smaller than t1
  // if (t0 > t1) {
  //  // if t0 is bigger than t1 swap them around
  //  real temp = t0;
  //  t0 = t1;
  //  t1 = temp;
  //}
  double2 vt0 = vmind2(vt0_0, vt1_0);
  double2 vt1 = vmaxd2(vt0_0, vt1_0);

  //// if t1 is less than zero, the object is in the ray's negative direction
  //// and consequently the ray misses the sphere
  // if (t1 < 0) {
  //  return false;
  //}
  double2 vt1_mask = _mm_cmpge_pd(vt1, vzero);

  //// if t0 is less than zero, the intersection point is at t1
  // if (t0 < 0) {
  //  real t = t1;
  //  if (t < tInOut) {
  //    tInOut = t;
  //    return true;
  //  }
  //}
  //// else the intersection point is at t0
  // else {
  //  real t = t0;
  //  if (t < tInOut) {
  //    tInOut = t;
  //    return true;
  //  }
  //}
  double2 vt0_mask = _mm_cmplt_pd(vt0, vzero);
  double2 vt = vseld2(vt1, vt0, vt0_mask);

  double2 vtprev_mask = _mm_cmplt_pd(vt, vtInOut);
  double2 vmask0 = _mm_and_pd(vtprev_mask, vt1_mask);
  double2 vmask_final = _mm_and_pd(vmask0, vrootMask);

  vtidInOut = vseld2(vtid, vtidInOut, vmask_final);
  vtInOut = vseld2(vt, vtInOut, vmask_final);

  return;
}

#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

//
// Do 1 ray - 2 spheres SIMD intersection test.
//
inline void SphereIsectD2(double2 &vtInOut,
                          double2 &vtidInOut, // store integer(53bit)
                          const double2 &vtid, const double2 &vcx, // xx
                          const double2 &vcy,                      // yy
                          const double2 &vcz,                      // zz
                          const double2 &vr,                       // rr
                          const double2 &rox,                      // rayOrg.xx
                          const double2 &roy,                      // rayOrg.yy
                          const double2 &roz,                      // rayOrg.zz
                          const double2 &rdx,                      // rayDir.xx
                          const double2 &rdy,                      // rayDir.yy
                          const double2 &rdz,                      // rayDir.zz
                          const double2 &va) // vdot(rayDir, rayDir).xx;
{
  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;
  const double2 veps = _mm_set_pd(kEPS, kEPS);

  // real3 center(v[0], v[1], v[2]);
  // real3 oc = rayOrg - center;
  double2 ocx = _mm_sub_pd(rox, vcx);
  double2 ocy = _mm_sub_pd(roy, vcy);
  double2 ocz = _mm_sub_pd(roz, vcz);

  double2 vzero = vzerod2();
  double2 vtwo = _mm_set_pd(2.0, 2.0);
  double2 vfour = _mm_set_pd(4.0, 4.0);

  // real b = 2.0f * vdot(rayDir, oc);
  double2 vb = _mm_mul_pd(vtwo, vdotd2(rdx, rdy, rdz, ocx, ocy, ocz));

  // real c = vdot(oc, oc) - radius * radius;
  double2 vc =
      _mm_sub_pd(vdotd2(ocx, ocy, ocz, ocx, ocy, ocz), _mm_mul_pd(vr, vr));

  // real disc = b * b - 4.0 * a * c;
  double2 vdisc =
      _mm_sub_pd(_mm_mul_pd(vb, vb), _mm_mul_pd(vfour, _mm_mul_pd(va, vc)));

  // if (disc > 0)
  double2 vrootMask = _mm_cmpge_pd(vdisc, vzerod2());

  { // early exit test
    if (_mm_movemask_pd(vrootMask) == 0) {
      return;
    }
  }

  //// compute q as described above
  // real distSqrt = sqrt(disc);
  // real q;
  // if (b < 0)
  //  q = (-b - distSqrt) / 2.0;
  // else
  //  q = (-b + distSqrt) / 2.0;
  double2 vdistSqrt = vsqrtd2(vdisc);
  double2 vcond_b = _mm_cmplt_pd(vb, vzero);
  double2 vhalf = _mm_set_pd(0.5, 0.5);

  // (b <0) ? -distSqrt : distSqrt
  double2 vdistSqrt_s =
      vseld2(_mm_sub_pd(vzero, vdistSqrt), vdistSqrt, vcond_b);
  double2 vq =
      _mm_mul_pd(_mm_add_pd(_mm_sub_pd(vzero, vb), vdistSqrt_s), vhalf);

  //// compute t0 and t1
  // real t0 = q / a;
  // real t1 = c / q;
  double2 vt0_0 = _mm_mul_pd(vq, vinvd2(va));
  double2 vt1_0 = _mm_mul_pd(vc, vinvd2(vq));

  //// make sure t0 is smaller than t1
  // if (t0 > t1) {
  //  // if t0 is bigger than t1 swap them around
  //  real temp = t0;
  //  t0 = t1;
  //  t1 = temp;
  //}
  double2 vt0 = vmind2(vt0_0, vt1_0);
  double2 vt1 = vmaxd2(vt0_0, vt1_0);

  //// if t1 is less than zero, the object is in the ray's negative direction
  //// and consequently the ray misses the sphere
  // if (t1 < 0) {
  //  return false;
  //}
  double2 vt1_mask = _mm_cmpge_pd(vt1, vzero);

  //// if t0 is less than zero, the intersection point is at t1
  // if (t0 < 0) {
  //  real t = t1;
  //  if (t < tInOut) {
  //    tInOut = t;
  //    return true;
  //  }
  //}
  //// else the intersection point is at t0
  // else {
  //  real t = t0;
  //  if (t < tInOut) {
  //    tInOut = t;
  //    return true;
  //  }
  //}
  double2 vt0_mask = _mm_cmplt_pd(vt0, vzero);
  double2 vt = vseld2(vt1, vt0, vt0_mask);

  double2 vtprev_mask = _mm_cmplt_pd(vt, vtInOut);
  double2 vmask0 = _mm_and_pd(vtprev_mask, vt1_mask);
  double2 vmask_final = _mm_and_pd(vmask0, vrootMask);

  vtidInOut = vseld2(vtid, vtidInOut, vmask_final);
  vtInOut = vseld2(vt, vtInOut, vmask_final);

  return;
}
#else
#error SIMD is not supported on this architecture. Disable ENABLE_SIMD_ISECTOR to solve this error.
////
//// Do 1 ray - 4 triangles SIMD intersection test.
////
// inline int4 SphereIsect4(float4 &vtInOut, int4 &vtidInOut, const int4 &vtid,
//                         const float4 &vcx, // xxxx
//                         const float4 &vcy, // yyyy
//                         const float4 &vcz, // zzzz
//                         const float4 &vr,  // rrrr
//                         const float4 &rox, // rayOrg.xxxx
//                         const float4 &roy, // rayOrg.yyyy
//                         const float4 &roz, // rayOrg.zzzz
//                         const float4 &rdx, // rayDir.xxxx
//                         const float4 &rdy, // rayDir.yyyy
//                         const float4 &rdz, // rayDir.zzzz
//                         const float4 &va)  // vdot(rayDir, rayDir).xxxx;
//{
//  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;
//  const float4 veps = vset1_f4(kEPS);
//
//  // real3 center(v[0], v[1], v[2]);
//  // real3 oc = rayOrg - center;
//  float4 ocx = vsub_f4(rox, vcx);
//  float4 ocy = vsub_f4(roy, vcy);
//  float4 ocz = vsub_f4(roz, vcz);
//
//  // real a = vdot(rayDir, rayDir);
//
//  float4 vzero = vset1_f4(0.0f);
//  float4 vtwo = vset1_f4(2.0f);
//  float4 vfour = vset1_f4(4.0f);
//  // real b = 2.0f * vdot(rayDir, oc);
//  float4 vb = vmul_f4(vtwo, vdot_f4(rdx, rdy, rdz, ocx, ocy, ocz));
//  // real c = vdot(oc, oc) - radius * radius;
//  float4 vc = vsub_f4(vdot_f4(ocx, ocy, ocz, ocx, ocy, ocz), vmul_f4(vr, vr));
//
//  // real disc = b * b - 4.0 * a * c;
//  float4 vdisc = vsub_f4(vmul_f4(vb, vb), vmul_f4(vfour, vmul_f4(va, vc)));
//
//  int mask = 0;
//
//  assert(0); // @todo
//
//  // if (disc > 0)
//  float4 vrootMask = vcmpge_f4(vdisc, vzero_f4());
//
//  //// compute q as described above
//  // real distSqrt = sqrt(disc);
//  // real q;
//  // if (b < 0)
//  //  q = (-b - distSqrt) / 2.0;
//  // else
//  //  q = (-b + distSqrt) / 2.0;
//  float4 vdistSqrt = vfastsqrt_f4(vdisc);
//  float4 vcond_b = vcmplt_f4(vb, vzero);
//  float4 vhalf = vset1_f4(0.5f);
//
//  // (b <0) ? -distSqrt : distSqrt
//  float4 vdistSqrt_s = vsel_f4(vsub_f4(vzero, vdistSqrt), vdistSqrt, vcond_b);
//  float4 vq = vmul_f4(vadd_f4(vsub_f4(vzero, vb), vdistSqrt_s), vhalf);
//
//  //// compute t0 and t1
//  // real t0 = q / a;
//  // real t1 = c / q;
//  float4 vt0_0 = vmul_f4(vq, vfastinv_f4(va));
//  float4 vt1_0 = vmul_f4(vc, vfastinv_f4(vq));
//
//  //// make sure t0 is smaller than t1
//  // if (t0 > t1) {
//  //  // if t0 is bigger than t1 swap them around
//  //  real temp = t0;
//  //  t0 = t1;
//  //  t1 = temp;
//  //}
//  float4 vt0 = vmin_f4(vt0_0, vt1_0);
//  float4 vt1 = vmax_f4(vt0_0, vt1_0);
//
//  //// if t1 is less than zero, the object is in the ray's negative direction
//  //// and consequently the ray misses the sphere
//  // if (t1 < 0) {
//  //  return false;
//  //}
//  float4 vt1_mask = vcmpge_f4(vzero, vt1);
//
//  //// if t0 is less than zero, the intersection point is at t1
//  // if (t0 < 0) {
//  //  real t = t1;
//  //  if (t < tInOut) {
//  //    tInOut = t;
//  //    return true;
//  //  }
//  //}
//  //// else the intersection point is at t0
//  // else {
//  //  real t = t0;
//  //  if (t < tInOut) {
//  //    tInOut = t;
//  //    return true;
//  //  }
//  //}
//  float4 vt0_mask = vcmplt_f4(vt0, vzero);
//  float4 vt = vsel_f4(vt1, vt0, vt0_mask);
//
//  float4 vtprev_mask = vcmplt_f4(vt, vtInOut);
//  float4 vmask0 = vand_f4(vtprev_mask, vt1_mask);
//  float4 vmask_final = vand_f4(vmask0, vrootMask);
//
//  int4 imask = vcasttoi4_f4(vmask_final); // conversion-free cast
//  vtidInOut = vsel_i4(vtid, vtidInOut, imask);
//
//  return imask;
//}
#endif
#endif

// http://wiki.cgsociety.org/index.php/Ray_Sphere_Intersection
inline bool SphereIsect(real &tInOut, const real *v, real radius,
                        const real3 &rayOrg, const real3 &rayDir) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  real3 center(v[0], v[1], v[2]);
  real3 oc = rayOrg - center;

  real a = vdot(rayDir, rayDir);
  real b = 2.0 * vdot(rayDir, oc);
  real c = vdot(oc, oc) - radius * radius;

  real disc = b * b - 4.0 * a * c;

  if (disc < 0) { // no roots
    return false;
  }

  // compute q as described above
  real distSqrt = sqrt(disc);
  real q;
  if (b < 0)
    q = (-b - distSqrt) / 2.0;
  else
    q = (-b + distSqrt) / 2.0;

  // compute t0 and t1
  real t0 = q / a;
  real t1 = c / q;

  // make sure t0 is smaller than t1
  if (t0 > t1) {
    // if t0 is bigger than t1 swap them around
    real temp = t0;
    t0 = t1;
    t1 = temp;
  }

  // if t1 is less than zero, the object is in the ray's negative direction
  // and consequently the ray misses the sphere
  if (t1 < 0) {
    return false;
  }

  // if t0 is less than zero, the intersection point is at t1
  if (t0 < 0) {
    real t = t1;
    if (t < tInOut) {
      tInOut = t;
      return true;
    }
  }
  // else the intersection point is at t0
  else {
    real t = t0;
    if (t < tInOut) {
      tInOut = t;
      return true;
    }
  }

  return false;
}

bool TestLeafNode(Intersection &isect, // [inout]
                  const ParticleNode &node,
                  const std::vector<unsigned int> &indices,
                  const Particles *particles, const Ray &ray) {

  bool hit = false;

  unsigned int numPoints = node.data[0];
  unsigned int offset = node.data[1];

  real t = isect.t; // current hit distance

  real3 rayorg;
  rayorg[0] = ray.origin()[0];
  rayorg[1] = ray.origin()[1];
  rayorg[2] = ray.origin()[2];

  real3 raydir;
  raydir[0] = ray.direction()[0];
  raydir[1] = ray.direction()[1];
  raydir[2] = ray.direction()[2];

#if ENABLE_SIMD_ISECTOR

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10

#define N (4) // Unroll factor

  int lanes = 2 * N;

  uint32_t numBlocks = numPoints / lanes;

  // Widen data from float precision to double precision since
  // HPC-ACE is optimized for double and float speed is identical to double.
  const double2 rox = _mm_set_pd(ray.origin()[0], ray.origin()[0]);
  const double2 roy = _mm_set_pd(ray.origin()[1], ray.origin()[1]);
  const double2 roz = _mm_set_pd(ray.origin()[2], ray.origin()[2]);

  const double2 rdx = _mm_set_pd(ray.direction()[0], ray.direction()[0]);
  const double2 rdy = _mm_set_pd(ray.direction()[1], ray.direction()[1]);
  const double2 rdz = _mm_set_pd(ray.direction()[2], ray.direction()[2]);

  double vd = ray.direction()[0] * ray.direction()[0] +
              ray.direction()[1] * ray.direction()[1] +
              ray.direction()[2] * ray.direction()[2];
  const double2 va = _mm_set_pd(vd, vd);

  double2 tInOut[N];
  double2 tidInOut[N]; // Use as integer(53bit)

  const double kInf = std::numeric_limits<double>::infinity();

  // Init
  for (int k = 0; k < N; k++) {
    tInOut[k] = _mm_set_pd(kInf, kInf);   // = no hit
    tidInOut[k] = _mm_set_pd(-1.0, -1.0); // = no hit
  }

  for (unsigned int i = 0; i < lanes * numBlocks; i += lanes) {

    for (int k = 0; k < N; k++) {

      int partIdx0 = indices[i + 2 * k + 0 + offset];
      int partIdx1 = indices[i + 2 * k + 1 + offset];

      const real *p0 = &particles->positions[3 * partIdx0];
      const real *p1 = &particles->positions[3 * partIdx1];
      double radius0 = GetRadius(particles, partIdx0);
      double radius1 = GetRadius(particles, partIdx1);

      double2 vx = _mm_set_pd(p0[0], p1[0]);
      double2 vy = _mm_set_pd(p0[1], p1[1]);
      double2 vz = _mm_set_pd(p0[2], p1[2]);
      double2 vr = _mm_set_pd(radius0, radius1);

      // double has enough precision to store integer(53bit)
      double2 tid = _mm_set_pd(i + 2 * k + 0, i + 2 * k + 1);

      SphereIsectD2(tInOut[k], tidInOut[k], tid, vx, vy, vz, vr, rox, roy, roz,
                    rdx, rdy, rdz, va);

      // if (SphereIsect(t, p, radius, rayorg, raydir)) {
      //  // Update isect state
      //  isect.t = t;

      //  // u and v are computed in later.
      //  isect.prim_id = partIdx;
      //  hit = true;
      //}
    }
  }

  // Merge result.
  for (int i = 0; i < N; i++) {
    double tidRet[2];
    double tRet[2];

    _mm_store_pd(tidRet, tidInOut[i]);
    _mm_store_pd(tRet, tInOut[i]);

    for (int k = 0; k < 2; k++) {
      // swap index
      if ((tidRet[1 - k] >= 0.0) &&
          (tRet[1 - k] <= isect.t)) { // && ((unsigned int)(tidRet[k]) !=
                                      // ray.prev_prim_id)) { //
                                      // +self-intersection check
        // printf("hit. tid = %d\n", (unsigned int)tidRet[k]);
        isect.t = tRet[1 - k];
        isect.prim_id = indices[(unsigned int)tidRet[1 - k] + offset];
        hit = true;
        t = isect.t;
      }
    }
  }

  // Remainder
  {
    real3 rayOrg;
    rayOrg[0] = ray.origin()[0];
    rayOrg[1] = ray.origin()[1];
    rayOrg[2] = ray.origin()[2];

    real3 rayDir;
    rayDir[0] = ray.direction()[0];
    rayDir[1] = ray.direction()[1];
    rayDir[2] = ray.direction()[2];

    for (unsigned int i = lanes * numBlocks; i < numPoints; i++) {
      int partIdx = indices[i + offset];
      // assert(partIdx < particles->numParticles);
      // assert(partIdx >= 0);

      const real *p = &particles->positions[3 * partIdx];
      real radius = GetRadius(particles, partIdx);

      if (SphereIsect(t, p, radius, rayorg, raydir)) {
        // Update isect state
        isect.t = t;

        // u and v are computed in later.
        isect.prim_id = partIdx;
        hit = true;
      }
    }
  }
#undef N

#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

#define N (1) // Unroll factor

  int lanes = 2 * N;

  uint32_t numBlocks = numPoints / lanes;

  // Widen data from float precision to double precision since
  // HPC-ACE is optimized for double and float speed is identical to double.
  const double2 rox = _mm_set_pd(ray.origin()[0], ray.origin()[0]);
  const double2 roy = _mm_set_pd(ray.origin()[1], ray.origin()[1]);
  const double2 roz = _mm_set_pd(ray.origin()[2], ray.origin()[2]);

  const double2 rdx = _mm_set_pd(ray.direction()[0], ray.direction()[0]);
  const double2 rdy = _mm_set_pd(ray.direction()[1], ray.direction()[1]);
  const double2 rdz = _mm_set_pd(ray.direction()[2], ray.direction()[2]);

  double vd = ray.direction()[0] * ray.direction()[0] +
              ray.direction()[1] * ray.direction()[1] +
              ray.direction()[2] * ray.direction()[2];
  const double2 va = _mm_set_pd(vd, vd);

  double2 tInOut[N];
  double2 tidInOut[N]; // Use as integer(53bit)

  const double kInf = std::numeric_limits<double>::infinity();

  // Init
  for (int k = 0; k < N; k++) {
    tInOut[k] = _mm_set_pd(kInf, kInf);   // = no hit
    tidInOut[k] = _mm_set_pd(-1.0, -1.0); // = no hit
  }

  for (unsigned int i = 0; i < lanes * numBlocks; i += lanes) {

    for (int k = 0; k < N; k++) {

      int partIdx0 = indices[i + 2 * k + 0 + offset];
      int partIdx1 = indices[i + 2 * k + 1 + offset];

      const real *p0 = &particles->positions[3 * partIdx0];
      const real *p1 = &particles->positions[3 * partIdx1];
      double radius0 = GetRadius(particles, partIdx0);
      double radius1 = GetRadius(particles, partIdx1);

      double2 vx = _mm_set_pd(p0[0], p1[0]);
      double2 vy = _mm_set_pd(p0[1], p1[1]);
      double2 vz = _mm_set_pd(p0[2], p1[2]);
      double2 vr = _mm_set_pd(radius0, radius1);

      // double has enough precision to store integer(53bit)
      double2 tid = _mm_set_pd(i + 2 * k + 0, i + 2 * k + 1);

      SphereIsectD2(tInOut[k], tidInOut[k], tid, vx, vy, vz, vr, rox, roy, roz,
                    rdx, rdy, rdz, va);

      // if (SphereIsect(t, p, radius, rayorg, raydir)) {
      //  // Update isect state
      //  isect.t = t;

      //  // u and v are computed in later.
      //  isect.prim_id = partIdx;
      //  hit = true;
      //}
    }
  }

  // Merge result.
  for (int i = 0; i < N; i++) {
    double tidRet[2];
    double tRet[2];

    _mm_store_pd(tidRet, tidInOut[i]);
    _mm_store_pd(tRet, tInOut[i]);

    for (int k = 0; k < 2; k++) {
      // swap index
      if ((tidRet[1 - k] >= 0.0) &&
          (tRet[1 - k] <= isect.t)) { // && ((unsigned int)(tidRet[k]) !=
                                      // ray.prev_prim_id)) { //
                                      // +self-intersection check
        // printf("hit. tid = %d\n", (unsigned int)tidRet[k]);
        isect.t = tRet[1 - k];
        isect.prim_id = indices[(unsigned int)tidRet[1 - k] + offset];
        hit = true;
        t = isect.t;
      }
    }
  }

  // Remainder
  {
    real3 rayOrg;
    rayOrg[0] = ray.origin()[0];
    rayOrg[1] = ray.origin()[1];
    rayOrg[2] = ray.origin()[2];

    real3 rayDir;
    rayDir[0] = ray.direction()[0];
    rayDir[1] = ray.direction()[1];
    rayDir[2] = ray.direction()[2];

    for (unsigned int i = lanes * numBlocks; i < numPoints; i++) {
      int partIdx = indices[i + offset];
      // assert(partIdx < particles->numParticles);
      // assert(partIdx >= 0);

      const real *p = &particles->positions[3 * partIdx];
      real radius = GetRadius(particles, partIdx);

      if (SphereIsect(t, p, radius, rayorg, raydir)) {
        // Update isect state
        isect.t = t;

        // u and v are computed in later.
        isect.prim_id = partIdx;
        hit = true;
      }
    }
  }
#undef N

#else

  for (unsigned int i = 0; i < numPoints; i++) {
    int partIdx = indices[i + offset];
    // assert(partIdx < particles->numParticles);
    // assert(partIdx >= 0);

    const real *p = &particles->positions[3 * partIdx];
    real radius = GetRadius(particles, partIdx);

    if (SphereIsect(t, p, radius, rayorg, raydir)) {
      // Update isect state
      isect.t = t;

      // u and v are computed in later.
      isect.prim_id = partIdx;
      hit = true;
    }
  }
#endif

#else // Scalar

  for (unsigned int i = 0; i < numPoints; i++) {
    int partIdx = indices[i + offset];
    // assert(partIdx < particles->numParticles);
    // assert(partIdx >= 0);

    const real *p = &particles->positions[3 * partIdx];
    real radius = GetRadius(particles, partIdx);

    if (SphereIsect(t, p, radius, rayorg, raydir)) {
      // Update isect state
      isect.t = t;

      // u and v are computed in later.
      isect.prim_id = partIdx;
      hit = true;
    }
  }

#endif // SIMD

  return hit;
}

#if 0 // disable for a while
inline real CalcBlobDensity(real3 pos, const Particles *particles,
                            const std::vector<unsigned int> &indices,
                            unsigned int offset, int numBlobs, real threshold) {
  //
  // @todo { Optimize intersection calclulation. }
  //
  real result = 0.0;
  for (unsigned int i = 0; i < numBlobs; i++) {

    int partIdx = indices[i + offset];

    const real *p = &particles->positions[3 * partIdx];
    real radius = GetRadius(particles, partIdx); // @fixme { 0.5 is HACK }

    real3 center(p[0], p[1], p[2]);
    real3 delta = pos - center;
    result += radius /
              (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
  }

  return threshold - result;
}

real3 CalcBlobNormal(real3 point, const Particles *particles,
                     const std::vector<unsigned int> &indices,
                     unsigned int offset, int numBlobs, real threshold) {

  //
  // Normal of blob = Weighted avarage of normal for each sphere.
  //
  real q = 0.0; // sum of weight.
  real3 normal = real3(0.0, 0.0, 0.0);

  for (unsigned int i = 0; i < numBlobs; i++) {

    int partIdx = indices[i + offset];

    const real *p = &particles->positions[3 * partIdx];
    real radius = GetRadius(particles, partIdx);

    real3 center(p[0], p[1], p[2]);
    real3 delta = point - center;
    real r = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
    real weight = 0.0;
    if (radius < r) {
      weight = radius / r;
    }

    delta = vnormalize(delta);

    normal[0] += weight * delta[0];
    normal[1] += weight * delta[1];
    normal[2] += weight * delta[2];
  }

  return vnormalize(normal);
}

bool TestLeafNodeBlob(Intersection &isect, // [inout]
                      const ParticleNode &node,
                      const std::vector<unsigned int> &indices,
                      const Particles *particles, const Ray &ray, real minT,
                      real maxT) {

  // @todo
  const real threshold = 0.5;
  const int interval = 25;

  bool hit = false;

  unsigned int numBlobs = node.data[0];
  unsigned int offset = node.data[1];

  real3 rayorg;
  rayorg[0] = ray.origin()[0];
  rayorg[1] = ray.origin()[1];
  rayorg[2] = ray.origin()[2];

  real3 raydir;
  raydir[0] = ray.direction()[0];
  raydir[1] = ray.direction()[1];
  raydir[2] = ray.direction()[2];

  //
  // Very simple raymatch-blob intersection.
  //
  real step = (maxT - minT) / real(interval);
  real t = minT;

  real3 p;
  p[0] = rayorg[0] + t * raydir[0];
  p[1] = rayorg[1] + t * raydir[1];
  p[2] = rayorg[2] + t * raydir[2];

  real left, right;

  left = CalcBlobDensity(p, particles, indices, offset, numBlobs, threshold);
  for (int i = 0; i < interval; i++) {
    t += step;
    p[0] += step * raydir[0];
    p[1] += step * raydir[1];
    p[2] += step * raydir[2];
    right = CalcBlobDensity(p, particles, indices, offset, numBlobs, threshold);

    if (left * right < 0.0) {
      real hitT = t + right * step / (left - right);
      if (hitT < isect.t) {
        isect.t = hitT;
        isect.prim_id = indices[offset]; // @todo { Blend result of each blob. }

        // @todo { Delay computing normal for better performance. }
        real3 hitP = rayorg + real3(hitT, hitT, hitT) * raydir;
        real3 hitN = CalcBlobNormal(hitP, particles, indices, offset, numBlobs,
                                    threshold);
        isect.normal[0] = hitN[0];
        isect.normal[1] = hitN[1];
        isect.normal[2] = hitN[2];
        isect.geometric[0] = hitN[0];
        isect.geometric[1] = hitN[1];
        isect.geometric[2] = hitN[2];
        return true;
      }
    }

    left = right;
  }

  return false;
}
#endif

void BuildIntersection(Intersection &isect, const Particles *particles,
                       Ray &ray) {
  isect.position[0] = ray.origin()[0] + isect.t * ray.direction()[0];
  isect.position[1] = ray.origin()[1] + isect.t * ray.direction()[1];
  isect.position[2] = ray.origin()[2] + isect.t * ray.direction()[2];

  // Same ID for particle data.
  isect.f0 = isect.prim_id;
  isect.f1 = isect.prim_id;
  isect.f2 = isect.prim_id;

#if USE_BLOB
// normal is already computed. nothing to do.
#else

  real3 center;

  if (particles->isDoublePrecisionPos) {

    center[0] = particles->dpositions[3 * isect.prim_id + 0];
    center[1] = particles->dpositions[3 * isect.prim_id + 1];
    center[2] = particles->dpositions[3 * isect.prim_id + 2];

  } else {
    center[0] = particles->positions[3 * isect.prim_id + 0];
    center[1] = particles->positions[3 * isect.prim_id + 1];
    center[2] = particles->positions[3 * isect.prim_id + 2];
  }

  real3 cp;
  cp[0] = isect.position[0] - center[0];
  cp[1] = isect.position[1] - center[1];
  cp[2] = isect.position[2] - center[2];

  real3 n = vnormalize(cp);

  isect.normal[0] = n[0];
  isect.normal[1] = n[1];
  isect.normal[2] = n[2];

  isect.geometric[0] = n[0];
  isect.geometric[1] = n[1];
  isect.geometric[2] = n[2];

  // Compute tangent and binormal
  real tu = (atan2(n[0], n[2]) + M_PI) * 0.5 * (1.0 / M_PI);
  real tv = acos(n[1]) / M_PI;
  real xz2 = n[0] * n[0] + n[2] * n[2];
  if (xz2 > 0.0) {
    real xz = sqrt(xz2);
    real inv = 1.0 / xz;

    isect.tangent[0] = -2.0 * M_PI * n[0];
    isect.tangent[1] = 2.0 * M_PI * n[2];
    isect.tangent[2] = 0.0;

    isect.binormal[0] = -M_PI * n[2] * inv * n[1];
    isect.binormal[1] = -M_PI * n[0] * inv * n[1];
    isect.binormal[2] = M_PI * xz;
  } else {
    if (n[1] > 0.0) {
      isect.tangent[0] = 0.0;
      isect.tangent[1] = 0.0;
      isect.tangent[2] = 1.0;
      isect.binormal[0] = 1.0;
      isect.binormal[1] = 0.0;
      isect.binormal[2] = 0.0;
    } else {
      isect.tangent[0] = 0.0;
      isect.tangent[1] = 0.0;
      isect.tangent[2] = 1.0;
      isect.binormal[0] = -1.0;
      isect.binormal[1] = 0.0;
      isect.binormal[2] = 0.0;
    }
  }
  

#endif

  isect.u = tu;
  isect.v = tv;
}

} // namespace

bool ParticleAccel::Traverse(Intersection &isect, Ray &ray) const {
  real hitT = std::numeric_limits<real>::max(); // far = no hit.

  int nodeStackIndex = 0;
  std::vector<int> nodeStack(512);
  nodeStack[0] = 0;

  // Init isect info as no hit
  isect.t = hitT;
  isect.u = 0.0;
  isect.v = 0.0;
  isect.prim_id = (unsigned int)(-1);

  int dirSign[3];
  dirSign[0] = ray.direction()[0] < 0.0 ? 1 : 0;
  dirSign[1] = ray.direction()[1] < 0.0 ? 1 : 0;
  dirSign[2] = ray.direction()[2] < 0.0 ? 1 : 0;

  // @fixme { Check edge case; i.e., 1/0 }
  real3 rayInvDir;
  rayInvDir[0] = 1.0 / ray.direction()[0];
  rayInvDir[1] = 1.0 / ray.direction()[1];
  rayInvDir[2] = 1.0 / ray.direction()[2];

  real3 rayOrg;
  rayOrg[0] = ray.origin()[0];
  rayOrg[1] = ray.origin()[1];
  rayOrg[2] = ray.origin()[2];

#if ENABLE_SIMD_ISECTOR
  const double2 rox = _mm_set_pd(ray.origin()[0], ray.origin()[0]);
  const double2 roy = _mm_set_pd(ray.origin()[1], ray.origin()[1]);
  const double2 roz = _mm_set_pd(ray.origin()[2], ray.origin()[2]);

  const double2 ridx = _mm_set_pd(rayInvDir[0], rayInvDir[0]);
  const double2 ridy = _mm_set_pd(rayInvDir[1], rayInvDir[1]);
  const double2 ridz = _mm_set_pd(rayInvDir[2], rayInvDir[2]);

  const double2 roidx = _mm_mul_pd(rox, ridx);
  const double2 roidy = _mm_mul_pd(roy, ridy);
  const double2 roidz = _mm_mul_pd(roz, ridz);

  const double2 vtmin = _mm_set_pd(-hitT, -hitT);
  double2 vtmax = _mm_set_pd(hitT, hitT);

#endif

  real minT = -hitT, maxT = hitT;
  while (nodeStackIndex >= 0) {
    int index = nodeStack[nodeStackIndex];
    const ParticleNode &node = nodes_[index];

    nodeStackIndex--;

    if (node.flag == 0) { // branch node

#if ENABLE_SIMD_ISECTOR
      double2 bboxes[2][3];
      bboxes[0][0] = _mm_set_pd(node.bmin[0][0], node.bmin[1][0]);
      bboxes[0][1] = _mm_set_pd(node.bmin[0][1], node.bmin[1][1]);
      bboxes[0][2] = _mm_set_pd(node.bmin[0][2], node.bmin[1][2]);
      bboxes[1][0] = _mm_set_pd(node.bmax[0][0], node.bmax[1][0]);
      bboxes[1][1] = _mm_set_pd(node.bmax[0][1], node.bmax[1][1]);
      bboxes[1][2] = _mm_set_pd(node.bmax[0][2], node.bmax[1][2]);

#if defined(__sparc__) && defined(__HPC_ACE__)
      uint64_t hitmask[2];
      IntersectRayAABBD2(hitmask, bboxes, rox, roy, roz, ridx, ridy, ridz,
                         roidx, roidy, roidz, dirSign, vtmin, vtmax);

      if (hitmask[0] || hitmask[1]) {

        int orderNear = dirSign[node.axis];
        int orderFar = 1 - orderNear;

        // Traverse near node first.
        // Note: mask bit in hitmask was saved in swapped order. so use near
        // index for far node and far index for near node.
        if (hitmask[orderNear]) {
          nodeStack[++nodeStackIndex] = node.data[orderFar];
        }

        if (hitmask[orderFar]) {
          nodeStack[++nodeStackIndex] = node.data[orderNear];
        }
      }
#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

      int hitmask =
          IntersectRayAABBD2(bboxes, rox, roy, roz, ridx, ridy, ridz, roidx,
                             roidy, roidz, dirSign, vtmin, vtmax);

      if (hitmask) {

        int mask[2] = {(hitmask >> 1) & 0x1, (hitmask)&0x1};

        int orderNear = dirSign[node.axis];
        int orderFar = 1 - orderNear;

        // Traverse near node first.
        if (mask[orderFar]) {
          nodeStack[++nodeStackIndex] = node.data[orderFar];
        }

        if (mask[orderNear]) {
          nodeStack[++nodeStackIndex] = node.data[orderNear];
        }
      }

#else
#error SIMD is not supported on this architecture. Disable ENABLE_SIMD_ISECTOR to solve this error.
#endif

#else
      bool hits[2];
      hits[0] = IntersectRayAABB(minT, maxT, hitT, node.bmin[0], node.bmax[0],
                                 rayOrg, rayInvDir, dirSign);
      hits[1] = IntersectRayAABB(minT, maxT, hitT, node.bmin[1], node.bmax[1],
                                 rayOrg, rayInvDir, dirSign);

      if (hits[0] || hits[1]) {

        int orderNear = dirSign[node.axis];
        int orderFar = 1 - orderNear;

        // Traverse near first.
        if (hits[orderFar]) {
          nodeStack[++nodeStackIndex] = node.data[orderFar];
        }

        if (hits[orderNear]) {
          nodeStack[++nodeStackIndex] = node.data[orderNear];
        }
      }
#endif

    } else { // leaf node

#if USE_BLOB
      if (TestLeafNodeBlob(isect, node, indices_, particles_, ray, minT,
                           maxT)) {
        hitT = isect.t;
      }
#else
      if (TestLeafNode(isect, node, indices_, particles_, ray)) {
        hitT = isect.t;

        // Shorten t max
        maxT = hitT;
#if ENABLE_SIMD_ISECTOR
        vtmax = _mm_set_pd(hitT, hitT);
#endif
      }
#endif
    }
  }

  assert(nodeStackIndex < kMaxStackDepth);

  if (isect.t < std::numeric_limits<real>::max()) {
    BuildIntersection(isect, particles_, ray);
    return true;
  }

  return false;
}

void ParticleAccel::BoundingBox(double bmin[3], double bmax[3]) {
  if (nodes_.empty()) {
    bmin[0] = std::numeric_limits<double>::max();
    bmin[1] = std::numeric_limits<double>::max();
    bmin[2] = std::numeric_limits<double>::max();
    bmax[0] = -std::numeric_limits<double>::max();
    bmax[1] = -std::numeric_limits<double>::max();
    bmax[2] = -std::numeric_limits<double>::max();
  } else {
    bmin[0] = bmin_[0];
    bmin[1] = bmin_[1];
    bmin[2] = bmin_[2];
    bmax[0] = bmax_[0];
    bmax[1] = bmax_[1];
    bmax[2] = bmax_[2];
  }
}
