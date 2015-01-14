/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <limits>
#include <functional>
#include <algorithm>

#include "accel_triangle.h"
#include "prefix_tree_util.h"

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10
#include <emmintrin.h>
#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2
#include <emmintrin.h>
#endif

using namespace lsgl::render;

#define MAX_LEAF_ELEMENTS (16)
#define MAX_TREE_DEPTH_32BIT                                                   \
  (22) // FYI, log2(1G/16) ~= 25.897, log2(1G/32) ~= 21

#define ENABLE_TRACE_PRINT (0)
#define ENABLE_DEBUG_PRINT (0)

#define ENABLE_SIMD_ISECTOR (0)

#define trace(f, ...)                                                          \
  {                                                                            \
    if (ENABLE_TRACE_PRINT)                                                    \
      printf(f, __VA_ARGS__);                                                  \
  }

#if ENABLE_DEBUG_PRINT
#define debug(f, ...)                                                          \
  { printf(f, __VA_ARGS__); }
#else
#define debug(f, ...)
#endif

namespace {

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10

typedef __m128d double2;

//#define FORCEINLINE __attribute__((always_inline))

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
  const double kEPS = std::numeric_limits<real>::epsilon() * 1024;

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

#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

typedef __m128d double2;

//#define FORCEINLINE __attribute__((always_inline))

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
  const double kEPS = std::numeric_limits<real>::epsilon() * 1024;

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

  // if (t < 0.0 || t > tInOut)
  //  return false;
  double2 tMask = _mm_and_pd(_mm_cmpge_pd(t, vzero), _mm_cmple_pd(t, tInOut));

  // if (u < 0.0 || u > 1.0)
  //  return false;
  double2 uMask = _mm_and_pd(_mm_cmpge_pd(u, vzero), _mm_cmple_pd(u, vone));

  // if (v < 0.0 || u + v > 1.0)
  //  return false;
  double2 uvMask =
      _mm_and_pd(_mm_cmpge_pd(v, vzero), _mm_cmple_pd(_mm_add_pd(u, v), vone));

  double2 finalMask =
      _mm_and_pd(_mm_and_pd(detMask, tMask), _mm_and_pd(uMask, uvMask));

  {
    // early exit.
    if (_mm_movemask_pd(finalMask) == 0) {
      return;
    }
  }

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

//
// SAH functions
//

struct BinBuffer {

  BinBuffer(int size) {
    binSize = size;
    bin.resize(2 * 3 * size);
    clear();
  }

  void clear() { memset(&bin[0], 0, sizeof(size_t) * 2 * 3 * binSize); }

  std::vector<size_t> bin; // (min, max) * xyz * binsize
  int binSize;
};

static inline double CalculateSurfaceArea(const real3 &min, const real3 &max) {
  real3 box = max - min;
  return 2.0 * (box[0] * box[1] + box[1] * box[2] + box[2] * box[0]);
}

static inline void GetBoundingBoxOfTriangle(real3 &bmin, real3 &bmax,
                                            const Mesh *mesh,
                                            unsigned int index) {
  unsigned int f0 = mesh->faces[3 * index + 0];
  unsigned int f1 = mesh->faces[3 * index + 1];
  unsigned int f2 = mesh->faces[3 * index + 2];

  real3 p[3];

  if (mesh->isDoublePrecisionPos) {
    p[0] = real3(mesh->dvertices[3 * f0 + 0], mesh->dvertices[3 * f0 + 1],
                 mesh->dvertices[3 * f0 + 2]);
    p[1] = real3(mesh->dvertices[3 * f1 + 0], mesh->dvertices[3 * f1 + 1],
                 mesh->dvertices[3 * f1 + 2]);
    p[2] = real3(mesh->dvertices[3 * f2 + 0], mesh->dvertices[3 * f2 + 1],
                 mesh->dvertices[3 * f2 + 2]);
  } else {
    p[0] = real3(mesh->vertices[3 * f0 + 0], mesh->vertices[3 * f0 + 1],
                 mesh->vertices[3 * f0 + 2]);
    p[1] = real3(mesh->vertices[3 * f1 + 0], mesh->vertices[3 * f1 + 1],
                 mesh->vertices[3 * f1 + 2]);
    p[2] = real3(mesh->vertices[3 * f2 + 0], mesh->vertices[3 * f2 + 1],
                 mesh->vertices[3 * f2 + 2]);
  }

  bmin = p[0];
  bmax = p[0];

  for (int i = 1; i < 3; i++) {
    bmin[0] = std::min(bmin[0], p[i][0]);
    bmin[1] = std::min(bmin[1], p[i][1]);
    bmin[2] = std::min(bmin[2], p[i][2]);

    bmax[0] = std::max(bmax[0], p[i][0]);
    bmax[1] = std::max(bmax[1], p[i][1]);
    bmax[2] = std::max(bmax[2], p[i][2]);
  }
}

static void ContributeBinBuffer(BinBuffer *bins, // [out]
                                const real3 &sceneMin, const real3 &sceneMax,
                                const Mesh *mesh, unsigned int *indices,
                                unsigned int leftIdx, unsigned int rightIdx) {
  static const real EPS = std::numeric_limits<real>::epsilon() * 1024;

  real binSize = (real)bins->binSize;

  // Calculate extent
  real3 sceneSize, sceneInvSize;
  sceneSize = sceneMax - sceneMin;
  for (int i = 0; i < 3; ++i) {
    assert(sceneSize[i] >= 0.0);

    if (sceneSize[i] > EPS) {
      sceneInvSize[i] = binSize / sceneSize[i];
    } else {
      sceneInvSize[i] = 0.0;
    }
  }

  // Clear bin data
  memset(&bins->bin[0], 0, sizeof(2 * 3 * bins->binSize));

  size_t idxBMin[3];
  size_t idxBMax[3];

  for (size_t i = leftIdx; i < rightIdx; i++) {

    //
    // Quantize the position into [0, BIN_SIZE)
    //
    // q[i] = (int)(p[i] - scene_bmin) / scene_size
    //
    real3 bmin;
    real3 bmax;

    GetBoundingBoxOfTriangle(bmin, bmax, mesh, indices[i]);

    real3 quantizedBMin = (bmin - sceneMin) * sceneInvSize;
    real3 quantizedBMax = (bmax - sceneMin) * sceneInvSize;

    // idx is now in [0, BIN_SIZE)
    for (size_t j = 0; j < 3; ++j) {
      idxBMin[j] = (unsigned int)floor(quantizedBMin[j]);
      idxBMax[j] = (unsigned int)floor(quantizedBMax[j]);

      if (idxBMin[j] >= binSize)
        idxBMin[j] = binSize - 1;
      if (idxBMax[j] >= binSize)
        idxBMax[j] = binSize - 1;

      assert(idxBMin[j] < binSize);
      assert(idxBMax[j] < binSize);

      // Increment bin counter
      bins->bin[0 * (bins->binSize * 3) + j * bins->binSize + idxBMin[j]] += 1;
      bins->bin[1 * (bins->binSize * 3) + j * bins->binSize + idxBMax[j]] += 1;
    }
  }
}

static inline double SAH(size_t ns1, real leftArea, size_t ns2, real rightArea,
                         real invS, real Taabb, real Ttri) {
  // const real Taabb = 0.2f;
  // const real Ttri = 0.8f;
  real T;

  T = 2.0f * Taabb + (leftArea * invS) * (real)(ns1) * Ttri +
      (rightArea * invS) * (real)(ns2) * Ttri;

  return T;
}

static bool FindCutFromBinBuffer(real *cutPos,     // [out] xyz
                                 int &minCostAxis, // [out]
                                 const BinBuffer *bins, const real3 &bmin,
                                 const real3 &bmax, size_t numTriangles,
                                 real costTaabb) // should be in [0.0, 1.0]
{
  const real eps = std::numeric_limits<real>::epsilon() * 1024;

  size_t left, right;
  real3 bsize, bstep;
  real3 bminLeft, bmaxLeft;
  real3 bminRight, bmaxRight;
  real saLeft, saRight, saTotal;
  real pos;
  real minCost[3];

  real costTtri = 1.0 - costTaabb;

  minCostAxis = 0;

  bsize = bmax - bmin;
  bstep = bsize * (1.0 / bins->binSize);
  saTotal = CalculateSurfaceArea(bmin, bmax);

  real invSaTotal = 0.0;
  if (saTotal > eps) {
    invSaTotal = 1.0 / saTotal;
  }

  for (int j = 0; j < 3; ++j) {

    //
    // Compute SAH cost for right side of each cell of the bbox.
    // Exclude both extreme side of the bbox.
    //
    //  i:      0    1    2    3
    //     +----+----+----+----+----+
    //     |    |    |    |    |    |
    //     +----+----+----+----+----+
    //

    real minCostPos = bmin[j] + 0.5 * bstep[j];
    minCost[j] = std::numeric_limits<real>::max();

    left = 0;
    right = numTriangles;
    bminLeft = bminRight = bmin;
    bmaxLeft = bmaxRight = bmax;

    for (int i = 0; i < bins->binSize - 1; ++i) {
      left += bins->bin[0 * (3 * bins->binSize) + j * bins->binSize + i];
      right -= bins->bin[1 * (3 * bins->binSize) + j * bins->binSize + i];

      assert(left <= numTriangles);
      assert(right <= numTriangles);

      //
      // Split pos bmin + (i + 1) * (bsize / BIN_SIZE)
      // +1 for i since we want a position on right side of the cell.
      //

      pos = bmin[j] + (i + 0.5) * bstep[j];
      bmaxLeft[j] = pos;
      bminRight[j] = pos;

      saLeft = CalculateSurfaceArea(bminLeft, bmaxLeft);
      saRight = CalculateSurfaceArea(bminRight, bmaxRight);

      real cost =
          SAH(left, saLeft, right, saRight, invSaTotal, costTaabb, costTtri);
      if (cost < minCost[j]) {
        //
        // Update the min cost
        //
        minCost[j] = cost;
        minCostPos = pos;
        // minCostAxis = j;
      }
    }

    cutPos[j] = minCostPos;
  }

  // cutAxis = minCostAxis;
  // cutPos = minCostPos;

  // Find min cost axis
  real cost = minCost[0];
  minCostAxis = 0;
  if (cost > minCost[1]) {
    minCostAxis = 1;
    cost = minCost[1];
  }
  if (cost > minCost[2]) {
    minCostAxis = 2;
    cost = minCost[2];
  }

  return true;
}

class SAHPred : public std::unary_function<unsigned int, bool> {
public:
  SAHPred(int axis, real pos, const Mesh *mesh)
      : axis_(axis), pos_(pos), mesh_(mesh) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;
    real pos = pos_;

    unsigned int i0 = mesh_->faces[3 * i + 0];
    unsigned int i1 = mesh_->faces[3 * i + 1];
    unsigned int i2 = mesh_->faces[3 * i + 2];

    real3 p0(mesh_->vertices[3 * i0 + 0], mesh_->vertices[3 * i0 + 1],
             mesh_->vertices[3 * i0 + 2]);
    real3 p1(mesh_->vertices[3 * i1 + 0], mesh_->vertices[3 * i1 + 1],
             mesh_->vertices[3 * i1 + 2]);
    real3 p2(mesh_->vertices[3 * i2 + 0], mesh_->vertices[3 * i2 + 1],
             mesh_->vertices[3 * i2 + 2]);

    real center = p0[axis] + p1[axis] + p2[axis];

    return (center < pos * 3.0);
  }

private:
  int axis_;
  real pos_;
  const Mesh *mesh_;
};

class SAHPredD : public std::unary_function<unsigned int, bool> {
public:
  SAHPredD(int axis, double pos, const Mesh *mesh)
      : axis_(axis), pos_(pos), mesh_(mesh) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;
    double pos = pos_;

    unsigned int i0 = mesh_->faces[3 * i + 0];
    unsigned int i1 = mesh_->faces[3 * i + 1];
    unsigned int i2 = mesh_->faces[3 * i + 2];

    double3 p0(mesh_->dvertices[3 * i0 + 0], mesh_->dvertices[3 * i0 + 1],
               mesh_->dvertices[3 * i0 + 2]);
    double3 p1(mesh_->dvertices[3 * i1 + 0], mesh_->dvertices[3 * i1 + 1],
               mesh_->dvertices[3 * i1 + 2]);
    double3 p2(mesh_->dvertices[3 * i2 + 0], mesh_->dvertices[3 * i2 + 1],
               mesh_->dvertices[3 * i2 + 2]);

    double center = p0[axis] + p1[axis] + p2[axis];

    return (center < pos * 3.0);
  }

private:
  int axis_;
  double pos_;
  const Mesh *mesh_;
};

template <typename T>
static void ComputeBoundingBox(real3 &bmin, real3 &bmax, const T *vertices,
                               unsigned int *faces, unsigned int *indices,
                               unsigned int leftIndex,
                               unsigned int rightIndex) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 1024;

  size_t i = leftIndex;
  size_t idx = indices[i];
  bmin[0] = vertices[3 * faces[3 * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * faces[3 * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * faces[3 * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * faces[3 * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * faces[3 * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * faces[3 * idx + 0] + 2] + kEPS;

  // Assume mesh are composed of all triangles
  for (i = leftIndex; i < rightIndex; i++) { // for each faces
    size_t idx = indices[i];
    for (int j = 0; j < 3; j++) { // for each face vertex
      size_t fid = faces[3 * idx + j];
      for (int k = 0; k < 3; k++) { // xyz
        real minval = vertices[3 * fid + k] - kEPS;
        real maxval = vertices[3 * fid + k] + kEPS;
        if (bmin[k] > minval)
          bmin[k] = minval;
        if (bmax[k] < maxval)
          bmax[k] = maxval;
      }
    }
  }
}

template <typename T>
static void
ComputeBoundingBox30(real3 &bmin, real3 &bmax, const T *vertices,
                     const uint32_t *faces, const std::vector<IndexKey30> &keys,
                     unsigned int leftIndex, unsigned int rightIndex) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 1024;

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

  size_t i = leftIndex;
  size_t idx = keys[i].index;
  bmin[0] = vertices[3 * faces[3 * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * faces[3 * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * faces[3 * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * faces[3 * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * faces[3 * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * faces[3 * idx + 0] + 2] + kEPS;

  // Assume mesh are composed of all triangles
  for (i = leftIndex; i < rightIndex; i++) { // for each faces
    size_t idx = keys[i].index;
    for (int j = 0; j < 3; j++) { // for each face vertex
      size_t fid = faces[3 * idx + j];
      for (int k = 0; k < 3; k++) { // xyz
        real minval = vertices[3 * fid + k] - kEPS;
        real maxval = vertices[3 * fid + k] + kEPS;
        if (bmin[k] > minval)
          bmin[k] = minval;
        if (bmax[k] < maxval)
          bmax[k] = maxval;
      }
    }
  }
}

#ifdef _OPENMP
template <typename T>
void ComputeBoundingBoxOMP(real3 &bmin, real3 &bmax, const T *vertices,
                           const uint32_t *faces, unsigned int *indices,
                           unsigned int leftIndex, unsigned int rightIndex) {

  // assert(leftIndex < rightIndex);
  // assert(rightIndex - leftIndex > 0);

  const real kEPS = std::numeric_limits<real>::epsilon() * 1024;

  bmin[0] = std::numeric_limits<real>::max();
  bmin[1] = std::numeric_limits<real>::max();
  bmin[2] = std::numeric_limits<real>::max();
  bmax[0] = -std::numeric_limits<real>::max();
  bmax[1] = -std::numeric_limits<real>::max();
  bmax[2] = -std::numeric_limits<real>::max();

  {
    size_t i = leftIndex;
    size_t idx = indices[i];
    bmin[0] = vertices[3 * faces[3 * idx + 0] + 0] - kEPS;
    bmin[1] = vertices[3 * faces[3 * idx + 0] + 1] - kEPS;
    bmin[2] = vertices[3 * faces[3 * idx + 0] + 2] - kEPS;
    bmax[0] = vertices[3 * faces[3 * idx + 0] + 0] + kEPS;
    bmax[1] = vertices[3 * faces[3 * idx + 0] + 1] + kEPS;
    bmax[2] = vertices[3 * faces[3 * idx + 0] + 2] + kEPS;
  }

  // We could use min and max reduction if OpenMP 3.1 ready compiler was
  // available.

  real local_bmin[3] = {bmin[0], bmin[1], bmin[2]};
  real local_bmax[3] = {bmax[0], bmax[1], bmax[2]};

  size_t n = rightIndex - leftIndex;

#pragma omp parallel firstprivate(local_bmin, local_bmax) if (n > (1024 * 124))
  {

#pragma omp for
    for (size_t i = leftIndex; i < rightIndex; i++) {

      size_t idx = indices[i];

      for (int k = 0; k < 3; k++) {
        real minval_x = vertices[3 * faces[3 * idx + k] + 0] - kEPS;
        real minval_y = vertices[3 * faces[3 * idx + k] + 1] - kEPS;
        real minval_z = vertices[3 * faces[3 * idx + k] + 2] - kEPS;

        real maxval_x = vertices[3 * faces[3 * idx + k] + 0] + kEPS;
        real maxval_y = vertices[3 * faces[3 * idx + k] + 1] + kEPS;
        real maxval_z = vertices[3 * faces[3 * idx + k] + 2] + kEPS;

        local_bmin[0] = std::min(local_bmin[0], minval_x);
        local_bmin[1] = std::min(local_bmin[1], minval_y);
        local_bmin[2] = std::min(local_bmin[2], minval_z);

        local_bmax[0] = std::max(local_bmax[0], maxval_x);
        local_bmax[1] = std::max(local_bmax[1], maxval_y);
        local_bmax[2] = std::max(local_bmax[2], maxval_z);
      }
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

template <typename T>
void MakeLeaf32(BVHNode &leaf, const T *vertices, const uint32_t *faces,
                real3 &bmin, real3 &bmax, const std::vector<IndexKey30> &keys,
                uint32_t leftIndex, uint32_t rightIndex) {

  // 1. Compute leaf AABB
  ComputeBoundingBox30(bmin, bmax, vertices, faces, keys, leftIndex,
                       rightIndex);

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
template <typename T>
size_t BuildTreeRecursive32(std::vector<BVHNode> &nodes, real3 &bmin,
                            real3 &bmax, const std::vector<IndexKey30> &keys,
                            const std::vector<NodeInfo32> &nodeInfos,
                            const T *vertices, const uint32_t *faces,
                            uint32_t rootIndex, uint32_t leftIndex,
                            uint32_t rightIndex, bool isLeaf, int depth) {
  InvalidateBoundingBox(bmin, bmax);

  uint32_t n = rightIndex - leftIndex;

  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= MAX_LEAF_ELEMENTS) || (depth > MAX_TREE_DEPTH_32BIT)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;
    // if (leftIndex == rightIndex) { // this would be OK. 1 tri in 1 leaf case.
    //  endIndex++;
    //}

    BVHNode leaf;
    MakeLeaf32(leaf, vertices, faces, bmin, bmax, keys, leftIndex, endIndex);

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

  BVHNode node;
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
      nodes, leftBMin, leftBMax, keys, nodeInfos, vertices, faces, midIndex,
      leftIndex, midIndex, isLeftLeaf, depth + 1);
  size_t rightChildIndex = BuildTreeRecursive32(
      nodes, rightBMin, rightBMax, keys, nodeInfos, vertices, faces,
      midIndex + 1, midIndex + 1, rightIndex, isRightLeaf, depth + 1);

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

#ifdef _OPENMP
template <typename T>
size_t BuildTreeRecursive32OMP(std::vector<BVHNode> &nodes, real3 &bmin,
                               real3 &bmax, const std::vector<IndexKey30> &keys,
                               const std::vector<NodeInfo32> &nodeInfos,
                               const T *vertices, const uint32_t *faces,
                               uint32_t rootIndex, uint32_t leftIndex,
                               uint32_t rightIndex, bool isLeaf, int depth) {
  InvalidateBoundingBox(bmin, bmax);

  uint32_t n = rightIndex - leftIndex;

  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= MAX_LEAF_ELEMENTS) || (depth > MAX_TREE_DEPTH_32BIT)) {
    // if (isLeaf || (n <= 0)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;

    BVHNode leaf;
    MakeLeaf32(leaf, vertices, faces, bmin, bmax, keys, leftIndex, endIndex);

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

  BVHNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  size_t offset;

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

  if (depth > 6) {

    // Enough number of tasks was launched. Switch to sequential code.
    std::vector<BVHNode> sub_nodes;

    leftChildIndex = BuildTreeRecursive32(
        sub_nodes, leftBMin, leftBMax, keys, nodeInfos, vertices, faces,
        midIndex, leftIndex, midIndex, isLeftLeaf, depth + 1);

    rightChildIndex = BuildTreeRecursive32(
        sub_nodes, rightBMin, rightBMax, keys, nodeInfos, vertices, faces,
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
                        leftBMax, keys, nodeInfos, vertices,                   \
                        faces) firstprivate(midIndex) if (depth < 10)
    leftChildIndex = BuildTreeRecursive32OMP(
        nodes, leftBMin, leftBMax, keys, nodeInfos, vertices, faces, midIndex,
        leftIndex, midIndex, isLeftLeaf, depth + 1);

#pragma omp task shared(leftIndex, rightChildIndex, nodes, rightBMin,          \
                        rightBMax, keys, nodeInfos, vertices,                  \
                        faces) firstprivate(midIndex) if (depth < 10)
    rightChildIndex = BuildTreeRecursive32OMP(
        nodes, rightBMin, rightBMax, keys, nodeInfos, vertices, faces,
        midIndex + 1, midIndex + 1, rightIndex, isRightLeaf, depth + 1);

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

size_t BVHAccel::BuildTree(const Mesh *mesh, real3 &bmin, real3 &bmax,
                           unsigned int leftIdx, unsigned int rightIdx,
                           int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  size_t offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
    stats_.maxTreeDepth = depth;
  }

  if (mesh->isDoublePrecisionPos) {
    ComputeBoundingBox<double>(bmin, bmax, mesh->dvertices, mesh->faces,
                               &indices_.at(0), leftIdx, rightIdx);
  } else {
    ComputeBoundingBox<float>(bmin, bmax, mesh->vertices, mesh->faces,
                              &indices_.at(0), leftIdx, rightIdx);
  }

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {
    // Create leaf node.
    BVHNode leaf;

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
  // Compute SAH and find best split axis and position
  //
  int minCutAxis = 0;
  real cutPos[3] = {0.0, 0.0, 0.0};

  BinBuffer bins(options_.binSize);
  ContributeBinBuffer(&bins, bmin, bmax, mesh, &indices_.at(0), leftIdx,
                      rightIdx);
  FindCutFromBinBuffer(cutPos, minCutAxis, &bins, bmin, bmax, n,
                       options_.costTaabb);

  debug("depth: %d, cutPos: (%f, %f, %f), cutAxis: %d\n", depth, cutPos[0],
        cutPos[1], cutPos[2], minCutAxis);

  // Try all 3 axis until good cut position avaiable.
  unsigned int midIdx;
  int cutAxis = minCutAxis;
  for (int axisTry = 0; axisTry < 1; axisTry++) {

    unsigned int *begin = &indices_[leftIdx];
    unsigned int *end = &indices_[rightIdx];
    unsigned int *mid = 0;

    // try minCutAxis first.
    cutAxis = (minCutAxis + axisTry) % 3;

    //
    // Split at (cutAxis, cutPos)
    // indices_ will be modified.
    //
    if (mesh->isDoublePrecisionPos) {
      mid =
          std::partition(begin, end, SAHPredD(cutAxis, cutPos[cutAxis], mesh));
    } else {
      mid = std::partition(begin, end, SAHPred(cutAxis, cutPos[cutAxis], mesh));
    }

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

  BVHNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch
  nodes_.push_back(node);

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;
  unsigned int leftChildIndex =
      BuildTree(mesh, leftBMin, leftBMax, leftIdx, midIdx, depth + 1);
  unsigned int rightChildIndex =
      BuildTree(mesh, rightBMin, rightBMax, midIdx, rightIdx, depth + 1);

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

bool BVHAccel::Build(const Mesh *mesh, const BVHBuildOptions &options) {
  options_ = options;
  stats_ = BVHBuildStatistics();

  assert(options_.binSize > 1);

  assert(mesh);

  size_t n = mesh->nfaces;
  trace("[BVHAccel] Input # of vertices = %lu\n", mesh->nvertices);
  trace("[BVHAccel] Input # of faces    = %lu\n", mesh->nfaces);

  //
  // 1. Create triangle indices(this will be permutated in BuildTree)
  //
  indices_.resize(n);
  for (size_t i = 0; i < n; i++) {
    indices_[i] = i;
  }

  //
  // 2. Build tree
  //
  real3 rootBMin, rootBMax;
  BuildTree(mesh, rootBMin, rootBMax, 0, n, 0);

  // Tree will be null if input triangle count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    trace("[BVHAccel] bound min = (%f, %f, %f)\n", rootBMin[0], rootBMin[1],
          rootBMin[2]);
    trace("[BVHAccel] bound max = (%f, %f, %f)\n", rootBMax[0], rootBMax[1],
          rootBMax[2]);
  }

  trace("[BVHAccel] # of nodes = %lu\n", nodes_.size());

  // Store root bbox.
  bmin_[0] = rootBMin[0];
  bmin_[1] = rootBMin[1];
  bmin_[2] = rootBMin[2];

  bmax_[0] = rootBMax[0];
  bmax_[1] = rootBMax[1];
  bmax_[2] = rootBMax[2];

  // Store pointer for later use.
  mesh_ = mesh;

  return true;
}

bool BVHAccel::Build32(const Mesh *mesh, const BVHBuildOptions &options) {

  assert(mesh);

  if (mesh->nfaces < 1024) {
    // Use non-optimized BVH builder.
    return Build(mesh, options);
  }

  options_ = options;
  stats_ = BVHBuildStatistics();

  assert(options_.binSize > 1);
  assert(mesh->isDoublePrecisionPos == false); // @todo

  assert(mesh);

  size_t n = mesh->nfaces;

  trace("[BVHAccel2] Input # of vertices = %lu\n", mesh->nvertices);
  trace("[BVHAccel2] Input # of faces    = %lu\n", mesh->nfaces);

  real3 rootBMin, rootBMax;

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

  {
    timerutil t;
    t.start();

    real3 bmin, bmax;

#ifdef _OPENMP
    ComputeBoundingBoxOMP(bmin, bmax, mesh->vertices, mesh->faces,
                          &indices_.at(0), 0, n);
#else
    ComputeBoundingBox(bmin, bmax, mesh->vertices, mesh->faces, &indices_.at(0),
                       0, n);
#endif

    t.end();
    printf("[2:scene bbox calculation] %d msec\n", (int)t.msec());

    printf(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    printf(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

    std::vector<uint32_t> codes(n);
    assert(sizeof(real) == sizeof(float));

    {
      timerutil t;
      t.start();

      if (mesh->isDoublePrecisionPos) {
        assert(0); // @todo
      } else {
        CalculateMortonCodesTriangleFloat30(&codes.at(0), mesh->vertices,
                                            mesh->faces, bmin, bmax, 0, n);
      }
      t.end();
      printf("[3:morton calculation] %d msec\n", (int)t.msec());

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

    { // sort
      timerutil t;
      t.start();
      std::vector<IndexKey30> temp(n);

#ifdef _OPENMP
#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10
      // OpenMP version doesn't work well. Use sequential method for now.
      RadixSort30(&keys.at(0), &keys.at(0) + n);
#else
      RadixSort30OMP(&keys.at(0), &keys.at(0) + n);
#pragma omp barrier
#endif
      t.end();
      printf("[4:radix sort] %d msec\n", (int)t.msec());

#else
      RadixSort30(&keys.at(0), &keys.at(0) + n);
#endif

      // Check
      // for (size_t i = 0; i < n-1; i++) {
      // std::string b = BitString32(keys[i].code);
      // printf("[%08d] i = %010d, c = %010d(%s)\n", i, keys[i].index,
      // keys[i].code, b.c_str());
      // assert(keys[i].code <= keys[i+1].code);
      //}
    }

    std::vector<NodeInfo32> nodeInfos(n - 1);
    {
      timerutil t;
      t.start();

#ifdef _OPENMP
//#pragma omp parallel forif (n > (1024 * 1024))
#pragma omp parallel for schedule(dynamic)
#endif
      for (size_t i = 0; i < n - 1; i++) {
        nodeInfos[i] = ConstructBinaryRadixTree30(&keys.at(0), i, n);
        // printf("I[%d].index   = %d\n", i, nodeInfos[i].index);
        // printf("I[%d].leftTy  = %d\n", i, nodeInfos[i].leftType);
        // printf("I[%d].rightTy = %d\n", i, nodeInfos[i].rightType);
      }

      t.end();
      printf("[5:Construct binary radix tree: %d msec\n", (int)t.msec());
    }

    {
      timerutil t;
      t.start();

      nodes_.clear();

      // Explicitly create root node and reserve storage here.
      BVHNode rootNode;
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

#ifdef _OPENMP

      size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

#pragma omp parallel shared(leftChildIndex, rightChildIndex, keys, nodeInfos)
      {
#pragma omp single
        {
          leftChildIndex = BuildTreeRecursive32OMP(
              nodes_, leftBMin, leftBMax, keys, nodeInfos, mesh->vertices,
              mesh->faces, midIndex, 0, midIndex, isLeftLeaf, 0);
        }

#pragma omp single
        {
          rightChildIndex = BuildTreeRecursive32OMP(
              nodes_, rightBMin, rightBMax, keys, nodeInfos, mesh->vertices,
              mesh->faces, midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);
        }
      }
#pragma omp barrier

      assert(leftChildIndex != (size_t)(-1));
      assert(rightChildIndex != (size_t)(-1));

#else

      size_t leftChildIndex = BuildTreeRecursive32(
          nodes_, leftBMin, leftBMax, keys, nodeInfos, mesh->vertices,
          mesh->faces, midIndex, 0, midIndex, isLeftLeaf, 0);
      size_t rightChildIndex = BuildTreeRecursive32(
          nodes_, rightBMin, rightBMax, keys, nodeInfos, mesh->vertices,
          mesh->faces, midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);

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
      printf("[6:Construct final AABB tree: %d msec\n", (int)t.msec());

      printf("  bmin = %f, %f, %f\n", rootBMin[0], rootBMin[1], rootBMin[2]);
      printf("  bmax = %f, %f, %f\n", rootBMax[0], rootBMax[1], rootBMax[2]);
    }

    {
      // Store sorted indices.
      assert(indices_.size() == keys.size());
      for (size_t i = 0; i < keys.size(); i++) {
        indices_[i] = keys[i].index;
      }
    }
  }

  // Tree will be null if input triangle count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    trace("[BVHAccel] bound min = (%f, %f, %f)\n", rootBMin[0], rootBMin[1],
          rootBMin[2]);
    trace("[BVHAccel] bound max = (%f, %f, %f)\n", rootBMax[0], rootBMax[1],
          rootBMax[2]);
  }

  trace("[BVHAccel] # of nodes = %lu\n", nodes_.size());

  // Store bbox
  bmin_[0] = rootBMin[0];
  bmin_[1] = rootBMin[1];
  bmin_[2] = rootBMin[2];

  bmax_[0] = rootBMax[0];
  bmax_[1] = rootBMax[1];
  bmax_[2] = rootBMax[2];

  // Store pointer for later use.
  mesh_ = mesh;

  return true;
}

bool BVHAccel::Dump(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "[BVHAccel] Cannot write a file: %s\n", filename);
    return false;
  }

  unsigned long long numNodes = nodes_.size();
  assert(nodes_.size() > 0);

  unsigned long long numIndices = indices_.size();

  int r = 0;
  r = fwrite(&numNodes, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&nodes_.at(0), sizeof(BVHNode), numNodes, fp);
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

  double2 vret =
      _mm_and_pd(_mm_cmpge_pd(tmax, vzerod2()), _mm_cmpge_pd(tmax, tmin));

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

    // printf("tmin, tmax = %f, %f\n", tmin, tmax);
    tminOut = tmin;
    tmaxOut = tmax;

    return true;
  }

  return false; // no hit
}

inline real3 vcross(real3 a, real3 b) {
  real3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

inline real vdot(real3 a, real3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double3 vcrossd(double3 a, double3 b) {
  double3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

inline double vdotd(double3 a, double3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline bool TriangleIsect(real &tInOut, real &uOut, real &vOut,
                          const vector3 &v0, const vector3 &v1,
                          const vector3 &v2, const real3 &rayOrg,
                          const real3 &rayDir, bool doubleSided) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 1024;

  real3 p0(v0[0], v0[1], v0[2]);
  real3 p1(v1[0], v1[1], v1[2]);
  real3 p2(v2[0], v2[1], v2[2]);
  real3 e1, e2;
  real3 p, s, q;

  e1 = p1 - p0;
  e2 = p2 - p0;

  p = vcross(rayDir, e2);

  real invDet;
  real det = vdot(e1, p);

  if (doubleSided) {
    if (std::abs(det) < kEPS) {
      return false;
    }
  } else {
    if (det < kEPS) { // single-sided
      return false;
    }
  }

  invDet = 1.0 / det;

  s = rayOrg - p0;
  q = vcross(s, e1);

  real u = vdot(s, p) * invDet;
  real v = vdot(q, rayDir) * invDet;
  real t = vdot(e2, q) * invDet;

  if (u < 0.0 || u > 1.0)
    return false;
  if (v < 0.0 || u + v > 1.0)
    return false;
  if (t < 0.0 || t > tInOut)
    return false;

  tInOut = t;
  uOut = u;
  vOut = v;

  return true;
}

inline bool TriangleIsectD(real &tInOut, real &uOut, real &vOut,
                           const vector3d &v0, const vector3d &v1,
                           const vector3d &v2, const double3 &rayOrg,
                           const double3 &rayDir, bool doubleSided) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 1024;

  double3 p0(v0[0], v0[1], v0[2]);
  double3 p1(v1[0], v1[1], v1[2]);
  double3 p2(v2[0], v2[1], v2[2]);
  double3 e1, e2;
  double3 p, s, q;

  e1 = p1 - p0;
  e2 = p2 - p0;

  p = vcrossd(rayDir, e2);

  double invDet;
  double det = vdotd(e1, p);

  if (doubleSided) {
    if (std::abs(det) < kEPS) {
      return false;
    }
  } else {
    if (det < kEPS) { // single-sided
      return false;
    }
  }

  invDet = 1.0 / det;

  s = rayOrg - p0;
  q = vcrossd(s, e1);

  double u = vdotd(s, p) * invDet;
  double v = vdotd(q, rayDir) * invDet;
  double t = vdotd(e2, q) * invDet;

  if (u < 0.0 || u > 1.0)
    return false;
  if (v < 0.0 || u + v > 1.0)
    return false;
  if (t < 0.0 || t > tInOut)
    return false;

  tInOut = t;
  uOut = u;
  vOut = v;

  return true;
}

static bool FORCEINLINE TestLeafNode(Intersection &isect, // [inout]
                                     const BVHNode &node,
                                     const std::vector<unsigned int> &indices,
                                     const Mesh *mesh, const Ray &ray) {

  bool hit = false;

  unsigned int numTriangles = node.data[0];
  unsigned int offset = node.data[1];

  real t = isect.t; // current hit distance

  bool doubleSided = (bool)ray.double_sided;

  if (mesh->isDoublePrecisionPos) {

    double3 rayOrg;
    rayOrg[0] = ray.origin()[0];
    rayOrg[1] = ray.origin()[1];
    rayOrg[2] = ray.origin()[2];

    double3 rayDir;
    rayDir[0] = ray.direction()[0];
    rayDir[1] = ray.direction()[1];
    rayDir[2] = ray.direction()[2];

    for (unsigned int i = 0; i < numTriangles; i++) {
      int faceIdx = indices[i + offset];

      // self-intersection check
      if (faceIdx == ray.prev_prim_id) {
        continue;
      }

      int f0 = mesh->faces[3 * faceIdx + 0];
      int f1 = mesh->faces[3 * faceIdx + 1];
      int f2 = mesh->faces[3 * faceIdx + 2];

      vector3d v0, v1, v2;
      v0[0] = mesh->dvertices[3 * f0 + 0];
      v0[1] = mesh->dvertices[3 * f0 + 1];
      v0[2] = mesh->dvertices[3 * f0 + 2];

      v1[0] = mesh->dvertices[3 * f1 + 0];
      v1[1] = mesh->dvertices[3 * f1 + 1];
      v1[2] = mesh->dvertices[3 * f1 + 2];

      v2[0] = mesh->dvertices[3 * f2 + 0];
      v2[1] = mesh->dvertices[3 * f2 + 1];
      v2[2] = mesh->dvertices[3 * f2 + 2];

      real u, v;
      if (TriangleIsectD(t, u, v, v0, v1, v2, rayOrg, rayDir, doubleSided)) {
        // Update isect state
        isect.t = t;
        isect.u = u;
        isect.v = v;
        isect.prim_id = faceIdx;
        hit = true;
      }
    }

  } else { // float

#if ENABLE_SIMD_ISECTOR

#if defined(__sparc__) && defined(__HPC_ACE__)

#define N (1) // unroll factor

    // make 0xFFF...FF or 0x000...00
    double doubleSideFlag = doubleSided ? 1.0 : 0.0;
    double2 doubleSidedMask = _mm_cmpgt_pd(
        _mm_set_pd(doubleSideFlag, doubleSideFlag), _mm_setzero_pd());

    int lanes = 2 * N;

    uint32_t numBlocks = numTriangles / lanes;

    // Widen data from float precision to double precision since
    // HPC-ACE is optimized for double and float speed is identical to double.
    double2 rox = _mm_set_pd(ray.origin()[0], ray.origin()[0]);
    double2 roy = _mm_set_pd(ray.origin()[1], ray.origin()[1]);
    double2 roz = _mm_set_pd(ray.origin()[2], ray.origin()[2]);

    double2 rdx = _mm_set_pd(ray.direction()[0], ray.direction()[0]);
    double2 rdy = _mm_set_pd(ray.direction()[1], ray.direction()[1]);
    double2 rdz = _mm_set_pd(ray.direction()[2], ray.direction()[2]);

    double2 tInOut[N];
    double2 tidInOut[N]; // Use as integer(53bit)
    double2 uOut[N];
    double2 vOut[N];

    const double kInf = std::numeric_limits<double>::infinity();

    // Init
    for (int k = 0; k < N; k++) {
      tInOut[k] = _mm_set_pd(kInf, kInf);   // = no hit
      tidInOut[k] = _mm_set_pd(-1.0, -1.0); // = no hit
      uOut[k] = _mm_setzero_pd();
      vOut[k] = _mm_setzero_pd();
    }

    // SIMD process
    for (unsigned int i = 0; i < lanes * numBlocks; i += lanes) {

      for (int k = 0; k < N; k++) {

        int faceIdx0 = indices[i + (2 * k + 0) + offset];
        int faceIdx1 = indices[i + (2 * k + 1) + offset];

        int f0_0 = mesh->faces[3 * faceIdx0 + 0];
        int f1_0 = mesh->faces[3 * faceIdx0 + 1];
        int f2_0 = mesh->faces[3 * faceIdx0 + 2];
        int f0_1 = mesh->faces[3 * faceIdx1 + 0];
        int f1_1 = mesh->faces[3 * faceIdx1 + 1];
        int f2_1 = mesh->faces[3 * faceIdx1 + 2];

        double2 v0x = _mm_set_pd(mesh->vertices[3 * f0_0 + 0],
                                 mesh->vertices[3 * f0_1 + 0]);
        double2 v0y = _mm_set_pd(mesh->vertices[3 * f0_0 + 1],
                                 mesh->vertices[3 * f0_1 + 1]);
        double2 v0z = _mm_set_pd(mesh->vertices[3 * f0_0 + 2],
                                 mesh->vertices[3 * f0_1 + 2]);

        double2 v1x = _mm_set_pd(mesh->vertices[3 * f1_0 + 0],
                                 mesh->vertices[3 * f1_1 + 0]);
        double2 v1y = _mm_set_pd(mesh->vertices[3 * f1_0 + 1],
                                 mesh->vertices[3 * f1_1 + 1]);
        double2 v1z = _mm_set_pd(mesh->vertices[3 * f1_0 + 2],
                                 mesh->vertices[3 * f1_1 + 2]);

        double2 v2x = _mm_set_pd(mesh->vertices[3 * f2_0 + 0],
                                 mesh->vertices[3 * f2_1 + 0]);
        double2 v2y = _mm_set_pd(mesh->vertices[3 * f2_0 + 1],
                                 mesh->vertices[3 * f2_1 + 1]);
        double2 v2z = _mm_set_pd(mesh->vertices[3 * f2_0 + 2],
                                 mesh->vertices[3 * f2_1 + 2]);

        // double has enough precision to store integer(53bit)
        double2 tid = _mm_set_pd(i + 2 * k + 0, i + 2 * k + 1);

        IsectD2(tInOut[k], tidInOut[k], uOut[k], vOut[k], v0x, v0y, v0z, v1x,
                v1y, v1z, v2x, v2y, v2z, rox, roy, roz, rdx, rdy, rdz, tid,
                doubleSidedMask);
      }
    }

    // Merge result.
    for (int i = 0; i < N; i++) {
      double tidRet[2];
      double tRet[2];
      double uRet[2];
      double vRet[2];

      _mm_store_pd(tidRet, tidInOut[i]);
      _mm_store_pd(tRet, tInOut[i]);
      _mm_store_pd(uRet, uOut[i]);
      _mm_store_pd(vRet, vOut[i]);

      for (int k = 0; k < 2; k++) {
        // swap index
        if ((tidRet[1 - k] >= 0.0) && (tRet[1 - k] <= isect.t) &&
            ((unsigned int)(tidRet[1 - k]) !=
             ray.prev_prim_id)) { // +self-intersection check
          // printf("hit. tid = %d\n", (unsigned int)tidRet[k]);
          isect.t = tRet[1 - k];
          isect.u = uRet[1 - k];
          isect.v = vRet[1 - k];
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

      for (unsigned int i = lanes * numBlocks; i < numTriangles; i++) {
        int faceIdx = indices[i + offset];

        // self-intersection check
        if (faceIdx == ray.prev_prim_id) {
          continue;
        }

        int f0 = mesh->faces[3 * faceIdx + 0];
        int f1 = mesh->faces[3 * faceIdx + 1];
        int f2 = mesh->faces[3 * faceIdx + 2];

        vector3 v0, v1, v2;
        v0[0] = mesh->vertices[3 * f0 + 0];
        v0[1] = mesh->vertices[3 * f0 + 1];
        v0[2] = mesh->vertices[3 * f0 + 2];

        v1[0] = mesh->vertices[3 * f1 + 0];
        v1[1] = mesh->vertices[3 * f1 + 1];
        v1[2] = mesh->vertices[3 * f1 + 2];

        v2[0] = mesh->vertices[3 * f2 + 0];
        v2[1] = mesh->vertices[3 * f2 + 1];
        v2[2] = mesh->vertices[3 * f2 + 2];

        real u, v;
        if (TriangleIsect(t, u, v, v0, v1, v2, rayOrg, rayDir, doubleSided)) {
          // Update isect state
          isect.t = t;
          isect.u = u;
          isect.v = v;
          isect.prim_id = faceIdx;
          hit = true;
        }
      }
    }

#undef N

#elif defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

#define N (1) // unroll factor

    // make 0xFFF...FF or 0x000...00
    double doubleSideFlag = doubleSided ? 1.0 : 0.0;
    double2 doubleSidedMask = _mm_cmpgt_pd(
        _mm_set_pd(doubleSideFlag, doubleSideFlag), _mm_setzero_pd());

    int lanes = 2 * N;

    uint32_t numBlocks = numTriangles / lanes;

    // Widen data from float precision to double precision since
    // HPC-ACE is optimized for double and float speed is identical to double.
    double2 rox = _mm_set_pd(ray.origin()[0], ray.origin()[0]);
    double2 roy = _mm_set_pd(ray.origin()[1], ray.origin()[1]);
    double2 roz = _mm_set_pd(ray.origin()[2], ray.origin()[2]);

    double2 rdx = _mm_set_pd(ray.direction()[0], ray.direction()[0]);
    double2 rdy = _mm_set_pd(ray.direction()[1], ray.direction()[1]);
    double2 rdz = _mm_set_pd(ray.direction()[2], ray.direction()[2]);

    double2 tInOut[N];
    double2 tidInOut[N]; // Use as integer(53bit)
    double2 uOut[N];
    double2 vOut[N];

    const double kInf = std::numeric_limits<double>::infinity();

    // Init
    for (int k = 0; k < N; k++) {
      tInOut[k] = _mm_set_pd(kInf, kInf);   // = no hit
      tidInOut[k] = _mm_set_pd(-1.0, -1.0); // = no hit
      uOut[k] = _mm_setzero_pd();
      vOut[k] = _mm_setzero_pd();
    }

    // SIMD process
    for (unsigned int i = 0; i < lanes * numBlocks; i += lanes) {

      for (int k = 0; k < N; k++) {

        int faceIdx0 = indices[i + (2 * k + 0) + offset];
        int faceIdx1 = indices[i + (2 * k + 1) + offset];

        int f0_0 = mesh->faces[3 * faceIdx0 + 0];
        int f1_0 = mesh->faces[3 * faceIdx0 + 1];
        int f2_0 = mesh->faces[3 * faceIdx0 + 2];
        int f0_1 = mesh->faces[3 * faceIdx1 + 0];
        int f1_1 = mesh->faces[3 * faceIdx1 + 1];
        int f2_1 = mesh->faces[3 * faceIdx1 + 2];

        double2 v0x = _mm_set_pd(mesh->vertices[3 * f0_0 + 0],
                                 mesh->vertices[3 * f0_1 + 0]);
        double2 v0y = _mm_set_pd(mesh->vertices[3 * f0_0 + 1],
                                 mesh->vertices[3 * f0_1 + 1]);
        double2 v0z = _mm_set_pd(mesh->vertices[3 * f0_0 + 2],
                                 mesh->vertices[3 * f0_1 + 2]);

        double2 v1x = _mm_set_pd(mesh->vertices[3 * f1_0 + 0],
                                 mesh->vertices[3 * f1_1 + 0]);
        double2 v1y = _mm_set_pd(mesh->vertices[3 * f1_0 + 1],
                                 mesh->vertices[3 * f1_1 + 1]);
        double2 v1z = _mm_set_pd(mesh->vertices[3 * f1_0 + 2],
                                 mesh->vertices[3 * f1_1 + 2]);

        double2 v2x = _mm_set_pd(mesh->vertices[3 * f2_0 + 0],
                                 mesh->vertices[3 * f2_1 + 0]);
        double2 v2y = _mm_set_pd(mesh->vertices[3 * f2_0 + 1],
                                 mesh->vertices[3 * f2_1 + 1]);
        double2 v2z = _mm_set_pd(mesh->vertices[3 * f2_0 + 2],
                                 mesh->vertices[3 * f2_1 + 2]);

        // double has enough precision to store integer(53bit)
        double2 tid = _mm_set_pd(i + 2 * k + 0, i + 2 * k + 1);

        IsectD2(tInOut[k], tidInOut[k], uOut[k], vOut[k], v0x, v0y, v0z, v1x,
                v1y, v1z, v2x, v2y, v2z, rox, roy, roz, rdx, rdy, rdz, tid,
                doubleSidedMask);
      }
    }

    // Merge result.
    for (int i = 0; i < N; i++) {
      double tidRet[2];
      double tRet[2];
      double uRet[2];
      double vRet[2];

      _mm_store_pd(tidRet, tidInOut[i]);
      _mm_store_pd(tRet, tInOut[i]);
      _mm_store_pd(uRet, uOut[i]);
      _mm_store_pd(vRet, vOut[i]);

      for (int k = 0; k < 2; k++) {
        // swap index
        if ((tidRet[1 - k] >= 0.0) && (tRet[1 - k] <= isect.t)) {
          uint32_t prim_id = indices[(unsigned int)tidRet[1 - k] + offset];
          if ((prim_id != ray.prev_prim_id)) { // +self-intersection check
            isect.t = tRet[1 - k];
            isect.u = uRet[1 - k];
            isect.v = vRet[1 - k];
            isect.prim_id = prim_id;
            hit = true;
            t = isect.t;
          }
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

      for (unsigned int i = lanes * numBlocks; i < numTriangles; i++) {
        int faceIdx = indices[i + offset];

        // self-intersection check
        if (faceIdx == ray.prev_prim_id) {
          continue;
        }

        int f0 = mesh->faces[3 * faceIdx + 0];
        int f1 = mesh->faces[3 * faceIdx + 1];
        int f2 = mesh->faces[3 * faceIdx + 2];

        vector3 v0, v1, v2;
        v0[0] = mesh->vertices[3 * f0 + 0];
        v0[1] = mesh->vertices[3 * f0 + 1];
        v0[2] = mesh->vertices[3 * f0 + 2];

        v1[0] = mesh->vertices[3 * f1 + 0];
        v1[1] = mesh->vertices[3 * f1 + 1];
        v1[2] = mesh->vertices[3 * f1 + 2];

        v2[0] = mesh->vertices[3 * f2 + 0];
        v2[1] = mesh->vertices[3 * f2 + 1];
        v2[2] = mesh->vertices[3 * f2 + 2];

        real u, v;
        if (TriangleIsect(t, u, v, v0, v1, v2, rayOrg, rayDir, doubleSided)) {
          // Update isect state
          isect.t = t;
          isect.u = u;
          isect.v = v;
          isect.prim_id = faceIdx;
          hit = true;
        }
      }
    }

#undef N

#endif

#else // !ENABLE_SIMD_ISECTOR

    real3 rayOrg;
    rayOrg[0] = ray.origin()[0];
    rayOrg[1] = ray.origin()[1];
    rayOrg[2] = ray.origin()[2];

    real3 rayDir;
    rayDir[0] = ray.direction()[0];
    rayDir[1] = ray.direction()[1];
    rayDir[2] = ray.direction()[2];

    for (unsigned int i = 0; i < numTriangles; i++) {
      int faceIdx = indices[i + offset];

      // self-intersection check
      if (faceIdx == ray.prev_prim_id) {
        continue;
      }

      int f0 = mesh->faces[3 * faceIdx + 0];
      int f1 = mesh->faces[3 * faceIdx + 1];
      int f2 = mesh->faces[3 * faceIdx + 2];

      vector3 v0, v1, v2;
      v0[0] = mesh->vertices[3 * f0 + 0];
      v0[1] = mesh->vertices[3 * f0 + 1];
      v0[2] = mesh->vertices[3 * f0 + 2];

      v1[0] = mesh->vertices[3 * f1 + 0];
      v1[1] = mesh->vertices[3 * f1 + 1];
      v1[2] = mesh->vertices[3 * f1 + 2];

      v2[0] = mesh->vertices[3 * f2 + 0];
      v2[1] = mesh->vertices[3 * f2 + 1];
      v2[2] = mesh->vertices[3 * f2 + 2];

      real u, v;
      if (TriangleIsect(t, u, v, v0, v1, v2, rayOrg, rayDir, doubleSided)) {
        // Update isect state
        isect.t = t;
        isect.u = u;
        isect.v = v;
        isect.prim_id = faceIdx;
        hit = true;
      }
    }
#endif
  }

  return hit;
}

void BuildIntersection(Intersection &isect, const Mesh *mesh, Ray &ray) {
  // face index
  const unsigned int *faces = mesh->faces;

  isect.f0 = faces[3 * isect.prim_id + 0];
  isect.f1 = faces[3 * isect.prim_id + 1];
  isect.f2 = faces[3 * isect.prim_id + 2];

  if (mesh->isDoublePrecisionPos) {
    const double *vertices = mesh->dvertices;

    vector3d p0, p1, p2;
    p0[0] = vertices[3 * isect.f0 + 0];
    p0[1] = vertices[3 * isect.f0 + 1];
    p0[2] = vertices[3 * isect.f0 + 2];
    p1[0] = vertices[3 * isect.f1 + 0];
    p1[1] = vertices[3 * isect.f1 + 1];
    p1[2] = vertices[3 * isect.f1 + 2];
    p2[0] = vertices[3 * isect.f2 + 0];
    p2[1] = vertices[3 * isect.f2 + 1];
    p2[2] = vertices[3 * isect.f2 + 2];

    // calc shading point.
    isect.position = ray.origin() + isect.t * ray.direction();

    // calc geometric normal.
    vector3d p10 = p1 - p0;
    vector3d p20 = p2 - p0;
    vector3d n = cross(p10, p20);
    n.normalize();

    isect.geometric = n;
    isect.normal = n;

  } else {
    const float *vertices = mesh->vertices;

    vector3 p0, p1, p2;
    p0[0] = vertices[3 * isect.f0 + 0];
    p0[1] = vertices[3 * isect.f0 + 1];
    p0[2] = vertices[3 * isect.f0 + 2];
    p1[0] = vertices[3 * isect.f1 + 0];
    p1[1] = vertices[3 * isect.f1 + 1];
    p1[2] = vertices[3 * isect.f1 + 2];
    p2[0] = vertices[3 * isect.f2 + 0];
    p2[1] = vertices[3 * isect.f2 + 1];
    p2[2] = vertices[3 * isect.f2 + 2];

    // calc shading point.
    isect.position = ray.origin() + isect.t * ray.direction();

    // calc geometric normal.
    vector3 p10 = p1 - p0;
    vector3 p20 = p2 - p0;
    vector3 n = cross(p10, p20);
    n.normalize();

    isect.geometric = n;
    isect.normal = n;
  }
}

} // namespace

bool BVHAccel::Traverse(Intersection &isect, Ray &ray) const {
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

  const double2 vtmin = _mm_set_pd(0.0, 0.0);
  double2 vtmax = _mm_set_pd(hitT, hitT);

#endif

  real minT, maxT;
  while (nodeStackIndex >= 0) {
    int index = nodeStack[nodeStackIndex];
    const BVHNode &node = nodes_[index];

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

        int mask[2] = {(hitmask >> 1) & 0x1, (hitmask) & 0x1};

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

      if (TestLeafNode(isect, node, indices_, mesh_, ray)) {
        hitT = isect.t;

#if ENABLE_SIMD_ISECTOR
        // Shorten t max
        vtmax = _mm_set_pd(hitT, hitT);
#endif
      }
    }
  }

  assert(nodeStackIndex < kMaxStackDepth);

  if (isect.t < std::numeric_limits<real>::max()) {
    BuildIntersection(isect, mesh_, ray);
    return true;
  }

  return false;
}

void BVHAccel::BoundingBox(double bmin[3], double bmax[3]) const {
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
