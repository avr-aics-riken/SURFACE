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

#include "accel_tetra.h"
#include "prefix_tree_util.h"
#include "tinymt64.h"

using namespace lsgl::render;

#define MAX_LEAF_ELEMENTS (8)
#define MAX_TREE_DEPTH_32BIT                                                   \
  (22) // FYI, log2(1G/16) ~= 25.897, log2(1G/32) ~= 21

#define ENABLE_TRACE_PRINT (0)
#define ENABLE_DEBUG_PRINT (0)

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

inline int fsign(const real x) {
  real eps = std::numeric_limits<real>::epsilon() * 16;
  return (x > eps ? 1 : (x < -eps ? -1 : 0));
}

//
// Simple Pluecker coordinate class
//
class Pluecker {
public:
  real3 d; // direction
  real3 c; // cross

  Pluecker(const real3 &v0, const real3 &v1) : d(v1 - v0), c(vcross(v1, v0)) {}
};

// Inner product
inline real operator*(const Pluecker &p0, const Pluecker &p1) {
  return vdot(p0.d, p1.c) + vdot(p1.d, p0.c);
}

// Ray - Tetrahedron intersection based on
// - Fast Ray-Tetrahedron Intersection using Plu ̈cker Coordinates
//   Nikos Platis and Theoharis Theoharis, JGT 2003
//   http://realtimecollisiondetection.net/files/PlatisRayTetra/RayTetraPluecker.cpp
//

// Computes the parametric distance tEnter and tLeave
// of enterPoint and leavePoint from orig
// (To enhance performance of the intersection algorithm, this code
// should in practice be incorporated to the function below.)

void ComputeParametricDist(const real3 &orig, const real3 &dir,
                           const real3 &enterPoint, const real3 &leavePoint,
                           double &tEnter, double &tLeave) {
  if (dir.x) {
    double invDirx = 1.0 / dir.x;
    tEnter = (enterPoint.x - orig.x) * invDirx;
    tLeave = (leavePoint.x - orig.x) * invDirx;
  } else if (dir.y) {
    double invDiry = 1.0 / dir.y;
    tEnter = (enterPoint.y - orig.y) * invDiry;
    tLeave = (leavePoint.y - orig.y) * invDiry;
  } else {
    double invDirz = 1.0 / dir.z;
    tEnter = (enterPoint.z - orig.z) * invDirz;
    tLeave = (leavePoint.z - orig.z) * invDirz;
  }
}

bool RayTetraPluecker(const real3 &orig, const real3 &dir, const real3 vert[],
                      int &enterFace, int &leaveFace, real3 &enterPoint,
                      real3 &leavePoint, double &uEnter1, double &uEnter2,
                      double &uLeave1, double &uLeave2, double &tEnter,
                      double &tLeave) {
  enterFace = -1;
  leaveFace = -1;

  double uAB = 0, uAC = 0, uDB = 0, uDC = 0, uBC = 0, uAD = 0;
  // Keep the compiler happy about uninitialized variables
  int signAB = -2, signAC = -2, signDB = -2, signDC = -2, signBC = -2,
      signAD = -2;

  // In the following: A,B,C,D=vert[i], i=0,1,2,3.
  real3 dest = orig + dir;
  Pluecker plRay(orig, dest);

  int nextSign = 0;

  // Examine face ABC
  uAB = plRay * Pluecker(vert[0], vert[1]);
  signAB = fsign(uAB);

  uAC = plRay * Pluecker(vert[0], vert[2]);
  signAC = fsign(uAC);

  if ((signAC == -signAB) || (signAC == 0) || (signAB == 0)) {
    // Face ABC may intersect with the ray
    uBC = plRay * Pluecker(vert[1], vert[2]);
    signBC = fsign(uBC);

    int signABC = signAB;
    if (signABC == 0) {
      signABC = -signAC;
      if (signABC == 0) {
        signABC = signBC;
      }
    }

    if ((signABC != 0) && ((signBC == signABC) || (signBC == 0))) {
      // Face ABC intersects with the ray
      double invVolABC = 1.0 / (uAB + uBC - uAC);
      if (signABC == 1) {
        enterFace = 3;
        uEnter1 = -uAC * invVolABC;
        uEnter2 = uAB * invVolABC;
        enterPoint = (1 - uEnter1 - uEnter2) * vert[0] + uEnter1 * vert[1] +
                     uEnter2 * vert[2];

        nextSign = -1;
      } else {
        leaveFace = 3;
        uLeave1 = -uAC * invVolABC;
        uLeave2 = uAB * invVolABC;
        leavePoint = (1 - uLeave1 - uLeave2) * vert[0] + uLeave1 * vert[1] +
                     uLeave2 * vert[2];

        nextSign = 1;
      }

      // Determine the other intersecting face between BAD, CDA, DCB
      // Examine face BAD
      uAD = plRay * Pluecker(vert[0], vert[3]);
      signAD = fsign(uAD);

      if ((signAD == nextSign) || (signAD == 0)) {
        // Face BAD may intersect with the ray
        uDB = plRay * Pluecker(vert[3], vert[1]);
        signDB = fsign(uDB);

        if ((signDB == nextSign) ||
            ((signDB == 0) && ((signAD != 0) || (signAB != 0)))) {
          // Face BAD intersects with the ray
          double invVolBAD = 1.0 / (uAD + uDB - uAB);
          if (nextSign == 1) {
            enterFace = 2;
            uEnter1 = uDB * invVolBAD;
            uEnter2 = -uAB * invVolBAD;
            enterPoint = (1 - uEnter1 - uEnter2) * vert[1] + uEnter1 * vert[0] +
                         uEnter2 * vert[3];

            ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                  tLeave);
            return true;
          } else {
            leaveFace = 2;
            uLeave1 = uDB * invVolBAD;
            uLeave2 = -uAB * invVolBAD;
            leavePoint = (1 - uLeave1 - uLeave2) * vert[1] + uLeave1 * vert[0] +
                         uLeave2 * vert[3];

            ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                  tLeave);
            return true;
          }
        }
      }

      // Face BAD does not intersect with the ray.
      // Determine the other intersecting face between CDA, DCB
      uDC = plRay * Pluecker(vert[3], vert[2]);
      signDC = fsign(uDC);

      if ((signDC == -nextSign) ||
          ((signDC == 0) && ((signAD != 0) || (signAC != 0)))) {
        // Face CDA intersects with the ray
        double invVolCDA = 1.0 / (uAC - uDC - uAD);
        if (nextSign == 1) {
          enterFace = 1;
          uEnter1 = uAC * invVolCDA;
          uEnter2 = -uDC * invVolCDA;
          enterPoint = (1 - uEnter1 - uEnter2) * vert[2] + uEnter1 * vert[3] +
                       uEnter2 * vert[0];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        } else {
          leaveFace = 1;
          uLeave1 = uAC * invVolCDA;
          uLeave2 = -uDC * invVolCDA;
          leavePoint = (1 - uLeave1 - uLeave2) * vert[2] + uLeave1 * vert[3] +
                       uLeave2 * vert[0];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        }
      } else {
        // Face DCB intersects with the ray
        if (signDB == -2) {
          uDB = plRay * Pluecker(vert[3], vert[1]);
        }

        double invVolDCB = 1.0 / (uDC - uBC - uDB);
        if (nextSign == 1) {
          enterFace = 0;
          uEnter1 = -uDB * invVolDCB;
          uEnter2 = uDC * invVolDCB;
          enterPoint = (1 - uEnter1 - uEnter2) * vert[3] + uEnter1 * vert[2] +
                       uEnter2 * vert[1];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        } else {
          leaveFace = 0;
          uLeave1 = -uDB * invVolDCB;
          uLeave2 = uDC * invVolDCB;
          leavePoint = (1 - uLeave1 - uLeave2) * vert[3] + uLeave1 * vert[2] +
                       uLeave2 * vert[1];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        }
      }
    }
  }

  // Examine face BAD
  uAD = plRay * Pluecker(vert[0], vert[3]);
  signAD = fsign(uAD);

  if ((signAD == -signAB) || (signAB == 0) || (signAD == 0)) {
    // Face BAD may intersect with the ray
    uDB = plRay * Pluecker(vert[3], vert[1]);
    signDB = fsign(uDB);

    int signBAD = -signAB;
    if (signBAD == 0) {
      signBAD = signAD;
      if (signBAD == 0) {
        signBAD = signDB;
      }
    }

    if ((signBAD != 0) && ((signDB == signBAD) || (signDB == 0))) {
      // Face BAD intersects with the ray
      double invVolBAD = 1.0 / (uAD + uDB - uAB);
      if (signBAD == 1) {
        enterFace = 2;
        uEnter1 = uDB * invVolBAD;
        uEnter2 = -uAB * invVolBAD;
        enterPoint = (1 - uEnter1 - uEnter2) * vert[1] + uEnter1 * vert[0] +
                     uEnter2 * vert[3];

        nextSign = -1;
      } else {
        leaveFace = 2;
        uLeave1 = uDB * invVolBAD;
        uLeave2 = -uAB * invVolBAD;
        leavePoint = (1 - uLeave1 - uLeave2) * vert[1] + uLeave1 * vert[0] +
                     uLeave2 * vert[3];

        nextSign = 1;
      }

      // Determine the other intersecting face between CDA, DCB
      uDC = plRay * Pluecker(vert[3], vert[2]);
      signDC = fsign(uDC);

      if ((signDC == -nextSign) ||
          ((signDC == 0) && ((signAD != 0) || (signAC != 0)))) {
        // Face CDA intersects with the ray
        double invVolCDA = 1.0 / (uAC - uDC - uAD);
        if (nextSign == 1) {
          enterFace = 1;
          uEnter1 = uAC * invVolCDA;
          uEnter2 = -uDC * invVolCDA;
          enterPoint = (1 - uEnter1 - uEnter2) * vert[2] + uEnter1 * vert[3] +
                       uEnter2 * vert[0];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        } else {
          leaveFace = 1;
          uLeave1 = uAC * invVolCDA;
          uLeave2 = -uDC * invVolCDA;
          leavePoint = (1 - uLeave1 - uLeave2) * vert[2] + uLeave1 * vert[3] +
                       uLeave2 * vert[0];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        }
      } else {
        // Face DCB intersects with the ray
        if (signBC == -2) {
          uBC = plRay * Pluecker(vert[1], vert[2]);
        }

        double invVolDCB = 1.0 / (uDC - uBC - uDB);
        if (nextSign == 1) {
          enterFace = 0;
          uEnter1 = -uDB * invVolDCB;
          uEnter2 = uDC * invVolDCB;
          enterPoint = (1 - uEnter1 - uEnter2) * vert[3] + uEnter1 * vert[2] +
                       uEnter2 * vert[1];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        } else {
          leaveFace = 0;
          uLeave1 = -uDB * invVolDCB;
          uLeave2 = uDC * invVolDCB;
          leavePoint = (1 - uLeave1 - uLeave2) * vert[3] + uLeave1 * vert[2] +
                       uLeave2 * vert[1];

          ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                                tLeave);
          return true;
        }
      }
    }
  }

  // Examine face CDA
  if ((-signAD == signAC) || (signAC == 0) || (signAD == 0)) {
    // Face CDA may intersect with the ray
    uDC = plRay * Pluecker(vert[3], vert[2]);
    signDC = fsign(uDC);

    int signCDA = signAC;
    if (signCDA == 0) {
      signCDA = -signAD;
      if (signCDA == 0) {
        signCDA = -signDC;
      }
    }

    if ((signCDA != 0) && ((-signDC == signCDA) || (signDC == 0))) {
      // Face CDA intersects with the ray
      // Face DCB also intersects with the ray
      double invVolCDA = 1.0 / (uAC - uDC - uAD);

      if (signBC == -2) {
        uBC = plRay * Pluecker(vert[1], vert[2]);
      }
      if (signDB == -2) {
        uDB = plRay * Pluecker(vert[3], vert[1]);
      }
      double invVolDCB = 1.0 / (uDC - uBC - uDB);

      if (signCDA == 1) {
        enterFace = 1;
        uEnter1 = uAC * invVolCDA;
        uEnter2 = -uDC * invVolCDA;
        enterPoint = (1 - uEnter1 - uEnter2) * vert[2] + uEnter1 * vert[3] +
                     uEnter2 * vert[0];

        leaveFace = 0;
        uLeave1 = -uDB * invVolDCB;
        uLeave2 = uDC * invVolDCB;
        leavePoint = (1 - uLeave1 - uLeave2) * vert[3] + uLeave1 * vert[2] +
                     uLeave2 * vert[1];

        ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                              tLeave);
        return true;
      } else {
        leaveFace = 1;
        uLeave1 = uAC * invVolCDA;
        uLeave2 = -uDC * invVolCDA;
        leavePoint = (1 - uLeave1 - uLeave2) * vert[2] + uLeave1 * vert[3] +
                     uLeave2 * vert[0];

        enterFace = 0;
        uEnter1 = -uDB * invVolDCB;
        uEnter2 = uDC * invVolDCB;
        enterPoint = (1 - uEnter1 - uEnter2) * vert[3] + uEnter1 * vert[2] +
                     uEnter2 * vert[1];

        ComputeParametricDist(orig, dir, enterPoint, leavePoint, tEnter,
                              tLeave);
        return true;
      }
    }
  }

  // Three faces do not intersect with the ray, the fourth will not.
  return false;
}

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

static inline void GetBoundingBoxOfTetra(real3 &bmin, real3 &bmax,
                                         const Tetrahedron *tetras,
                                         unsigned int index) {
  unsigned int f0 = tetras->faces[4 * index + 0];
  unsigned int f1 = tetras->faces[4 * index + 1];
  unsigned int f2 = tetras->faces[4 * index + 2];
  unsigned int f3 = tetras->faces[4 * index + 3];

  real3 p[4];

  if (tetras->isDoublePrecisionPos) {
    p[0] = real3(tetras->dvertices[3 * f0 + 0], tetras->dvertices[3 * f0 + 1],
                 tetras->dvertices[3 * f0 + 2]);
    p[1] = real3(tetras->dvertices[3 * f1 + 0], tetras->dvertices[3 * f1 + 1],
                 tetras->dvertices[3 * f1 + 2]);
    p[2] = real3(tetras->dvertices[3 * f2 + 0], tetras->dvertices[3 * f2 + 1],
                 tetras->dvertices[3 * f2 + 2]);
    p[3] = real3(tetras->dvertices[3 * f3 + 0], tetras->dvertices[3 * f3 + 1],
                 tetras->dvertices[3 * f3 + 2]);
  } else {
    p[0] = real3(tetras->vertices[3 * f0 + 0], tetras->vertices[3 * f0 + 1],
                 tetras->vertices[3 * f0 + 2]);
    p[1] = real3(tetras->vertices[3 * f1 + 0], tetras->vertices[3 * f1 + 1],
                 tetras->vertices[3 * f1 + 2]);
    p[2] = real3(tetras->vertices[3 * f2 + 0], tetras->vertices[3 * f2 + 1],
                 tetras->vertices[3 * f2 + 2]);
    p[3] = real3(tetras->vertices[3 * f3 + 0], tetras->vertices[3 * f3 + 1],
                 tetras->vertices[3 * f3 + 2]);
  }

  bmin = p[0];
  bmax = p[0];

  for (int i = 1; i < 4; i++) {
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
                                const Tetrahedron *tetras,
                                unsigned int *indices, unsigned int leftIdx,
                                unsigned int rightIdx) {
  static const real EPS = std::numeric_limits<real>::epsilon() * 16;

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

    GetBoundingBoxOfTetra(bmin, bmax, tetras, indices[i]);

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
                                 const real3 &bmax, size_t numPrims,
                                 real costTaabb) // should be in [0.0, 1.0]
{
  const real eps = std::numeric_limits<real>::epsilon() * 16;

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
    right = numPrims;
    bminLeft = bminRight = bmin;
    bmaxLeft = bmaxRight = bmax;

    for (int i = 0; i < bins->binSize - 1; ++i) {
      left += bins->bin[0 * (3 * bins->binSize) + j * bins->binSize + i];
      right -= bins->bin[1 * (3 * bins->binSize) + j * bins->binSize + i];

      assert(left <= numPrims);
      assert(right <= numPrims);

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
  SAHPred(int axis, real pos, const Tetrahedron *tetras)
      : axis_(axis), pos_(pos), tetras_(tetras) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;
    real pos = pos_;

    unsigned int i0 = tetras_->faces[4 * i + 0];
    unsigned int i1 = tetras_->faces[4 * i + 1];
    unsigned int i2 = tetras_->faces[4 * i + 2];
    unsigned int i3 = tetras_->faces[4 * i + 3];

    real3 p0(tetras_->vertices[3 * i0 + 0], tetras_->vertices[3 * i0 + 1],
             tetras_->vertices[3 * i0 + 2]);
    real3 p1(tetras_->vertices[3 * i1 + 0], tetras_->vertices[3 * i1 + 1],
             tetras_->vertices[3 * i1 + 2]);
    real3 p2(tetras_->vertices[3 * i2 + 0], tetras_->vertices[3 * i2 + 1],
             tetras_->vertices[3 * i2 + 2]);
    real3 p3(tetras_->vertices[3 * i3 + 0], tetras_->vertices[3 * i3 + 1],
             tetras_->vertices[3 * i3 + 2]);

    real center = p0[axis] + p1[axis] + p2[axis] + p3[axis];

    return (center < pos * 4.0);
  }

private:
  int axis_;
  real pos_;
  const Tetrahedron *tetras_;
};

class SAHPredD : public std::unary_function<unsigned int, bool> {
public:
  SAHPredD(int axis, double pos, const Tetrahedron *tetras)
      : axis_(axis), pos_(pos), tetras_(tetras) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;
    double pos = pos_;

    unsigned int i0 = tetras_->faces[4 * i + 0];
    unsigned int i1 = tetras_->faces[4 * i + 1];
    unsigned int i2 = tetras_->faces[4 * i + 2];
    unsigned int i3 = tetras_->faces[4 * i + 3];

    double3 p0(tetras_->dvertices[3 * i0 + 0], tetras_->dvertices[3 * i0 + 1],
               tetras_->dvertices[3 * i0 + 2]);
    double3 p1(tetras_->dvertices[3 * i1 + 0], tetras_->dvertices[3 * i1 + 1],
               tetras_->dvertices[3 * i1 + 2]);
    double3 p2(tetras_->dvertices[3 * i2 + 0], tetras_->dvertices[3 * i2 + 1],
               tetras_->dvertices[3 * i2 + 2]);
    double3 p3(tetras_->dvertices[3 * i3 + 0], tetras_->dvertices[3 * i3 + 1],
               tetras_->dvertices[3 * i3 + 2]);

    double center = p0[axis] + p1[axis] + p2[axis] + p3[axis];

    return (center < pos * 4.0);
  }

private:
  int axis_;
  double pos_;
  const Tetrahedron *tetras_;
};

template <typename T>
static void ComputeBoundingBox(real3 &bmin, real3 &bmax, const T *vertices,
                               unsigned int *faces, unsigned int *indices,
                               unsigned int leftIndex,
                               unsigned int rightIndex) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 16;

  bmin[0] = std::numeric_limits<real>::max();
  bmin[1] = std::numeric_limits<real>::max();
  bmin[2] = std::numeric_limits<real>::max();
  bmax[0] = -std::numeric_limits<real>::max();
  bmax[1] = -std::numeric_limits<real>::max();
  bmax[2] = -std::numeric_limits<real>::max();

  if (rightIndex <= leftIndex) {
    return;
  }

  size_t i = leftIndex;
  size_t idx = indices[i];
  bmin[0] = vertices[3 * faces[4 * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * faces[4 * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * faces[4 * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * faces[4 * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * faces[4 * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * faces[4 * idx + 0] + 2] + kEPS;

  // Assume tetras are composed of all tetrahedrons
  for (i = leftIndex; i < rightIndex; i++) { // for each faces
    size_t idx = indices[i];
    for (int j = 0; j < 4; j++) { // for each tetra vertex
      size_t fid = faces[4 * idx + j];
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
  const real kEPS = std::numeric_limits<real>::epsilon() * 16;

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
  bmin[0] = vertices[3 * faces[4 * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * faces[4 * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * faces[4 * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * faces[4 * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * faces[4 * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * faces[4 * idx + 0] + 2] + kEPS;

  // Assume tetras are composed of all triangles
  for (i = leftIndex; i < rightIndex; i++) { // for each faces
    size_t idx = keys[i].index;
    for (int j = 0; j < 4; j++) { // for each tetra vertex
      size_t fid = faces[4 * idx + j];
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

  const real kEPS = std::numeric_limits<real>::epsilon() * 16;

  bmin[0] = std::numeric_limits<real>::max();
  bmin[1] = std::numeric_limits<real>::max();
  bmin[2] = std::numeric_limits<real>::max();
  bmax[0] = -std::numeric_limits<real>::max();
  bmax[1] = -std::numeric_limits<real>::max();
  bmax[2] = -std::numeric_limits<real>::max();

  {
    size_t i = leftIndex;
    size_t idx = indices[i];
    bmin[0] = vertices[3 * faces[4 * idx + 0] + 0] - kEPS;
    bmin[1] = vertices[3 * faces[4 * idx + 0] + 1] - kEPS;
    bmin[2] = vertices[3 * faces[4 * idx + 0] + 2] - kEPS;
    bmax[0] = vertices[3 * faces[4 * idx + 0] + 0] + kEPS;
    bmax[1] = vertices[3 * faces[4 * idx + 0] + 1] + kEPS;
    bmax[2] = vertices[3 * faces[4 * idx + 0] + 2] + kEPS;
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

      for (int k = 0; k < 4; k++) { // for each tetra vertex
        real minval_x = vertices[3 * faces[4 * idx + k] + 0] - kEPS;
        real minval_y = vertices[3 * faces[4 * idx + k] + 1] - kEPS;
        real minval_z = vertices[3 * faces[4 * idx + k] + 2] - kEPS;

        real maxval_x = vertices[3 * faces[4 * idx + k] + 0] + kEPS;
        real maxval_y = vertices[3 * faces[4 * idx + k] + 1] + kEPS;
        real maxval_z = vertices[3 * faces[4 * idx + k] + 2] + kEPS;

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
      for (int k = 0; k < 3; k++) { // xyz

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
void MakeLeaf32(TetraNode &leaf, const T *vertices, const uint32_t *faces,
                real3 &bmin, real3 &bmax, const std::vector<IndexKey30> &keys,
                uint32_t leftIndex, uint32_t rightIndex) {

  // 1. Compute leaf AABB
  ComputeBoundingBox30(bmin, bmax, vertices, faces, keys, leftIndex,
                       rightIndex);

  // 2. Create leaf node. `n' may be null(create empty leaf node for that case.)
  int n = rightIndex - leftIndex;

  leaf.bmin[0] = bmin[0];
  leaf.bmin[1] = bmin[1];
  leaf.bmin[2] = bmin[2];

  leaf.bmax[0] = bmax[0];
  leaf.bmax[1] = bmax[1];
  leaf.bmax[2] = bmax[2];

  leaf.flag = 1; // leaf
  leaf.data[0] = n;
  leaf.data[1] = (uint32_t)leftIndex;

  return;
}

//
// Build BVH tree from bottom to up manner
//
template <typename T>
size_t BuildTreeRecursive32(std::vector<TetraNode> &nodes, real3 &bmin,
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
    //printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;
    // if (leftIndex == rightIndex) { // this would be OK. 1 tri in 1 leaf case.
    //  endIndex++;
    //}

    TetraNode leaf;
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

  TetraNode node;
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

  node.bmin[0] = bmin[0];
  node.bmin[1] = bmin[1];
  node.bmin[2] = bmin[2];

  node.bmax[0] = bmax[0];
  node.bmax[1] = bmax[1];
  node.bmax[2] = bmax[2];

  nodes[offset] = node;

  return offset;
}

#ifdef _OPENMP
template <typename T>
size_t BuildTreeRecursive32OMP(std::vector<TetraNode> &nodes, real3 &bmin,
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
    //printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;

    TetraNode leaf;
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

  TetraNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  size_t offset;

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

  if (depth > 6) {

    // Enough number of tasks was launched. Switch to sequential code.
    std::vector<TetraNode> sub_nodes;

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

  node.bmin[0] = bmin[0];
  node.bmin[1] = bmin[1];
  node.bmin[2] = bmin[2];

  node.bmax[0] = bmax[0];
  node.bmax[1] = bmax[1];
  node.bmax[2] = bmax[2];

#pragma omp critical
  { nodes[offset] = node; }

  return offset;
}
#endif

} // namespace

//
// --
//

size_t TetraAccel::BuildTree(const Tetrahedron *tetras, unsigned int leftIdx,
                             unsigned int rightIdx, int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  assert(tetras);
  assert(tetras->faces);

  size_t offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
    stats_.maxTreeDepth = depth;
  }

  real3 bmin, bmax;
  if (tetras->isDoublePrecisionPos) {
    ComputeBoundingBox<double>(bmin, bmax, tetras->dvertices, tetras->faces,
                               &indices_.at(0), leftIdx, rightIdx);
  } else {
    ComputeBoundingBox<float>(bmin, bmax, tetras->vertices, tetras->faces,
                              &indices_.at(0), leftIdx, rightIdx);
  }

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {
    // Create leaf node.
    TetraNode leaf;

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
  ContributeBinBuffer(&bins, bmin, bmax, tetras, &indices_.at(0), leftIdx,
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
    if (tetras->isDoublePrecisionPos) {
      mid = std::partition(begin, end,
                           SAHPredD(cutAxis, cutPos[cutAxis], tetras));
    } else {
      mid =
          std::partition(begin, end, SAHPred(cutAxis, cutPos[cutAxis], tetras));
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

  TetraNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch
  nodes_.push_back(node);

  // Recurively split tree.
  unsigned int leftChildIndex = BuildTree(tetras, leftIdx, midIdx, depth + 1);
  unsigned int rightChildIndex = BuildTree(tetras, midIdx, rightIdx, depth + 1);

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

bool TetraAccel::Build(const Tetrahedron *tetras,
                       const TetraBuildOptions &options) {
  options_ = options;
  stats_ = TetraBuildStatistics();

  assert(options_.binSize > 1);

  assert(tetras);

  size_t n = tetras->numTetrahedrons;
  trace("[TetraAccel] Input # of tetrahedrons    = %lu\n",
        tetras->numTetrahedrons);

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
  BuildTree(tetras, 0, n, 0);

  // Tree will be null if input triangle count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    real3 bmin(&nodes_[0].bmin[0]);
    real3 bmax(&nodes_[0].bmax[0]);
    trace("[TetraAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
    trace("[TetraAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);
  }

  trace("[TetraAccel] # of nodes = %lu\n", nodes_.size());

  // Store pointer for later use.
  tetras_ = tetras;

  return true;
}

bool TetraAccel::Build32(const Tetrahedron *tetras,
                         const TetraBuildOptions &options) {

  options_ = options;
  stats_ = TetraBuildStatistics();

  assert(options_.binSize > 1);
  assert(tetras->isDoublePrecisionPos == false); // @todo

  assert(tetras);

  size_t n = tetras->numTetrahedrons;

  if (n < 1024) {
    // Use non-optimized BVH builder.
    return Build(tetras, options);
  }

  trace("[LSGL] [TetraAccel2] Input # of tetrahedrons    = %lu\n",
        tetras->numTetrahedrons);

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
  trace("[LSGL] [1:indexing] %d msec\n", (int)t.msec());

  {
    timerutil t;
    t.start();

    real3 bmin, bmax;

#ifdef _OPENMP
    ComputeBoundingBoxOMP(bmin, bmax, tetras->vertices, tetras->faces,
                          &indices_.at(0), 0, n);
#else
    ComputeBoundingBox(bmin, bmax, tetras->vertices, tetras->faces,
                       &indices_.at(0), 0, n);
#endif

    t.end();
    trace("[LSGL] [2:scene bbox calculation] %d msec\n", (int)t.msec());

    trace("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    trace("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

    std::vector<uint32_t> codes(n);
    assert(sizeof(real) == sizeof(float));

    {
      timerutil t;
      t.start();

      if (tetras->isDoublePrecisionPos) {
        assert(0); // @todo
      } else {
        CalculateMortonCodesTetraFloat30(&codes.at(0), tetras->vertices,
                                         tetras->faces, bmin, bmax, 0, n);
      }
      t.end();
      trace("[LSGL] [3:morton calculation] %d msec\n", (int)t.msec());

      //for (size_t i = 0; i < n; i++) {
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
      trace("[LSGL] [4:radix sort] %d msec\n", (int)t.msec());

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
      trace("[LSGL] [5:Construct binary radix tree: %d msec\n", (int)t.msec());
    }

    {
      timerutil t;
      t.start();

      nodes_.clear();

      // Explicitly create root node and reserve storage here.
      TetraNode rootNode;
      nodes_.push_back(rootNode);

      bool isLeftLeaf =
          (nodeInfos[0].leftType == NODE_TYPE_LEAF) ? true : false;
      bool isRightLeaf =
          (nodeInfos[0].rightType == NODE_TYPE_LEAF) ? true : false;
      uint32_t midIndex = nodeInfos[0].childIndex;

      real3 leftBMin, leftBMax;
      real3 rightBMin, rightBMax;

      //printf("root: midIndex = %d, range (%d, %d), flag = %d/%d\n", midIndex,
      //0, n-1, isLeftLeaf, isRightLeaf);

#ifdef _OPENMP

      size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

#pragma omp parallel shared(leftChildIndex, rightChildIndex, keys, nodeInfos)
      {
#pragma omp single
        {
          leftChildIndex = BuildTreeRecursive32OMP(
              nodes_, leftBMin, leftBMax, keys, nodeInfos, tetras->vertices,
              tetras->faces, midIndex, 0, midIndex, isLeftLeaf, 0);
        }

#pragma omp single
        {
          rightChildIndex = BuildTreeRecursive32OMP(
              nodes_, rightBMin, rightBMax, keys, nodeInfos, tetras->vertices,
              tetras->faces, midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);
        }
      }
#pragma omp barrier

      assert(leftChildIndex != (size_t)(-1));
      assert(rightChildIndex != (size_t)(-1));

#else

      size_t leftChildIndex = BuildTreeRecursive32(
          nodes_, leftBMin, leftBMax, keys, nodeInfos, tetras->vertices,
          tetras->faces, midIndex, 0, midIndex, isLeftLeaf, 0);
      size_t rightChildIndex = BuildTreeRecursive32(
          nodes_, rightBMin, rightBMax, keys, nodeInfos, tetras->vertices,
          tetras->faces, midIndex + 1, midIndex + 1, n - 1, isRightLeaf, 0);

#endif

      //printf("  leftbmin = %f, %f, %f\n", leftBMin[0], leftBMin[1],
      //leftBMin[2]);
      //printf("  leftbmax = %f, %f, %f\n", leftBMax[0], leftBMax[1],
      //leftBMax[2]);
      //printf("  rightbmin = %f, %f, %f\n", rightBMin[0], rightBMin[1],
      //rightBMin[2]);
      //printf("  rightbmax = %f, %f, %f\n", rightBMax[0], rightBMax[1],
      //rightBMax[2]);
      //printf("  leftIndex = %d, rightIndex = %d\n", leftChildIndex, rightChildIndex);
      //printf("  isLeaf = %d, %d\n", isLeftLeaf, isRightLeaf);

      real3 rootBMin, rootBMax;

      MergeBoundingBox(rootBMin, rootBMax, leftBMin, leftBMax, rightBMin,
                       rightBMax);

      rootNode.axis = GetSplitAxis(keys[0].code);
      rootNode.flag = 0; // = branch

      rootNode.data[0] = leftChildIndex;
      rootNode.data[1] = rightChildIndex;

      rootNode.bmin[0] = rootBMin[0];
      rootNode.bmin[1] = rootBMin[1];
      rootNode.bmin[2] = rootBMin[2];

      rootNode.bmax[0] = rootBMax[0];
      rootNode.bmax[1] = rootBMax[1];
      rootNode.bmax[2] = rootBMax[2];

      // @atomic update
      nodes_[0] = rootNode;

      t.end();
      // printf("root: midIndex = %d, range (%d, %d), flag = %d/%d\n", midIndex,
      // 0, n-1, isLeftLeaf, isRightLeaf);
      // printf("node.total = %d", nodes_.size());
      // printf("[%d] -> (%d, %d)\n", 0, leftChildIndex, rightChildIndex);
      trace("[LSGL] [6:Construct final AABB tree: %d msec\n", (int)t.msec());

      trace("[LSGL]   bmin = %f, %f, %f\n", rootBMin[0], rootBMin[1], rootBMin[2]);
      trace("[LSGL]   bmax = %f, %f, %f\n", rootBMax[0], rootBMax[1], rootBMax[2]);
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
    real3 bmin(&nodes_[0].bmin[0]);
    real3 bmax(&nodes_[0].bmax[0]);
    trace("[TetraAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
    trace("[TetraAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);
  }

  trace("[TetraAccel] # of nodes = %lu\n", nodes_.size());

  // Store pointer for later use.
  tetras_ = tetras;

  return true;
}

bool TetraAccel::Dump(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "[TetraAccel] Cannot write a file: %s\n", filename);
    return false;
  }

  unsigned long long numNodes = nodes_.size();
  assert(nodes_.size() > 0);

  unsigned long long numIndices = indices_.size();

  int r = 0;
  r = fwrite(&numNodes, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&nodes_.at(0), sizeof(TetraNode), numNodes, fp);
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

  //
  // Hit include (tmin == tmax) edge case(hit 2D plane).
  //
  if ((tmax > 0.0) && (tmin <= tmax) && (tmin <= maxT)) {

    // printf("tmin, tmax = %f, %f\n", tmin, tmax);
    tminOut = tmin;
    tmaxOut = tmax;

    return true;
  }

  return false; // no hit
}

bool TestLeafNode(Intersection &isect, // [inout]
                  const TetraNode &node,
                  const std::vector<unsigned int> &indices,
                  const Tetrahedron *tetras, const Ray &ray) {

  bool hit = false;

  unsigned int numTetras = node.data[0];
  unsigned int offset = node.data[1];

  real t = isect.t; // current hit distance

  bool doubleSided = (bool)ray.double_sided;

  if (tetras->isDoublePrecisionPos) {

    assert(0); // @todo

  } else { // float

    real3 rayOrg;
    rayOrg[0] = ray.origin()[0];
    rayOrg[1] = ray.origin()[1];
    rayOrg[2] = ray.origin()[2];

    real3 rayDir;
    rayDir[0] = ray.direction()[0];
    rayDir[1] = ray.direction()[1];
    rayDir[2] = ray.direction()[2];

    for (unsigned int i = 0; i < numTetras; i++) {
      int faceIdx = indices[i + offset];

      // self-intersection check
      if (faceIdx == ray.prev_prim_id) {
        continue;
      }
      int f0 = tetras->faces[4 * faceIdx + 0];
      int f1 = tetras->faces[4 * faceIdx + 1];
      int f2 = tetras->faces[4 * faceIdx + 2];
      int f3 = tetras->faces[4 * faceIdx + 3];

      real3 vtx[4];
      vtx[0][0] = tetras->vertices[3 * f0 + 0];
      vtx[0][1] = tetras->vertices[3 * f0 + 1];
      vtx[0][2] = tetras->vertices[3 * f0 + 2];

      vtx[1][0] = tetras->vertices[3 * f1 + 0];
      vtx[1][1] = tetras->vertices[3 * f1 + 1];
      vtx[1][2] = tetras->vertices[3 * f1 + 2];

      vtx[2][0] = tetras->vertices[3 * f2 + 0];
      vtx[2][1] = tetras->vertices[3 * f2 + 1];
      vtx[2][2] = tetras->vertices[3 * f2 + 2];

      vtx[3][0] = tetras->vertices[3 * f3 + 0];
      vtx[3][1] = tetras->vertices[3 * f3 + 1];
      vtx[3][2] = tetras->vertices[3 * f3 + 2];

      int enterF, leaveF;
      real3 enterP, leaveP;
      double enterU, leaveU;
      double enterV, leaveV;
      double enterT, leaveT;
      if (RayTetraPluecker(rayOrg, rayDir, vtx, enterF, leaveF, enterP, leaveP,
                           enterU, leaveU, enterV, leaveV, enterT, leaveT)) {

        if ((enterT >= 0.0) && (enterT < t)) {
          // Update isect state.
          // @todo { Record leaving point. }
          isect.t = enterT;
          isect.u = enterU;
          isect.v = enterV;
          isect.prim_id = faceIdx;
          isect.subface_id = enterF; // 0,1,2 or 3
          t = enterT;
          hit = true;
        }
      }
    }
  }

  return hit;
}

void BuildIntersection(Intersection &isect, const Tetrahedron *tetras,
                       Ray &ray) {
  // face index
  const unsigned int *faces = tetras->faces;

  isect.f0 = faces[4 * isect.prim_id + 0];
  isect.f1 = faces[4 * isect.prim_id + 1];
  isect.f2 = faces[4 * isect.prim_id + 2];
  isect.f3 = faces[4 * isect.prim_id + 3];

  if (tetras->isDoublePrecisionPos) {
    assert(0); // @todo
#if 0
    const double *vertices = tetras->dvertices;

    real3d p0, p1, p2;
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
    real3d p10 = p1 - p0;
    real3d p20 = p2 - p0;
    real3d n = cross(p10, p20);
    n.normalize();

    isect.geometric = n;
    isect.normal = n;
#endif

  } else {
    const float *vertices = tetras->vertices;

    real3 p[4];
    p[0][0] = vertices[3 * isect.f0 + 0];
    p[0][1] = vertices[3 * isect.f0 + 1];
    p[0][2] = vertices[3 * isect.f0 + 2];
    p[1][0] = vertices[3 * isect.f1 + 0];
    p[1][1] = vertices[3 * isect.f1 + 1];
    p[1][2] = vertices[3 * isect.f1 + 2];
    p[2][0] = vertices[3 * isect.f2 + 0];
    p[2][1] = vertices[3 * isect.f2 + 1];
    p[2][2] = vertices[3 * isect.f2 + 2];
    p[3][0] = vertices[3 * isect.f3 + 0];
    p[3][1] = vertices[3 * isect.f3 + 1];
    p[3][2] = vertices[3 * isect.f3 + 2];

    // calc shading point.
    isect.position = ray.origin() + isect.t * ray.direction();

    // calc geometric normal.
    real3 p0, p1, p2;

    // @fixme { Validation required! }
    int face_table[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 1, 2}};

    p0 = p[face_table[isect.subface_id][0]];
    p1 = p[face_table[isect.subface_id][1]];
    p2 = p[face_table[isect.subface_id][2]];

    real3 p10 = p1 - p0;
    real3 p20 = p2 - p0;
    real3 n = cross(p10, p20);
    // printf("n = %f, %f, %f\n", n[0], n[1], n[2]);
    n.normalize();

    isect.geometric = n;
    isect.normal = n;
  }
}

} // namespace

bool TetraAccel::Traverse(Intersection &isect, Ray &ray) const {
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

  real minT, maxT;
  while (nodeStackIndex >= 0) {
    int index = nodeStack[nodeStackIndex];
    const TetraNode &node = nodes_[index];

    nodeStackIndex--;

    if (node.flag == 0) { // branch node

      bool hit = IntersectRayAABB(minT, maxT, hitT, node.bmin, node.bmax,
                                  rayOrg, rayInvDir, dirSign);

      if (hit) {

        int orderNear = dirSign[node.axis];
        int orderFar = 1 - orderNear;

        // Traverse near first.
        nodeStack[++nodeStackIndex] = node.data[orderFar];
        nodeStack[++nodeStackIndex] = node.data[orderNear];
      }

    } else { // leaf node

      if (TestLeafNode(isect, node, indices_, tetras_, ray)) {
        hitT = isect.t;
      }
    }
  }

  assert(nodeStackIndex < kMaxStackDepth);

  if (isect.t < std::numeric_limits<real>::max()) {
    BuildIntersection(isect, tetras_, ray);
    return true;
  }

  return false;
}

void TetraAccel::BoundingBox(double bmin[3], double bmax[3]) const {
  if (nodes_.empty()) {
    bmin[0] = std::numeric_limits<double>::max();
    bmin[1] = std::numeric_limits<double>::max();
    bmin[2] = std::numeric_limits<double>::max();
    bmax[0] = -std::numeric_limits<double>::max();
    bmax[1] = -std::numeric_limits<double>::max();
    bmax[2] = -std::numeric_limits<double>::max();
  } else {
    bmin[0] = nodes_[0].bmin[0];
    bmin[1] = nodes_[0].bmin[1];
    bmin[2] = nodes_[0].bmin[2];
    bmax[0] = nodes_[0].bmax[0];
    bmax[1] = nodes_[0].bmax[1];
    bmax[2] = nodes_[0].bmax[2];
  }
}
