/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <vector>

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <limits>
#include <functional>
#include <algorithm>

#include "render_accel_line.h"

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

using namespace lsgl::render;

namespace {

#if RENDER_USE_DOUBLE_PRECISION
const real kEPS = 2.0e-15;
#else
const real kEPS = 1.0e-6f;
#endif

typedef struct {
  real t;
  real u;
  real v;
  bool hitCap;
  unsigned int prim_id;
} CylinderIntersection;

inline real GetRadius(const Lines *lines, unsigned int index) {
  if (lines->radius) {
    return lines->radius[index];
  } else {
    return lines->constantRadius;
  }
}

inline real3 GetPosition(const Lines *lines, unsigned int index) {
  return real3(lines->positions[3 * index + 0], lines->positions[3 * index + 1],
               lines->positions[3 * index + 2]);
}

inline real GetMin(real a, real b) { return (a < b) ? a : b; }

inline real GetMax(real a, real b) { return (a > b) ? a : b; }

inline void GetBoundingBoxOfLine(real3 &bmin, real3 &bmax, const Lines *lines,
                                 unsigned int index) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 16.0;

  real3 p0 = GetPosition(lines, lines->segments[2 * index + 0]);
  real3 p1 = GetPosition(lines, lines->segments[2 * index + 1]);

  real r0 = GetRadius(lines, lines->segments[2 * index + 0]);
  real r1 = GetRadius(lines, lines->segments[2 * index + 1]);

  real3 points[4] = {
      p0 - real3(r0, r0, r0), p0 + real3(r0, r0, r0), p1 - real3(r1, r1, r1),
      p1 + real3(r1, r1, r1),
  };

  bmin = points[0];
  bmax = points[0];

  for (int i = 1; i < 4; i++) {
    real3 tp = points[i];
    for (int j = 0; j < 3; j++) {
      bmin[j] = GetMin(bmin[j], tp[j]);
      bmax[j] = GetMax(bmax[j], tp[j]);
    }
  }

  // add epsilon
  for (int j = 0; j < 3; j++) {
    bmin[j] -= bmin[j]*kEPS;
    bmax[j] += bmax[j]*kEPS;
  }
}

class PartPred : public std::unary_function<unsigned int, bool> {
public:
  PartPred(int axis, real pos, const Lines *lines)
      : axis_(axis), pos_(pos), lines_(lines) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;

    real3 p0 = GetPosition(lines_, lines_->segments[2 * i + 0]);
    real3 p1 = GetPosition(lines_, lines_->segments[2 * i + 1]);

    real3 p = (p0 + p1) * 0.5;

    real center = p[axis];

    return (center < pos_);
  }

private:
  int axis_;
  real pos_;
  const Lines *lines_;
};

} // namespace

static void ComputeBoundingBox(real3 &bmin, real3 &bmax, const Lines *lines,
                               unsigned int *indices, unsigned int leftIndex,
                               unsigned int rightIndex) {
  size_t i = leftIndex;
  size_t idx = indices[i];
  GetBoundingBoxOfLine(bmin, bmax, lines, idx);

  for (i = leftIndex; i < rightIndex; i++) {
    idx = indices[i];

    real3 cmin;
    real3 cmax;
    GetBoundingBoxOfLine(cmin, cmax, lines, idx); // Get child Bounding Box.

    for (int k = 0; k < 3; k++) { // xyz
      bmin[k] = GetMin(bmin[k], cmin[k]);
      bmax[k] = GetMax(bmax[k], cmax[k]);
    }
  }
}

//
// --
//

size_t LineAccel::BuildTree(const Lines *lines, unsigned int leftIdx,
                            unsigned int rightIdx, int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  size_t offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
    stats_.maxTreeDepth = depth;
  }

  real3 bmin, bmax;
  ComputeBoundingBox(bmin, bmax, lines, &indices_.at(0), leftIdx, rightIdx);

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {
    // Create leaf node.
    LineNode leaf;

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
    unsigned int *end = begin + (rightIdx - leftIdx); // &indices_[rightIdx];
    unsigned int *mid = 0;

    // Try longest axis first.
    cutAxis = (longestAxis + axisTry) % 3;
    real cutPos = bmin[cutAxis] + 0.5 * bsize[cutAxis]; // spatial median

    //
    // Split at (cutAxis, cutPos)
    // indices_ will be modified.
    //
    mid = std::partition(begin, end, PartPred(cutAxis, cutPos, lines));

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

  LineNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch
  nodes_.push_back(node);

  // Recurively split tree.
  unsigned int leftChildIndex = BuildTree(lines, leftIdx, midIdx, depth + 1);
  unsigned int rightChildIndex = BuildTree(lines, midIdx, rightIdx, depth + 1);

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

bool LineAccel::Build(const Lines *lines, const LineBuildOptions &options) {
  options_ = options;
  stats_ = LineBuildStatistics();

  assert(lines);

  // @todo { double precision }
  assert(lines->isDoublePrecisionPos == false);

  size_t n = lines->numLines;
  trace("[LineAccel] Input # of lines = %lu\n", n);
  trace("[LineAccel] Cap = %d\n", (int)options.cap);

  assert(n > 0);

  //
  // 1. Create line indices(this will be permutated in BuildTree)
  //
  indices_.resize(n);
  for (size_t i = 0; i < n; i++) {
    indices_[i] = i;
  }

  //
  // 2. Build tree
  //
  BuildTree(lines, 0, n, 0);

  // Tree will be null if input line count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    real3 bmin(&nodes_[0].bmin[0]);
    real3 bmax(&nodes_[0].bmax[0]);
    trace("[LineAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
    trace("[LineAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);
  }

  trace("[LineAccel] # of nodes = %lu\n", nodes_.size());

  // Save ptr
  lines_ = lines;

  return true;
}

bool LineAccel::Dump(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "[LineAccel] Cannot write a file: %s\n", filename);
    return false;
  }

  unsigned long long numNodes = nodes_.size();
  assert(nodes_.size() > 0);

  unsigned long long numIndices = indices_.size();

  int r = 0;
  r = fwrite(&numNodes, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&nodes_.at(0), sizeof(LineNode), numNodes, fp);
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

    tminOut = tmin;
    tmaxOut = tmax;

    return true;
  }

  return false; // no hit
}

inline real3 vcross(real3 a, real3 b) {
  real3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[1] * b[2] - a[2] * b[1];
  c[2] = a[1] * b[2] - a[2] * b[1];
  return c;
}

inline real vdot(real3 a, real3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline real3 vneg(real3 a) {
  real3 c(-a[0], -a[1], -a[2]);
  return c;
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

int solve2e(real root[], real A, real B, real C) {
  if (fabs(A) <= std::numeric_limits<real>::epsilon()) {
    real x = -C / B;
    root[0] = x;
    return 1;
  } else {
    real D = B * B - A * C;
    if (D < 0) {
      return 0;
    } else if (D == 0) {
      real x = -B / A;
      root[0] = x;
      return 1;
    } else {
      real x1 = (fabs(B) + sqrt(D)) / A;
      if (B >= 0.0) {
        x1 = -x1;
      }
      real x2 = C / (A * x1);
      if (x1 > x2) {
        real tmp = x1;
        x1 = x2;
        x2 = tmp;
      }

      root[0] = x1;
      root[1] = x2;
      return 2;
    }
  }
}

inline bool CylinderIsect(real *pt, real *pu, real *pv, bool &hitCap,
                          const real3 &p0, real r0, const real3 &p1, real r1,
                          const real3 &rayOrg, const real3 &rayDir, bool cap) {

  real tmax = *pt;
  real rr = std::max<real>(r0, r1);
  real3 ORG = rayOrg;
  real3 n = rayDir;
  real3 d = p1 - p0;
  real3 m = ORG - p0;

  real md = vdot(m, d);
  real nd = vdot(n, d);
  real dd = vdot(d, d);

  hitCap = false;
  real capT = (std::numeric_limits<real>::max)(); // far

  if (cap) {
    real3 dN0 = vnormalize(p0 - p1);
    real3 dN1 = vneg(dN0);
    real3 rd = vnormalize(rayDir);

    if (fabs(vdot(rayDir, dN0)) > kEPS) {

      // test with 2 planes
      real p0D = -vdot(p0, dN0); // plane D
      real p1D = -vdot(p1, dN1); // plane D

      real p0T = -(vdot(rayOrg, dN0) + p0D) / vdot(rd, dN0);
      real p1T = -(vdot(rayOrg, dN1) + p1D) / vdot(rd, dN1);

      real3 q0 = rayOrg + p0T * rd;
      real3 q1 = rayOrg + p1T * rd;

      real qp0Sqr = vdot(q0 - p0, q0 - p0);
      real qp1Sqr = vdot(q1 - p1, q1 - p1);
      // printf("p0T = %f, p1T = %f, q0Sqr = %f, rr^2 = %f\n", p0T, p1T, q0Sqr,
      // rr*rr);

      if (p0T > 0.0 && p0T < tmax && (qp0Sqr < rr * rr)) {
        // hit p0's plane
        hitCap = true;
        capT = p0T;
        *pt = capT;
        *pu = sqrt(qp0Sqr);
        *pv = 0;
      }

      if (p1T > 0.0 && p1T < tmax && p1T < capT && (qp1Sqr < rr * rr)) {
        hitCap = true;
        capT = p1T;
        *pt = capT;
        *pu = sqrt(qp1Sqr);
        *pv = 1.0;
      }
    }
  }

  if (md <= 0.0 && nd <= 0.0)
    return hitCap;
  if (md >= dd && nd >= 0.0)
    return hitCap;

  real nn = vdot(n, n);
  real mn = vdot(m, n);
  real A = dd * nn - nd * nd;
  real k = vdot(m, m) - rr * rr;
  real C = dd * k - md * md;
  real B = dd * mn - nd * md;

  real root[2] = {};
  int nRet = solve2e(root, A, B, C);
  if (nRet) {
    real t = root[0];
    if (0 <= t && t <= tmax && t <= capT) {
      real s = md + t * nd;
      s /= dd;
      if (0 <= s && s <= 1) {
        hitCap = false;
        *pt = t;
        *pu = 0;
        *pv = s;

        return true;
      }
    }
  }
  return hitCap;
}

bool TestLeafNode(CylinderIntersection &isect, // [inout]
                  const LineNode &node,
                  const std::vector<unsigned int> &indices, const Lines *lines,
                  const Ray &ray, bool cap) {

  bool hit = false;

  unsigned int numLines = node.data[0];
  unsigned int offset = node.data[1];

  real t = isect.t; // current hit distance

  for (unsigned int i = 0; i < numLines; i++) {
    unsigned int index = indices[i + offset];

    real3 p0 = GetPosition(lines, lines->segments[2 * index + 0]);
    real3 p1 = GetPosition(lines, lines->segments[2 * index + 1]);

    real r0 = GetRadius(lines, lines->segments[2 * index + 0]);
    real r1 = GetRadius(lines, lines->segments[2 * index + 1]);

    real u;
    real v;

    real3 rayOrg;
    rayOrg[0] = ray.origin()[0];
    rayOrg[1] = ray.origin()[1];
    rayOrg[2] = ray.origin()[2];

    real3 rayDir;
    rayDir[0] = ray.direction()[0];
    rayDir[1] = ray.direction()[1];
    rayDir[2] = ray.direction()[2];

    bool hitCap = false;
    if (CylinderIsect(&t, &u, &v, hitCap, p0, r0, p1, r1, rayOrg, rayDir,
                      cap)) {
      // Update isect state
      isect.t = t;
      isect.u = u;
      isect.v = v;
      isect.hitCap = hitCap;

      // u and v are computed in later.
      isect.prim_id = index;
      hit = true;
    }
  }

  return hit;
}

#if 0
//
// Z-up
//
inline void XYZToUV(real &theta, real &phi, const real3 &v) {
  phi = atan2(v[1], v[0]);
  if (phi < 0.0) {
    phi += 2.0 * M_PI; // [0, 2PI]
  }
  if (phi < 0.0)
    phi = 0.0;
  if (phi > 2.0 * M_PI)
    phi = 2.0 * M_PI;

  real z = v[2];
  if (z < -1.0)
    z = -1.0;
  if (z > 1.0)
    z = 1.0;
  theta = acos(z);
}
#endif

void BuildIntersection(Intersection &isectRet,
                       const CylinderIntersection &isect, const Lines *lines,
                       Ray &ray) {

  isectRet.t = isect.t;
  isectRet.u = isect.u;
  isectRet.v = isect.v;
  isectRet.prim_id = isect.prim_id;

  real v = isect.v;
  unsigned int index = isect.prim_id;
  real3 p0 = GetPosition(lines, lines->segments[2 * index + 0]);
  real3 p1 = GetPosition(lines, lines->segments[2 * index + 1]);

  real3 center = p0 + real3(v, v, v) * (p1 - p0);

  isectRet.position[0] = ray.origin()[0] + isect.t * ray.direction()[0];
  isectRet.position[1] = ray.origin()[1] + isect.t * ray.direction()[1];
  isectRet.position[2] = ray.origin()[2] + isect.t * ray.direction()[2];

  real3 n;
  if (isect.hitCap) {
    real3 p(isectRet.position[0], isectRet.position[1], isectRet.position[2]);
    real3 c = 0.5 * (p1 - p0) + p0;
    n = vnormalize(p1 - p0);

    if (vdot((p - c), n) > 0.0) {
      // hit p1's plane
    } else {
      // hit p0's plane
      n = vneg(n);
    }
  } else {
    n[0] = isectRet.position[0] - center[0];
    n[1] = isectRet.position[1] - center[1];
    n[2] = isectRet.position[2] - center[2];
    n = vnormalize(n);
  }
  isectRet.normal[0] = n[0];
  isectRet.normal[1] = n[1];
  isectRet.normal[2] = n[2];

  // Same ID for line data.
  isectRet.f0 = isect.prim_id;
  isectRet.f1 = isect.prim_id;
  isectRet.f2 = isect.prim_id;
}

} // namespace

bool LineAccel::Traverse(Intersection &isectRet, Ray &ray) const {
  real hitT = (std::numeric_limits<real>::max)(); // far = no hit.

  int nodeStackIndex = 0;
  std::vector<int> nodeStack(512);
  nodeStack[0] = 0;
  CylinderIntersection isect;

  bool ret = false;

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
    const LineNode &node = nodes_[index];

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

      if (TestLeafNode(isect, node, indices_, lines_, ray, options_.cap)) {
        hitT = isect.t;
      }
    }
  }

  assert(nodeStackIndex < kMaxStackDepth);

  if (isect.t < (std::numeric_limits<real>::max)()) {
    BuildIntersection(isectRet, isect, lines_, ray);
    return true;
  }

  return false;
}

void LineAccel::BoundingBox(double bmin[3], double bmax[3]) {
  if (nodes_.empty()) {
    bmin[0] = (std::numeric_limits<double>::max)();
    bmin[1] = (std::numeric_limits<double>::max)();
    bmin[2] = (std::numeric_limits<double>::max)();
    bmax[0] = -(std::numeric_limits<double>::max)();
    bmax[1] = -(std::numeric_limits<double>::max)();
    bmax[2] = -(std::numeric_limits<double>::max)();
  } else {
    bmin[0] = nodes_[0].bmin[0];
    bmin[1] = nodes_[0].bmin[1];
    bmin[2] = nodes_[0].bmin[2];
    bmax[0] = nodes_[0].bmax[0];
    bmax[1] = nodes_[0].bmax[1];
    bmax[2] = nodes_[0].bmax[2];
  }
}
