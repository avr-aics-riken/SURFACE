/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2016 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_ACCEL_SOLID_H__
#define __LSGL_RENDER_ACCEL_SOLID_H__

#include <vector>

#include "render_common.h"
#include "render_prim_solid.h"
#include "render_intersection.h"
#include "render_ray.h"

namespace lsgl {
namespace render {

/// Solid BVH node
class SolidNode {
public:
  SolidNode(){};
  ~SolidNode(){};

  real bmin[3];
  real bmax[3];

  int flag; // 1 = leaf node, 0 = branch node
  int axis;

  // leaf
  //   data[0] = npoints
  //   data[1] = index
  //
  // branch
  //   data[0] = child[0]
  //   data[1] = child[1]
  unsigned int data[2];
};

/// Solid BVH build option.
struct SolidBuildOptions {
  bool debugPrint;
  real costTaabb;
  int minLeafPrimitives;
  int maxTreeDepth;
  int binSize;

  // Set default value: Taabb = 0.2
  SolidBuildOptions()
      : costTaabb(0.2), minLeafPrimitives(8), maxTreeDepth(256), binSize(64) {}
};

/// Solid BVH build statistics.
struct SolidBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  SolidBuildStatistics()
      : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// Solid BVH traversal statistics.
struct SolidTraversalStatistics {
  unsigned long long numRays;
  unsigned long long numLeafTests;
  unsigned long long numPrimIsectTests;
  unsigned long long numNodeTraversals;

  SolidTraversalStatistics()
      : numRays(0), numLeafTests(0), numPrimIsectTests(0), numNodeTraversals(0) {}
};


/// Solid Acceleration class
class SolidAccel {
public:
  SolidAccel(){};
  ~SolidAccel(){};

  /// Build BVH for input tetra.
  bool Build(const Solid *solids, const SolidBuildOptions &options);

  /// Build BVH for input tetra. Optimied for ~2G input
  bool Build32(const Solid *solids, const SolidBuildOptions &options);

  /// Get statistics of built BVH tree. Valid after Build()
  SolidBuildStatistics GetStatistics() const { return buildStats_; }

  /// Dump built BVH to the file.
  bool Dump(const char *filename);

  /// Traverse into BVH along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  // @todo
  /// Traverse into BVH along ray and find all hit points
  // bool TraverseMultiHit(std::vector<Intersection> &isects, Ray &ray) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]) const;

  // These are valid Only if ENABLE_TRAVERSAL_STATISTICS macro was defiend in `render_accel_solid.cc`
  // Since logging traversal information affects the performance in rendering.
  void ResetTraversalStatistics() const;
  void ReportTraversalStatistics() const;

private:
  /// Builds BVH tree recursively.
  size_t BuildTree(const Solid *solids, unsigned int leftIdx,
                   unsigned int rightIdx, int depth);

  SolidBuildOptions options_;
  std::vector<SolidNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G triangles.
  SolidBuildStatistics buildStats_;
  mutable SolidTraversalStatistics traversalStats_;
  const Solid *solids_;
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_ACCEL_SOLID_H__
