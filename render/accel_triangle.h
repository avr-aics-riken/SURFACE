/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __ACCEL_TRIANGLE_H__
#define __ACCEL_TRIANGLE_H__

#include <vector>

#include "render_common.h"
#include "prim_mesh.h"
#include "intersection.h"

namespace lsgl {
namespace render {

/// BVH node
class BVHNode {
public:
  BVHNode() {};
  ~BVHNode() {};

  real bmin[2][3];
  real bmax[2][3];

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

/// BVH build option.
struct BVHBuildOptions {
  bool debugPrint;
  real costTaabb;
  int minLeafPrimitives;
  int maxTreeDepth;
  int binSize;

  // Set default value: Taabb = 0.2
  BVHBuildOptions()
      : costTaabb(0.2), minLeafPrimitives(16), maxTreeDepth(256), binSize(64) {}
};

/// BVH build statistics.
struct BVHBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  BVHBuildStatistics() : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// BVH acceleration class
class BVHAccel {
public:
  BVHAccel() {};
  ~BVHAccel() {};

  /// Build BVH for input mesh.
  bool Build(const Mesh *mesh, const BVHBuildOptions &options);

  /// Build BVH for input mesh. Optimied for ~2G input
  bool Build32(const Mesh *mesh, const BVHBuildOptions &options);

  /// Get statistics of built BVH tree. Valid after Build()
  BVHBuildStatistics GetStatistics() const { return stats_; }

  /// Dump built BVH to the file.
  bool Dump(const char *filename);

  /// Traverse into BVH along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]) const;

private:
  /// Builds BVH tree recursively.
  size_t BuildTree(const Mesh *mesh, real3 &bmin, real3 &bmax,
                   unsigned int leftIdx, unsigned int rightIdx, int depth);

  BVHBuildOptions options_;
  std::vector<BVHNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G triangles.
  BVHBuildStatistics stats_;
  const Mesh *mesh_;

  double bmin_[3], bmax_[3]; // Bounding box of whole tree
};

} // namespace render
} // namespace lsgl

#endif // __BVH_ACCEL_H__
