/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_ACCEL_TETRA_H__
#define __LSGL_RENDER_ACCEL_TETRA_H__

#include <vector>

#include "render_common.h"
#include "render_prim_tetra.h"
#include "render_intersection.h"
#include "render_ray.h"

namespace lsgl {
namespace render {

/// Tetra BVH node
class TetraNode {
public:
  TetraNode(){};
  ~TetraNode(){};

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

/// Tetra BVH build option.
struct TetraBuildOptions {
  bool debugPrint;
  real costTaabb;
  int minLeafPrimitives;
  int maxTreeDepth;
  int binSize;

  // Set default value: Taabb = 0.2
  TetraBuildOptions()
      : costTaabb(0.2), minLeafPrimitives(16), maxTreeDepth(256), binSize(64) {}
};

/// Tetra BVH build statistics.
struct TetraBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  TetraBuildStatistics()
      : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// Tetra Acceleration class
class TetraAccel {
public:
  TetraAccel(){};
  ~TetraAccel(){};

  /// Build BVH for input tetra.
  bool Build(const Tetrahedron *tetras, const TetraBuildOptions &options);

  /// Build BVH for input tetra. Optimied for ~2G input
  bool Build32(const Tetrahedron *tetras, const TetraBuildOptions &options);

  /// Get statistics of built BVH tree. Valid after Build()
  TetraBuildStatistics GetStatistics() const { return stats_; }

  /// Dump built BVH to the file.
  bool Dump(const char *filename);

  /// Traverse into BVH along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  // @todo
  /// Traverse into BVH along ray and find all hit points
  // bool TraverseMultiHit(std::vector<Intersection> &isects, Ray &ray) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]) const;

private:
  /// Builds BVH tree recursively.
  size_t BuildTree(const Tetrahedron *tetras, unsigned int leftIdx,
                   unsigned int rightIdx, int depth);

  TetraBuildOptions options_;
  std::vector<TetraNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G triangles.
  TetraBuildStatistics stats_;
  const Tetrahedron *tetras_;
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_ACCEL_TETRA_H__
