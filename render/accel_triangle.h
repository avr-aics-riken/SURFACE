/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_ACCEL_TRIANGLE_H__
#define __LSGL_RENDER_ACCEL_TRIANGLE_H__

#include <vector>

#include "render_common.h"
#include "prim_mesh.h"
#include "render_intersection.h"

namespace lsgl {
namespace render {

/// Triangle BVH node
class TriangleNode {
public:
  TriangleNode(){};
  ~TriangleNode(){};

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

/// Triangle BVH build option.
struct TriangleBuildOptions {
  bool debugPrint;
  real costTaabb;
  int minLeafPrimitives;
  int maxTreeDepth;
  int binSize;

  // Set default value: Taabb = 0.2
  TriangleBuildOptions()
      : costTaabb(0.2), minLeafPrimitives(16), maxTreeDepth(256), binSize(64) {}
};

/// Triangle BVH build statistics.
struct TriangleBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  TriangleBuildStatistics()
      : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// Triangle BVH acceleration class
class TriangleAccel {
public:
  TriangleAccel(){};
  ~TriangleAccel(){};

  /// Build Triangle BVH for input mesh.
  bool Build(const Mesh *mesh, const TriangleBuildOptions &options);

  /// Build Triangle BVH for input mesh. Optimied for ~2G input
  bool Build32(const Mesh *mesh, const TriangleBuildOptions &options);

  /// Get statistics of built Triangle BVH tree. Valid after Build()
  TriangleBuildStatistics GetStatistics() const { return stats_; }

  /// Dump built Triangle BVH to the file.
  bool Dump(const char *filename);

  /// Traverse into Triangle BVH along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  /// Get the bounding box of Triangle BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]) const;

private:
  /// Builds Triangle BVH tree recursively.
  size_t BuildTree(const Mesh *mesh, real3 &bmin, real3 &bmax,
                   unsigned int leftIdx, unsigned int rightIdx, int depth);

  TriangleBuildOptions options_;
  std::vector<TriangleNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G triangles.
  TriangleBuildStatistics stats_;
  const Mesh *mesh_;

  double bmin_[3], bmax_[3]; // Bounding box of whole tree
  double epsScale_; // Epsilon scale(depends on the size of scene bounding box)
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_ACCEL_TRIANGLE_H__
