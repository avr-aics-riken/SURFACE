/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LINE_ACCEL_H__
#define __LINE_ACCEL_H__

#include <vector>

#include "render_common.h"
#include "prim_line.h"
#include "ray.h"
#include "intersection.h"

namespace lsgl {
namespace render {

/// Represents node data for line BVH accelerator.
class LineNode {
public:
  LineNode() {};
  ~LineNode() {};

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

/// Line BVH build option.
struct LineBuildOptions {
  bool debugPrint;
  int minLeafPrimitives;
  int maxTreeDepth;

  // Set default value: Taabb = 0.2
  LineBuildOptions() : minLeafPrimitives(16), maxTreeDepth(512) {}
};

///< Line BVH build statistics.
struct LineBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  LineBuildStatistics() : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// Line BVH accelerator class
class LineAccel {
public:
  LineAccel() {};
  ~LineAccel() {};

  /// Build BVH for input lines.
  bool Build(const Lines *particles, const LineBuildOptions &options);

  /// Get statistics of built BVH tree. Valid after Build()
  LineBuildStatistics GetStatistics() const { return stats_; }

  /// Dump built BVH to the file.
  bool Dump(const char *filename);

  /// Traverse into BVH along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]);

private:
  /// Builds BVH tree recursively.
  size_t BuildTree(const Lines *particles, unsigned int leftIdx,
                   unsigned int rightIdx, int depth);

  LineBuildOptions options_;
  std::vector<LineNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G lines.
  LineBuildStatistics stats_;
  const Lines *lines_;
};

} // namespace render
} // namespace lsgl

#endif // __LINE_ACCEL_H__
