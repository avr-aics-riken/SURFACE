//
// Sparse volume data structure.
//
#ifndef __LSGL_RENDER_ACCEL_VOLUME_H__
#define __LSGL_RENDER_ACCEL_VOLUME_H__

#include <vector>

#include "render_common.h"
#include "prim_volume.h"
#include "ray.h"
#include "intersection.h"

namespace lsgl {
namespace render {

class BlockNode {
public:
  BlockNode() {};
  ~BlockNode() {};

  // 2-ary BVH
  real bmin[2][3];
  real bmax[2][3];

  // leaf
  //   data[0] = npoints
  //   data[1] = index
  //
  // branch
  //   data[0] = child[0]
  //   data[1] = child[1]
  unsigned int data[2];

  short flag; // 1 = leaf node, 0 = branch node
  short axis;
};

/// SparseVolume build option.
struct SparseVolumeBuildOptions {
  bool debugPrint;
  int minLeafPrimitives;
  int maxTreeDepth;

  // Set default value: Taabb = 0.2
  SparseVolumeBuildOptions() : minLeafPrimitives(16), maxTreeDepth(256) {}
};

/// SparseVolume build statistics.
struct SparseVolumeBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  SparseVolumeBuildStatistics()
      : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// SparseVolume BVH acceleration class.
class SparseVolumeAccel {
public:
  SparseVolumeAccel() {};
  ~SparseVolumeAccel() {};

  /// Build volume BVH for input volume blocks.
  bool Build(const SparseVolume *sparseVolume, const SparseVolumeBuildOptions &options);

  /// Get statistics of built SparseVolume tree. Valid after Build()
  SparseVolumeBuildStatistics GetStatistics() const { return stats_; }

  /// Dump built volume BVH to the file.
  bool Dump(const char *filename);

  /// Traverse into SparseVolume along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]) const;

private:
  /// Builds volume block tree recursively.
  size_t BuildTree(const SparseVolume *sparseVolume,
                   unsigned int leftIdx, unsigned int rightIdx, int depth);

  SparseVolumeBuildOptions options_;
  std::vector<BlockNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G paritcles.
  SparseVolumeBuildStatistics stats_;
  const SparseVolume *sparseVolume_;

  double bmin_[3], bmax_[3]; // Bounding box of whole tree
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_ACCEL_PARTICLE_H__
