//
// Sparse volume data structure.
//
#ifndef __LSGL_RENDER_ACCEL_VOLUME_H__
#define __LSGL_RENDER_ACCEL_VOLUME_H__

#include <vector>

#include "render_common.h"
#include "render_prim_volume.h"
#include "render_bvh_tree.h"
#include "render_ray.h"
#include "render_intersection.h"

namespace lsgl {
namespace render {

//
// We use BVH approach to look up volume data for sparse volume.
// `BlockNode` represents non-empty region.
// It is possible that `BlockNode` can be overlap each other, but
// we generally assume each `BlockNode` does not overlap.
//
// @todo { Consider DDA style ray traversal }
//

class BlockNode {
public:
  BlockNode(){};
  ~BlockNode(){};

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

/// SparseVolume BVH acceleration class.
class SparseVolumeAccel {
public:
  SparseVolumeAccel() : sparseVolume_(0){};
  ~SparseVolumeAccel(){};

  /// Build volume BVH for input volume blocks.
  bool Build(const SparseVolume *sparseVolume);

  /// Sample sparse volume at specified position.
  /// Volume elements are limited to 4 elements(RGBA)
  /// Note that input range for `position` must within [bmin, bmax]^3
  bool Sample(double value[4], const double position[3]) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]) const;

private:
  const SparseVolume *sparseVolume_;
  BVHTree tree_;

  double bmin_[3], bmax_[3]; // Bounding box of whole tree
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_ACCEL_PARTICLE_H__
