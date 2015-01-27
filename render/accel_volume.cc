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

#include "accel_volume.h"

using namespace lsgl::render;

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

bool SparseVolumeAccel::Build(const SparseVolume *sparseVolume) {

  assert(sparseVolume);

  for (size_t i = 0; i < sparseVolume->blocks.size(); i++) {
    const VolumeBlock &block = sparseVolume->blocks[i];

    BVHData nodeData;
    nodeData.bmin[0] = (double)block.offset[0];
    nodeData.bmin[1] = (double)block.offset[1];
    nodeData.bmin[2] = (double)block.offset[2];
    nodeData.bmax[0] = (double)(block.offset[0] + block.extent[0]);
    nodeData.bmax[1] = (double)(block.offset[1] + block.extent[1]);
    nodeData.bmax[2] = (double)(block.offset[2] + block.extent[2]);
    nodeData.nodeID = i;

    tree_.AddNode(nodeData);
  }

  tree_.BuildTree();

  sparseVolume_ = sparseVolume; // retain pointer.

  return true;
}

bool SparseVolumeAccel::Sample(double value[4],
                               const double position[3]) const {

  std::vector<BVHNodeLocator> locaters;
  bool hit = tree_.Locate(locaters, position);

  if (hit) {
    // @note { We don't allow overlapping of volume block. }
    BVHNodeLocator locator = locaters[0];

    // @fixme
    value[0] = 1.0f;
    value[1] = 1.0f;
    value[2] = 1.0f;
    value[3] = 1.0f;
  } else {
    value[0] = 0.0f;
    value[1] = 0.0f;
    value[2] = 0.0f;
    value[3] = 0.0f;
  }

  return hit;
}

void SparseVolumeAccel::BoundingBox(double bmin[3], double bmax[3]) const {
  if (!sparseVolume_ || sparseVolume_->blocks.empty()) {
    bmin[0] = std::numeric_limits<double>::max();
    bmin[1] = std::numeric_limits<double>::max();
    bmin[2] = std::numeric_limits<double>::max();
    bmax[0] = -std::numeric_limits<double>::max();
    bmax[1] = -std::numeric_limits<double>::max();
    bmax[2] = -std::numeric_limits<double>::max();
  } else {
    tree_.BoundingBox(bmin, bmax);
  }
}
