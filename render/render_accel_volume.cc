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

#include "render_texture.h" // LSGL_RENDER_TEXTURE3D_FORMAT_XXXX
#include "render_accel_volume.h"

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

  StackVector<BVHNodeLocator, 32> locaters;
  bool hit = tree_.Locate(locaters, position);

  value[0] = 0.0f;
  value[1] = 0.0f;
  value[2] = 0.0f;
  value[3] = 1.0f; // In OpenGL, alpha is 1.0 by default.

  if (hit) {
    // @note { We don't allow overlapping of volume block. }
    BVHNodeLocator locator = locaters[0];

    // fetch texel color
    int blockID = locator.nodeID;
    int components =
        (sparseVolume_->components >= 4) ? 4 : sparseVolume_->components;
    int type = (sparseVolume_->components >= 4) ? 4 : sparseVolume_->components;
    const VolumeBlock &block = sparseVolume_->blocks[blockID];

    // find local offset.
    double x = position[0] - block.offset[0];
    double y = position[1] - block.offset[1];
    double z = position[2] - block.offset[2];
    int ix = (std::max)(
        (std::min)(block.offset[0] + block.extent[0] - 1, (int)x), 0);
    int iy = (std::max)(
        (std::min)(block.offset[1] + block.extent[1] - 1, (int)y), 1);
    int iz = (std::max)(
        (std::min)(block.offset[2] + block.extent[2] - 1, (int)z), 2);

    if (sparseVolume_->format == LSGL_RENDER_TEXTURE3D_FORMAT_BYTE) { // uint8
      const unsigned char *voxel = block.data;
      assert(voxel);

      size_t idx =
          iz * block.extent[0] * block.extent[1] + iy * block.extent[0] + ix;
      for (int c = 0; c < components; c++) {
        value[c] = voxel[components * idx + c] / 255.0f;
      }
    } else if (sparseVolume_->format == LSGL_RENDER_TEXTURE3D_FORMAT_FLOAT) {
      const float *voxel = reinterpret_cast<const float *>(block.data);
      assert(voxel);
      size_t idx =
          iz * block.extent[0] * block.extent[1] + iy * block.extent[0] + ix;
      for (int c = 0; c < components; c++) {
        value[c] = voxel[components * idx + c];
      }
    } else if (sparseVolume_->format == LSGL_RENDER_TEXTURE3D_FORMAT_DOUBLE) {
      const double *voxel = reinterpret_cast<const double *>(block.data);
      assert(voxel);
      size_t idx =
          iz * block.extent[0] * block.extent[1] + iy * block.extent[0] + ix;
      for (int c = 0; c < components; c++) {
        value[c] = (float)voxel[components * idx + c];
      }
    }
  }

  return hit;
}

void SparseVolumeAccel::BoundingBox(double bmin[3], double bmax[3]) const {
  if (!sparseVolume_ || sparseVolume_->blocks.empty()) {
    bmin[0] = (std::numeric_limits<double>::max)();
    bmin[1] = (std::numeric_limits<double>::max)();
    bmin[2] = (std::numeric_limits<double>::max)();
    bmax[0] = -(std::numeric_limits<double>::max)();
    bmax[1] = -(std::numeric_limits<double>::max)();
    bmax[2] = -(std::numeric_limits<double>::max)();
  } else {
    tree_.BoundingBox(bmin, bmax);
  }
}
