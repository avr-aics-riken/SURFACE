#ifndef __LSGL_RENDER_PRIM_VOLUME_H__
#define __LSGL_RENDER_PRIM_VOLUME_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

// offset + extent must be less than Volume.size
typedef struct {
	int offset[3];
	int extent[3];
	int id;
} VolumeBlock;

typedef struct {
  std::vector<VolumeBlock> blocks;
  int size[3];
  int components;
  int type;
  unsigned char* data;	// Opaque pointer representation of volume data.
} SparseVolume;

// NonUniformVolume.spacing.size() must be equal to NonUniformVolume.size.
typedef struct {
  std::vector<float> spacingX;
  std::vector<float> spacingY;
  std::vector<float> spacingZ;
  int size[3];
  int components;
  int type;
  unsigned char* data;	// Opaque pointer representation of volume data.
} NonUniformVolume;

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_PRIM_VOLUME_H__
