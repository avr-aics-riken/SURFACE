#ifndef __LSGL_RENDER_PRIM_VOLUME_H__
#define __LSGL_RENDER_PRIM_VOLUME_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

/// \struct 
/// Volume block structure.
///
/// \note
/// Volume data is reprented in the integer coordinate.
/// Global transform(translation, rotation and scaling) is done in the scene graph level.
///
/// \note
/// offset + extent must be less than SparseVolume.globalDim
typedef struct {
  int offset[3];              ///< Block offset in world space.
  int extent[3];              ///< Block extent in world space(a multiple of size[]).
  int size[3];                ///< Cell size(Actual volume size)
  int id;                     ///< Block ID
  int level;                  ///< LoD level
  const unsigned char *data;  ///< Pointer to actual volume data(voxel data)
} VolumeBlock;

/// \struct 
/// Represents sparse volume primitive(contains a list of volume blocks)
typedef struct {
  std::vector<VolumeBlock> blocks;
  int globalDim[3]; // global dim(the extent containing all volume blocks).
  int components;   // # of components in each voxel(up to 4).
                    // Assume all volume blocks have same # of components.
  int format;       // Voxel format(LSGL_RENDER_TEXTURE3D_FORMAT_BYTE, ...)
                    // (Assume all volume block have same format)
  int maxLevel;     // Maximum LoD level in this volume.
} SparseVolume;

// NonUniformVolume.spacing.size() must be equal to NonUniformVolume.size.
typedef struct {
  std::vector<float> spacingX;
  std::vector<float> spacingY;
  std::vector<float> spacingZ;
  int size[3];
  int components;
  int type;
  unsigned char *data; // Opaque pointer representation of volume data.
} NonUniformVolume;

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_PRIM_VOLUME_H__
