/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_PARTICLE_H__
#define __LSGL_RENDER_PARTICLE_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

typedef struct {
  size_t numParticles;
  const float *positions;   /// float precision. [xyz] * numParticles
  const double *dpositions; /// float precision. [xyz] * numParticles
  bool isDoublePrecisionPos;
  real *radius;        /// particle radius
  int *materialIDs;    /// particle mat ID
  real constantRadius; /// Valid if radius == NULL
} Particles;

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_PARTICLE_H__
