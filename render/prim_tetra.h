/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_PRIM_TETRA_H__
#define __LSGL_RENDER_PRIM_TETRA_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

//
// Tetra primitive is defined as 4 vertices.
//
//             +----+----+----+----+----+      +----+
// vertices    + v0 | v1 | v2 | v3 | v4 | .... | vN |
//             +----+----+----+----+----+      +----+
//
//             +----------------+----------------+      +--------------------+
// faces       | i0, i1, i2, i3 | i0, i1, i2, i3 | .... | iN0, iN1, iN2, iN3 |
//             +----------------+----------------+      +--------------------+
//
typedef struct {
  size_t numTetrahedrons;
  const float *vertices;   /// float precision. [xyz] * numTetrahedrons
  const double *dvertices; /// float precision. [xyz] * numParticles
  bool isDoublePrecisionPos;
  uint32_t *faces; /// 1 tetrahedron face = 4 vertex indices. up to 2^32 - 1 primitives.
} Tetrahedron;

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_PRIM_TETRA_H__
