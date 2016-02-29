/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2016 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

// @todo { Move Tetra to Sold. }
#ifndef __LSGL_RENDER_PRIM_SOLID_H__
#define __LSGL_RENDER_PRIM_SOLID_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

//
// Pyramid primitive is defined as 5 vertices.
//
//             +----+----+----+----+----+      +----+
// vertices    + v0 | v1 | v2 | v3 | v4 | .... | vN |
//             +----+----+----+----+----+      +----+
//
//             +--------------------+--------------------+      +-------------------------+
// indices     | i0, i1, i2, i3, i4 | i0, i1, i2, i3, i4 | .... | iN0, iN1, iN2, iN3, iN4 |
//             +--------------------+--------------------+      +-------------------------+
//
// In similar fasion, Prism and Hexahedron are defined as 6 vertices and 8 vertices, respectively.
//
typedef struct {
  int    numVertsPerSolid;  /// 5 = Pyramid, 6 = Prism, 8 = Hexahedron.
  size_t numSolids;
  const float *vertices;   /// float precision. [xyz] * numSolids * numPlanes
  const double *dvertices; /// float precision. [xyz] * numSolids * numPlanes
  bool isDoublePrecisionPos;
  uint32_t *indices; /// primitive indices. up to 2^32 - 1 primitives.
} Solid;

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_PRIM_SOLID_H__
