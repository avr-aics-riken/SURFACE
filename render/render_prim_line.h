/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LINE_H__
#define __LINE_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

//
// Line primitive is defined as 2 vertices.
//
//             +----+----+----+----+----+      +----+
// vertices    + v0 | v1 | v2 | v3 | v4 | .... | vN |
//             +----+----+----+----+----+      +----+
//
//             +--------+-------+      +---------+
// segmentss   | i0, i1 | i0, i1| .... | iN0, iN1|
//             +--------+-------+      +---------+
//
typedef struct {
  size_t numLines;
  const float *positions;   /// float precision. [xyz] * numLines * 2
  const double *dpositions; /// double precision. [xyz] * numLines * 2.
  bool isDoublePrecisionPos;
  const real *radius;  /// particle radius
  real constantRadius; /// Valid if radius == NULL
  uint32_t *segments;  /// Line segments(2 * numLines)
} Lines;

} // namespace render
} // namespace lsgl

#endif // __LINE_H__
