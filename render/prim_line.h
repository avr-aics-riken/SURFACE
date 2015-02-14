/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LINE_H__
#define __LINE_H__

#include <cstdio>
#include "render_common.h"

namespace lsgl {
namespace render {

typedef struct {
  size_t numLines;
  const float *positions;   /// float precision. [xyz] * numLines * 2
  const double *dpositions; /// double precision. [xyz] * numLines * 2.
  bool isDoublePrecisionPos;
  const real *radius;        /// particle radius
  real constantRadius; /// Valid if radius == NULL
} Lines;

} // namespace render
} // namespace lsgl

#endif // __LINE_H__
