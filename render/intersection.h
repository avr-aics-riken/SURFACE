/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_INTERSECTION_H__
#define __LSGL_INTERSECTION_H__

#include "render_common.h"

namespace lsgl {
namespace render {

struct Intersection {
  float t;
  float u; // varycentric coord
  float v; // varycentric coord
  unsigned int prim_id;

  unsigned int subface_id; // For tetra primitive

  real3 position;
  real3 geometric; // geometric_normal
  real3 normal;    // shading normal
  real3 tangent;
  real3 binormal;

  unsigned int f0;
  unsigned int f1;
  unsigned int f2;
  unsigned int f3; // for tetra primitive

#ifdef ENABLE_TRAVERSE_PROFILING
  unsigned int numTriangleTests_;
  unsigned int numTraversals_;
#endif

  unsigned char *renderElement; // Opeque pointer to render element in the
                                // scene.

  void clear() {
    renderElement = NULL;
    prim_id = (unsigned int)(-1);
  }
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_INTERSECTION_H__
