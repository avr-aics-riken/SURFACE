/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_INTERSECTION_H__
#define __LSGL_INTERSECTION_H__

#include "vector3.h"

namespace lsgl {
namespace render {

struct Intersection {
  float t;
  float u; // varycentric coord
  float v; // varycentric coord
  unsigned int prim_id;

  vector3 position;
  vector3 geometric; // geometric_normal
  vector3 normal;    // shading normal

  unsigned int f0;
  unsigned int f1;
  unsigned int f2;

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
