/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
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
  unsigned int f3;
  unsigned int f4;
  unsigned int f5;
  unsigned int f6;
  unsigned int f7; // for tetra and solid primitive
  
  float d0;
  float d1;
  float d2;
  float d3;
  float d4;
  float d5;
  float d6;
  float d7;
  
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
