/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_PRIMITIVE_MESH_H__
#define __LSGL_PRIMITIVE_MESH_H__

#include "render_ray.h"
#include "render_intersection.h"

namespace lsgl {
namespace render {

///
/// Mesh class
///
class Mesh {
public:
  Mesh()
      : nvertices(0), vertices(0), dvertices(0), normals(0), nfaces(0),
        faces(0), isDoublePrecisionPos(false) {}

  ~Mesh() { Release(); }

  size_t nvertices;
  const float *vertices;
  float *normals;

  // double precision
  const double *dvertices;
  bool isDoublePrecisionPos;

  size_t nfaces;
  unsigned int *faces;

  bool CalculateNormals(float limit_angle = 30);

  inline void Release() {
    if (vertices)
      delete[] vertices;
    if (dvertices)
      delete[] dvertices;
    if (normals)
      delete[] normals;
    if (faces)
      delete[] faces;

    nvertices = 0;
    nfaces = 0;
  }
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_MESH_H__
