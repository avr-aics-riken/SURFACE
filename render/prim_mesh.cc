/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "render_common.h"
#include "prim_mesh.h"
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265
#endif

using namespace lsgl::render;

static float radians(float x) { return x * M_PI / 180.0f; }

static inline real3 select_normal(const real3 &nv, const real3 &nt,
                                    float limit_cos) {
  float lv = nv.length()*nv.length();
  float lt = nt.length()*nt.length();
  int n = 0;
  if (lv != 0)
    n |= 1;
  if (lt != 0)
    n |= 2;

  switch (n) {
  case 1:
    return nv;
  case 2:
    return nt;
  default: {
    if (dot(nt, nv) <= limit_cos) {
      return nt;
    } else {
      return nv;
    }
  }
  }
}

bool Mesh::CalculateNormals(float limit_angle) {
  if (nvertices <= 0)
    return true;
  if (nfaces <= 0)
    return true;
  std::vector<real3> nv(nvertices);
  std::vector<real3> nt(nfaces);
  for (size_t i = 0; i < nvertices; i++) {
    nv[i] = real3(0, 0, 0);
  }
  for (size_t i = 0; i < nfaces; i++) {
    unsigned int i0 = faces[3 * i + 0];
    unsigned int i1 = faces[3 * i + 1];
    unsigned int i2 = faces[3 * i + 2];
    real3 p0(vertices[3*i0+0], vertices[3*i0+1], vertices[3*i0+2]);
    real3 p1(vertices[3*i1+0], vertices[3*i1+1], vertices[3*i1+2]);
    real3 p2(vertices[3*i2+0], vertices[3*i2+1], vertices[3*i2+2]);
    real3 e1 = p1 - p0;
    real3 e2 = p2 - p0;
    real3 n = cross(e1, e2);

    nv[i0] = nv[i0] + n;
    nv[i1] = nv[i1] + n;
    nv[i2] = nv[i2] + n;

	real3 nn = n;
    nt[i] = nn.normalize();
  }

  for (size_t i = 0; i < nvertices; i++) {
    nv[i].normalize();
  }

  /*
                isFacevaryingnormal = true;
                if(normals){
                        delete[] normals;
                }
                normals = new float[3*3*nfaces];
                */

  if (0) { // facevaryign normal
    float limit_cos = cos(radians(limit_angle));

    if (!normals) {
      normals = new float[3 * 3 * nfaces];
    }
    for (size_t i = 0; i < nfaces; i++) {
      real3 n = nt[i];
      unsigned int i0 = faces[3 * i + 0];
      unsigned int i1 = faces[3 * i + 1];
      unsigned int i2 = faces[3 * i + 1];

      real3 n0 = nv[i0];
      real3 n1 = nv[i1];
      real3 n2 = nv[i2];

      real3 o0 = select_normal(n0, n, limit_cos);
      real3 o1 = select_normal(n1, n, limit_cos);
      real3 o2 = select_normal(n2, n, limit_cos);

      /*
                                real3 o0 = n;
                                real3 o1 = n;
                                real3 o2 = n;
                                */

      normals[9 * i + 0] = o0[0];
      normals[9 * i + 1] = o0[1];
      normals[9 * i + 2] = o0[2];

      normals[9 * i + 3] = o1[0];
      normals[9 * i + 4] = o1[1];
      normals[9 * i + 5] = o1[2];

      normals[9 * i + 6] = o2[0];
      normals[9 * i + 7] = o2[1];
      normals[9 * i + 8] = o2[2];
    }
  } else {
    if (!normals) {
      normals = new float[3 * nvertices];
    }
    /*
                        for(size_t i=0;i<nvertices;i++){
                                nv[i] =
       normalize(real3(rand(),rand(),rand()));
                        }
                        */
    for (size_t i = 0; i < nvertices; i++) {
      normals[3 * i + 0] = nv[i][0];
      normals[3 * i + 1] = nv[i][1];
      normals[3 * i + 2] = nv[i][2];
    }
    /*
                        for(size_t i=0;i<100;i++){
                                printf("%f %f
       %f\n",normals[3*i+0],normals[3*i+1],normals[3*i+2]);
                        }
                        */
  }
  return true;
}
