#include "harness.h"
#include "../render/accel_triangle.h"
#include "../render/timerutil.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

using namespace lsgl::render;

double BM_BuildTriangleBVH(const float *points, size_t n)
{
  BVHBuildOptions options;

  BVHAccel *accel = new BVHAccel();

  uint32_t *faces = new uint32_t[3*n];
  for (size_t i = 0; i < 3 * n; i++) {
    faces[i] = i;
  }

  Mesh mesh;

  mesh.dvertices = NULL;
  mesh.vertices = points;
  mesh.isDoublePrecisionPos = false;
  mesh.nvertices = n;
  mesh.faces = faces;
  mesh.nfaces =  n;

  double elap;

  {
    lsgl::timerutil t;
    t.start();

    accel->Build(&mesh, options);
    t.end();

    printf("[DBG] Mesh accel built time (%lld points) : %d msec\n", n, (int) t.msec());

    elap = 0.001 * t.msec();

    double bmin[3], bmax[3];
    accel->BoundingBox(bmin, bmax);

    printf("[DBG]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    printf("[DBG]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  }

  delete accel;

  mesh.vertices = 0; // remote points in outside of this function.

  return elap;

}

double BM_BuildTriangle32BVH(const float *points, size_t n)
{
  BVHBuildOptions options;

  BVHAccel *accel = new BVHAccel();

  uint32_t *faces = new uint32_t[3*n];
  for (size_t i = 0; i < 3 * n; i++) {
    faces[i] = i;
  }

  Mesh mesh;

  mesh.dvertices = NULL;
  mesh.vertices = points;
  mesh.isDoublePrecisionPos = false;
  mesh.nvertices = n;
  mesh.faces = faces;
  mesh.nfaces =  n;


  double elapsedSec = 0.0;

  {
    lsgl::timerutil t;
    t.start();

    accel->Build32(&mesh, options);
    t.end();

    printf("[DBG] Mesh accel built time (%lld points) : %d msec\n", n, (int) t.msec());

    elapsedSec = 0.001 * t.msec();

    double bmin[3], bmax[3];
    accel->BoundingBox(bmin, bmax);

    printf("[DBG]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    printf("[DBG]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  }

  delete accel;

  mesh.vertices = 0; // remote points in outside of this function.

  return elapsedSec;
}
