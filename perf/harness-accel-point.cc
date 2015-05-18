//
// Generate random points
//
#include "harness.h"
#include "../render/render_accel_particle.h"
#include "../render/render_timerutil.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

using namespace lsgl::render;

double BM_BuildPointBVH(const float *points, size_t n)
{
  ParticleBuildOptions options;

  ParticleAccel *accel = new ParticleAccel();

  Particles part;

  part.radius = NULL;
  part.constantRadius = 0.01; // @fixme {}

  part.dpositions = NULL;
  part.positions = points;
  part.isDoublePrecisionPos = false;
  part.numParticles = n;

  {
    lsgl::timerutil t;
    t.start();

    accel->Build(&part, options);
    t.end();

    printf("[DBG] Particle accel built time (%lld points) : %d msec\n", n, (int) t.msec());

    double bmin[3], bmax[3];
    accel->BoundingBox(bmin, bmax);

    printf("[DBG]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    printf("[DBG]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  }

  delete accel;

}

double BM_BuildPoint32BVH(const float *points, size_t n)
{
  ParticleBuildOptions options;

  ParticleAccel *accel = new ParticleAccel();

  Particles part;

  part.radius = NULL;
  part.constantRadius = 0.01; // @fixme {}

  part.dpositions = NULL;
  part.positions = points;
  part.isDoublePrecisionPos = false;
  part.numParticles = n;

  double elapsedSec = 0.0;

  {
    lsgl::timerutil t;
    t.start();

    accel->Build32(&part, options);
    t.end();

    printf("[DBG] Particle accel built time (%lld points) : %d msec\n", n, (int) t.msec());

    elapsedSec = 0.001 * t.msec();

    double bmin[3], bmax[3];
    accel->BoundingBox(bmin, bmax);

    printf("[DBG]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    printf("[DBG]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  }

  delete accel;

  return elapsedSec;
}
