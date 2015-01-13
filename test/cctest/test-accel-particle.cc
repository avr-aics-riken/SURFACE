#include "gtest/gtest.h"

#include "../render/accel_particle.h"

using namespace lsgl::render;

TEST(LSGLRenderTest, SameKeyPrefixParticles) {

  int ns[] = {8, 16, 32};

  int niter = sizeof(ns) / sizeof(int);
  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    float *points = new float[3*npoints];
    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = 0.5;
      points[3*j+1] = 0.5;
      points[3*j+2] = 0.5;
    }

    ParticleBuildOptions options;

    ParticleAccel *accel = new ParticleAccel();

    Particles part;

    part.radius = NULL;
    part.constantRadius = 0.01; // @fixme {}

    part.dpositions = NULL;
    part.positions = points;
    part.isDoublePrecisionPos = false;
    part.numParticles = npoints;

    accel->Build(&part, options);

    delete accel;
    delete [] points;
  }

}

TEST(LSGLRenderTest, FastAccelBuilder30) {

  int ns[] = {8, 16, 32};

  int niter = sizeof(ns) / sizeof(int);
  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    float *points = new float[3*npoints];
    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = 0.5;
      points[3*j+1] = 0.5;
      points[3*j+2] = 0.5;
    }

    ParticleBuildOptions options;

    ParticleAccel *accel = new ParticleAccel();

    Particles part;

    part.radius = NULL;
    part.constantRadius = 0.01; // @fixme {}

    part.dpositions = NULL;
    part.positions = points;
    part.isDoublePrecisionPos = false;
    part.numParticles = npoints;

    accel->Build32(&part, options);

    delete accel;
    delete [] points;
  }

}

