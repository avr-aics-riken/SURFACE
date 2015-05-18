//
// Generate random triangles
//
#include "harness.h"
#include "../render/tinymt64.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

void GenerateRandomTrianglesFloat(
  float* vertices, size_t n, double bmin[3], double bmax[3], int seed)
{

  assert(vertices);

  tinymt64_t rng;
  rng.mat1 = 0xfa051f40;
  rng.mat2 = 0xffd0fff4;
  rng.tmat = 0x58d02ffeffbfffbcULL;
  tinymt64_init(&rng, seed);

  float f = tinymt64_generate_double(&rng);

  //#ifdef _OPENMP
  //#pragma omp parallel for
  //#endif
  for (size_t i = 0; i < n; i++) {

    // [0, 1.0)
    double u0 = tinymt64_generate_double(&rng);
    double u1 = tinymt64_generate_double(&rng);
    double u2 = tinymt64_generate_double(&rng);

    // [bmin, bmax)
    double px = (bmax[0] - bmin[0]) * u0 + bmin[0];
    double py = (bmax[1] - bmin[1]) * u1 + bmin[1];
    double pz = (bmax[2] - bmin[2]) * u2 + bmin[2];

    for (size_t j = 0; j < 3; j++) {

      // random offset. [-0.025, 0.025)
      double v0 = 0.05 * tinymt64_generate_double(&rng) - 0.025;
      double v1 = 0.05 * tinymt64_generate_double(&rng) - 0.025;
      double v2 = 0.05 * tinymt64_generate_double(&rng) - 0.025;

      vertices[3*(3*i+j)+0] = px + v0;
      vertices[3*(3*i+j)+1] = py + v1;
      vertices[3*(3*i+j)+2] = pz + v2;

    }
  }
  
}
