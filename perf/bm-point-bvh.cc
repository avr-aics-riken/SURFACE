#include <cstdio>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "harness.h"

namespace {

int gSeed = 0;

void BuildPointBVH(size_t n) {

  double bmin[3] = {-1.0, -1.0, -1.0};
  double bmax[3] = { 1.0,  1.0,  1.0};

  float *pts = new float[n * 3];
  GenerateRandomPointsFloat(pts, n, bmin, bmax, gSeed);

  double elap = BM_BuildPointBVH(pts, n);
  delete [] pts;

}

void BuildPoint32BVH(size_t n) {

  double bmin[3] = {-1.0, -1.0, -1.0};
  double bmax[3] = { 1.0,  1.0,  1.0};

  printf("Generating %lld random points...\n", n);
  float *pts = new float[n * 3];
  GenerateRandomPointsFloat(pts, n, bmin, bmax, gSeed);

  printf("Start building particle BVH .\n");
  double elap = BM_BuildPoint32BVH(pts, n);
  printf("End building particle BVH . \n", elap);

  int cores = 1;
#if _OPENMP
  cores = omp_get_max_threads();
#endif

  printf("{\"cores\": %d, \"n\" : %lld, \"secs\" : %f}\n", cores, n, elap);
  delete [] pts;

}

}

int main(int argc, const char **argv) {

  gSeed = 123456;

  if (argc < 2) {
    printf("Usage: %s numPoints(in kilo)\n", argv[0]);
  }

  size_t n = 100*1024*1024;

  if (argc > 1) {
    n = atoi(argv[1]) * 1024; // kilo
  }

#ifdef _OPENMP
  printf("OMP Enabled: %d threads\n", omp_get_max_threads());
#endif

  BuildPointBVH(n);
  //BuildPoint32BVH(n);

  return EXIT_SUCCESS;

}
