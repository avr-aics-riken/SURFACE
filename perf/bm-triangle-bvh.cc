#include <cstdio>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "harness.h"

namespace {

int gSeed = 0;

void BuildTriangleBVH(size_t n) {

  double bmin[3] = {-1.0, -1.0, -1.0};
  double bmax[3] = { 1.0,  1.0,  1.0};

  float *pts = new float[n * 9];
  GenerateRandomTrianglesFloat(pts, n, bmin, bmax, gSeed);

  double elap = BM_BuildTriangleBVH(pts, n);
  delete [] pts;

}

void BuildTriangle32BVH(size_t n) {

  double bmin[3] = {-1.0, -1.0, -1.0};
  double bmax[3] = { 1.0,  1.0,  1.0};

  printf("Generating %lld random triangles...\n", n);
  float *pts = new float[n * 9];
  GenerateRandomTrianglesFloat(pts, n, bmin, bmax, gSeed);

  printf("Start building BVH .\n");
  double elap = BM_BuildTriangle32BVH(pts, n);
  printf("End building BVH . \n", elap);

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
    printf("Usage: %s numTris(in kilo)\n", argv[0]);
  }

  size_t n = 10*1024;

  if (argc > 1) {
    n = atoi(argv[1]) * 1024; // kilo
  }

#ifdef _OPENMP
  printf("OMP Enabled: %d threads\n", omp_get_max_threads());
#endif

  //BuildTriangleBVH(n);
  BuildTriangle32BVH(n);

  return EXIT_SUCCESS;

}
