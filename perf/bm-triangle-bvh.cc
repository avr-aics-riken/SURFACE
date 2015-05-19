#include <hayai.hpp>

#include <cstdio>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "harness.h"

namespace {

int gSeed = 0;

}

//void BuildTriangleBVH(size_t n) {
//
//  double bmin[3] = {-1.0, -1.0, -1.0};
//  double bmax[3] = { 1.0,  1.0,  1.0};
//
//  float *pts = new float[n * 9];
//  GenerateRandomTrianglesFloat(pts, n, bmin, bmax, gSeed);
//
//  double elap = BM_BuildTriangleBVH(pts, n);
//  delete [] pts;
//
//}

BENCHMARK_P(BVH, BuildTriangle32BVH, 2, 3,
  (int numThreads, std::size_t numPoints))
{
  
  size_t n = numPoints;

#if _OPENMP
  omp_set_num_threads(numThreads);
#endif

  double bmin[3] = {-1.0, -1.0, -1.0};
  double bmax[3] = { 1.0,  1.0,  1.0};

  printf("Generating %lld random triangles...\n", n);
  float *pts = new float[n * 9];
  GenerateRandomTrianglesFloat(pts, n, bmin, bmax, gSeed);

  printf("Start building BVH .\n");
  double elap = BM_BuildTriangle32BVH(pts, n);
  printf("End building BVH . \n");

  int cores = 1;
#if _OPENMP
  cores = omp_get_max_threads();
#endif

  printf("{\"cores\": %d, \"n\" : %lld, \"secs\" : %f}\n", cores, n, elap);
  delete [] pts;

}

BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (1, 1024));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (1, 1024*1024));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (1, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (2, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (3, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (4, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (5, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (6, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (7, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (8, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (9, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (10, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (12, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (13, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (14, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (15, 1024ULL*1024ULL*10ULL));
BENCHMARK_P_INSTANCE(BVH, BuildTriangle32BVH, (16, 1024ULL*1024ULL*10ULL));
