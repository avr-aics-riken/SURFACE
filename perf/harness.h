//
// Benchmark harness
//

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

void GenerateRandomPointsFloat(
  float* points, size_t n, double bmin[3], double bmax[3], int seed);
void GenerateRandomTrianglesFloat(
  float* points, size_t n, double bmin[3], double bmax[3], int seed);

extern double BM_BuildTriangleBVH(const float *vertices, size_t n);
extern double BM_BuildTriangle32BVH(const float *vertices, size_t n);
extern double BM_BuildPointBVH(const float *points, size_t n);
extern double BM_BuildPoint32BVH(const float *points, size_t n);
extern void   BM_MortonCodes30(const float *points, size_t n, double bmin[3], double bmax[3]);
extern void   BM_MortonCodes30SIMD(const float *points, size_t n, double bmin[3], double bmax[3]);
extern void   BM_RadixSortMortonCodes30(void* keys); // std::vector<IndexKey30>
extern void   BM_MergeSortMortonCodes30(void* keys); // std::vector<IndexKey30>
extern void   BM_ParallelRadixSortMortonCodes30(void* keys); // std::vector<IndexKey30>

extern char* log_json(const char* title, const char* variant, int64_t cpu_ns);

#ifdef __cplusplus
}
#endif

