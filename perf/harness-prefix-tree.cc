#include "harness.h"
#include "../render/render_prefix_tree_util.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <vector>

using namespace lsgl::render;

void BM_MortonCodes30(const float *points, size_t n, double bmin[3], double bmax[3])
{
  std::vector<uint32_t> codes(n);
  real3 bbox_min(bmin[0], bmin[1], bmin[2]);
  real3 bbox_max(bmax[0], bmax[1], bmax[2]);
  CalculateMortonCodes30(&codes.at(0), points, bbox_min, bbox_max, 0, n);
}

void BM_MortonCodes30SIMD(const float *points, size_t n, double bmin[3], double bmax[3])
{
  std::vector<uint32_t> codes(n);
  real3 bbox_min(bmin[0], bmin[1], bmin[2]);
  real3 bbox_max(bmax[0], bmax[1], bmax[2]);
  CalculateMortonCodes30SIMD(&codes.at(0), points, bbox_min, bbox_max, 0, n);
}

void BM_RadixSortMortonCodes30(void* arg)
{
  std::vector<IndexKey30>* keys = reinterpret_cast<std::vector<IndexKey30>* >(arg);

  RadixSortByMortionCode30MSB(&keys->at(0), &keys->at(0) + keys->size());
}

void BM_MergeSortMortonCodes30(void* arg)
{
  std::vector<IndexKey30>* keys = reinterpret_cast<std::vector<IndexKey30>* >(arg);

  std::vector<IndexKey30> tmp(keys->size());

  size_t n = keys->size();
  MergeSort30(&keys->at(0), &tmp.at(0), 0, n-1);
}

void BM_ParallelRadixSortMortonCodes30(void* arg)
{
  std::vector<IndexKey30>* keys = reinterpret_cast<std::vector<IndexKey30>* >(arg);

  std::vector<IndexKey30> tmp(keys->size());

  size_t n = keys->size();
  RadixSort30(&keys->at(0), &keys->at(0) + n);
}
