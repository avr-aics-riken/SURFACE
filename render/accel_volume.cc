#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <limits>
#include <functional>
#include <algorithm>

#include "accel_volume.h"

using namespace lsgl::render;

#define MAX_LEAF_ELEMENTS (8)
#define MAX_TREE_DEPTH_32BIT                                                   \
  (22) // FYI, log2(1G/16) ~= 25.897, log2(1G/32) ~= 21

#define ENABLE_TRACE_PRINT (0)
#define ENABLE_DEBUG_PRINT (0)

#define trace(f, ...)                                                          \
  {                                                                            \
    if (ENABLE_TRACE_PRINT)                                                    \
      printf(f, __VA_ARGS__);                                                  \
  }

#if ENABLE_DEBUG_PRINT
#define debug(f, ...)                                                          \
  { printf(f, __VA_ARGS__); }
#else
#define debug(f, ...)
#endif

namespace {

template <typename T>
static void ComputeBoundingBox(real3 &bmin, real3 &bmax, const T *vertices,
                               unsigned int *faces, unsigned int *indices,
                               unsigned int leftIndex,
                               unsigned int rightIndex) {
  const real kEPS = std::numeric_limits<real>::epsilon() * 1024;

  bmin[0] = std::numeric_limits<real>::max();
  bmin[1] = std::numeric_limits<real>::max();
  bmin[2] = std::numeric_limits<real>::max();
  bmax[0] = -std::numeric_limits<real>::max();
  bmax[1] = -std::numeric_limits<real>::max();
  bmax[2] = -std::numeric_limits<real>::max();

  if (rightIndex <= leftIndex) {
    return;
  }

  assert(0); // @todo {}

  size_t i = leftIndex;
  size_t idx = indices[i];
  bmin[0] = vertices[3 * faces[4 * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * faces[4 * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * faces[4 * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * faces[4 * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * faces[4 * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * faces[4 * idx + 0] + 2] + kEPS;

  for (i = leftIndex; i < rightIndex; i++) { // for each faces
    //size_t idx = indices[i];
    //for (int j = 0; j < 4; j++) { // for each tetra vertex
    //  size_t fid = faces[4 * idx + j];
    //  for (int k = 0; k < 3; k++) { // xyz
    //    real minval = vertices[3 * fid + k] - kEPS;
    //    real maxval = vertices[3 * fid + k] + kEPS;
    //    if (bmin[k] > minval)
    //      bmin[k] = minval;
    //    if (bmax[k] < maxval)
    //      bmax[k] = maxval;
    //  }
    //}
  }
}

} // namespace

//
// --
//

size_t SparseVolumeAccel::BuildTree(const SparseVolume *sparseVolume, unsigned int leftIdx,
                             unsigned int rightIdx, int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  assert(sparseVolume);

  size_t offset = nodes_.size();

  if (stats_.maxTreeDepth < depth) {
    stats_.maxTreeDepth = depth;
  }

  real3 bmin, bmax;
  assert(0);
  // @todo { bbox }
  //ComputeBoundingBox<float>(bmin, bmax, tetras->vertices, tetras->faces,
  //                            &indices_.at(0), leftIdx, rightIdx);

  debug(" bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  debug(" bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {
    // Create leaf node.
	assert(0);
	// @todo { bbox }
    //TetraNode leaf;

    //leaf.bmin[0] = bmin[0];
    //leaf.bmin[1] = bmin[1];
    //leaf.bmin[2] = bmin[2];

    //leaf.bmax[0] = bmax[0];
    //leaf.bmax[1] = bmax[1];
    //leaf.bmax[2] = bmax[2];

    //assert(leftIdx < std::numeric_limits<unsigned int>::max());

    //leaf.flag = 1; // leaf
    //leaf.data[0] = n;
    //leaf.data[1] = (unsigned int)leftIdx;

    //nodes_.push_back(leaf);

    //stats_.numLeafNodes++;

    return offset;
  }

  //
  // Create branch node.
  //

  assert(0); // @todo
  //TetraNode node;
  //node.axis = cutAxis;
  //node.flag = 0; // 0 = branch
  //nodes_.push_back(node);

  // Recurively split tree.
  unsigned int midIdx = 0; // @fixme
  unsigned int leftChildIndex = BuildTree(sparseVolume, leftIdx, midIdx, depth + 1);
  unsigned int rightChildIndex = BuildTree(sparseVolume, midIdx, rightIdx, depth + 1);

  nodes_[offset].data[0] = leftChildIndex;
  nodes_[offset].data[1] = rightChildIndex;

  //nodes_[offset].bmin[0] = bmin[0];
  //nodes_[offset].bmin[1] = bmin[1];
  //nodes_[offset].bmin[2] = bmin[2];

  //nodes_[offset].bmax[0] = bmax[0];
  //nodes_[offset].bmax[1] = bmax[1];
  //nodes_[offset].bmax[2] = bmax[2];

  stats_.numBranchNodes++;

  return offset;
}

bool SparseVolumeAccel::Build(const SparseVolume *sparseVolume,
                       const SparseVolumeBuildOptions &options) {
  options_ = options;
  stats_ = SparseVolumeBuildStatistics();

  assert(sparseVolume);

  size_t n = sparseVolume->blocks.size();
  trace("[SparseVolumeAccel] Input # of volume blocks    = %lu\n", n);

  //
  // 1. Create indices(this will be permutated in BuildTree)
  //
  indices_.resize(n);
  for (size_t i = 0; i < n; i++) {
    indices_[i] = i;
  }

  //
  // 2. Build tree
  //
  BuildTree(sparseVolume, 0, n, 0);

  // Tree will be null if input triangle count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    real3 bmin(nodes_[0].bmin[0][0], nodes_[0].bmin[0][1], nodes_[0].bmin[0][2]);
    real3 bmax(nodes_[0].bmax[0][0], nodes_[0].bmax[0][1], nodes_[1].bmax[0][2]);
    trace("[SparseVolumeAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
    trace("[SparseVolumeAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);
  }

  trace("[SparseVolumeAccel] # of nodes = %lu\n", nodes_.size());

  // Store pointer for later use.
  sparseVolume_ = sparseVolume;

  return true;
}

bool SparseVolumeAccel::Dump(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "[SparseVolumeAccel] Cannot write a file: %s\n", filename);
    return false;
  }

  assert(0); // @todo

  //unsigned long long numNodes = nodes_.size();
  //assert(nodes_.size() > 0);

  //unsigned long long numIndices = indices_.size();

  //int r = 0;
  //r = fwrite(&numNodes, sizeof(unsigned long long), 1, fp);
  //assert(r == 1);

  //r = fwrite(&nodes_.at(0), sizeof(TetraNode), numNodes, fp);
  //assert(r == numNodes);

  //r = fwrite(&numIndices, sizeof(unsigned long long), 1, fp);
  //assert(r == 1);

  //r = fwrite(&indices_.at(0), sizeof(unsigned int), numIndices, fp);
  //assert(r == numIndices);

  //fclose(fp);

  return false;
}

namespace {

const int kMaxStackDepth = 512;

inline bool IntersectRayAABB(real &tminOut, // [out]
                             real &tmaxOut, // [out]
                             real maxT, real const bmin[3], real const bmax[3],
                             real3 rayOrg, real3 rayInvDir, int rayDirSign[3]) {
  real tmin, tmax;

  const real min_x = rayDirSign[0] ? bmax[0] : bmin[0];
  const real min_y = rayDirSign[1] ? bmax[1] : bmin[1];
  const real min_z = rayDirSign[2] ? bmax[2] : bmin[2];
  const real max_x = rayDirSign[0] ? bmin[0] : bmax[0];
  const real max_y = rayDirSign[1] ? bmin[1] : bmax[1];
  const real max_z = rayDirSign[2] ? bmin[2] : bmax[2];

  // X
  const double tmin_x = (min_x - rayOrg[0]) * rayInvDir[0];
  const double tmax_x = (max_x - rayOrg[0]) * rayInvDir[0];

  // Y
  const double tmin_y = (min_y - rayOrg[1]) * rayInvDir[1];
  const double tmax_y = (max_y - rayOrg[1]) * rayInvDir[1];

  tmin = (tmin_x > tmin_y) ? tmin_x : tmin_y;
  tmax = (tmax_x < tmax_y) ? tmax_x : tmax_y;

  // Z
  const double tmin_z = (min_z - rayOrg[2]) * rayInvDir[2];
  const double tmax_z = (max_z - rayOrg[2]) * rayInvDir[2];

  tmin = (tmin > tmin_z) ? tmin : tmin_z;
  tmax = (tmax < tmax_z) ? tmax : tmax_z;

  //
  // Hit include (tmin == tmax) edge case(hit 2D plane).
  //
  if ((tmax > 0.0) && (tmin <= tmax) && (tmin <= maxT)) {

    // printf("tmin, tmax = %f, %f\n", tmin, tmax);
    tminOut = tmin;
    tmaxOut = tmax;

    return true;
  }

  return false; // no hit
}

} // namespace

bool SparseVolumeAccel::Traverse(Intersection &isect, Ray &ray) const {
  real hitT = std::numeric_limits<real>::max(); // far = no hit.

  int nodeStackIndex = 0;
  std::vector<int> nodeStack(512);
  nodeStack[0] = 0;

  // Init isect info as no hit
  isect.t = hitT;
  isect.u = 0.0;
  isect.v = 0.0;
  isect.prim_id = (unsigned int)(-1);

  int dirSign[3];
  dirSign[0] = ray.direction()[0] < 0.0 ? 1 : 0;
  dirSign[1] = ray.direction()[1] < 0.0 ? 1 : 0;
  dirSign[2] = ray.direction()[2] < 0.0 ? 1 : 0;

  // @fixme { Check edge case; i.e., 1/0 }
  real3 rayInvDir;
  rayInvDir[0] = 1.0 / ray.direction()[0];
  rayInvDir[1] = 1.0 / ray.direction()[1];
  rayInvDir[2] = 1.0 / ray.direction()[2];

  real3 rayOrg;
  rayOrg[0] = ray.origin()[0];
  rayOrg[1] = ray.origin()[1];
  rayOrg[2] = ray.origin()[2];

  real minT, maxT;
  while (nodeStackIndex >= 0) {
    int index = nodeStack[nodeStackIndex];
	assert(0); // @todo
    //const TetraNode &node = nodes_[index];

    //nodeStackIndex--;

    //if (node.flag == 0) { // branch node

    //  bool hit = IntersectRayAABB(minT, maxT, hitT, node.bmin, node.bmax,
    //                              rayOrg, rayInvDir, dirSign);

    //  if (hit) {

    //    int orderNear = dirSign[node.axis];
    //    int orderFar = 1 - orderNear;

    //    // Traverse near first.
    //    nodeStack[++nodeStackIndex] = node.data[orderFar];
    //    nodeStack[++nodeStackIndex] = node.data[orderNear];
    //  }

    //} else { // leaf node

    //  //if (TestLeafNode(isect, node, indices_, tetras_, ray)) {
    //  //  hitT = isect.t;
    //  //}
    //}
  }

  assert(nodeStackIndex < kMaxStackDepth);

  if (isect.t < std::numeric_limits<real>::max()) {
    //BuildIntersection(isect, tetras_, ray);
    return true;
  }

  return false;
}

void SparseVolumeAccel::BoundingBox(double bmin[3], double bmax[3]) const {
  if (nodes_.empty()) {
    bmin[0] = std::numeric_limits<double>::max();
    bmin[1] = std::numeric_limits<double>::max();
    bmin[2] = std::numeric_limits<double>::max();
    bmax[0] = -std::numeric_limits<double>::max();
    bmax[1] = -std::numeric_limits<double>::max();
    bmax[2] = -std::numeric_limits<double>::max();
  } else {
    bmin[0] = nodes_[0].bmin[0][0];
    bmin[1] = nodes_[0].bmin[0][1];
    bmin[2] = nodes_[0].bmin[0][2];
    bmax[0] = nodes_[0].bmax[0][0];
    bmax[1] = nodes_[0].bmax[0][1];
    bmax[2] = nodes_[0].bmax[0][2];
  }
}
