/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdint.h>

#include "render_bvh_tree.h"
#include "render_matrix.h"

using namespace lsgl::render;

#define TREE_MAXDEPTH (256)

static const double kEPS = 2.0e-15;

typedef struct {
  int sign[3];
  double org[3];
  double invdir[3];
  double dir[3];
} Ray;

// Sorter
static bool sorter_isect(const BVHNodeIntersection &a,
                         const BVHNodeIntersection &b) {
  return a.tMin < b.tMin;
}

static void SetupRay(Ray &ray) {
  ray.sign[0] = (ray.dir[0] < 0.0) ? 1 : 0;
  ray.sign[1] = (ray.dir[1] < 0.0) ? 1 : 0;
  ray.sign[2] = (ray.dir[2] < 0.0) ? 1 : 0;

  if (fabs(ray.dir[0]) > kEPS) {
    ray.invdir[0] = 1.0 / ray.dir[0];
  } else {
    if (ray.dir[0] < 0.0)
      ray.invdir[0] = -FLT_MAX;
    else
      ray.invdir[0] = FLT_MAX;
  }

  if (fabs(ray.dir[1]) > kEPS) {
    ray.invdir[1] = 1.0 / ray.dir[1];
  } else {
    if (ray.dir[1] < 0.0)
      ray.invdir[1] = -FLT_MAX;
    else
      ray.invdir[1] = FLT_MAX;
  }

  if (fabs(ray.dir[2]) > kEPS) {
    ray.invdir[2] = 1.0 / ray.dir[2];
  } else {
    if (ray.dir[2] < 0.0)
      ray.invdir[2] = -FLT_MAX;
    else
      ray.invdir[2] = FLT_MAX;
  }
}

static inline void SwapNode(std::vector<BVHData> &nodes, uint64_t i,
                            uint64_t j) {
  BVHData tmp;

  tmp = nodes[i];
  nodes[i] = nodes[j];
  nodes[j] = tmp;
}

template <typename T>
static size_t PartitionNodes(std::vector<T> &nodes, size_t indexLeft,
                             size_t indexRight, double splitPos,
                             size_t splitAxis) {
  assert(indexLeft < indexRight);

  size_t i;

  // Save left most element.
  T leftmost = nodes[indexLeft];

  size_t pivotIndex = indexLeft;

  for (i = indexLeft + 1; i < indexRight; i++) {

    //
    // An idea for further optiomisation: Calculate only bounding min and
    // max for axis 'split_axis'
    //
    if (nodes[i].bmax[splitAxis] < splitPos) {

      pivotIndex++;

      SwapNode(nodes, pivotIndex, i);

      assert(pivotIndex <= indexRight);
    }
  }

  assert(pivotIndex <= indexRight);

  nodes[indexLeft] = nodes[pivotIndex];
  nodes[pivotIndex] = leftmost;

  //
  // 2. Seek the partion index which satisfies
  //      nodes[i]->bbox.bmax[splitAxis] < splitPos &&
  //      nodes[i]->bbox.bmax[splitAxis+1] >= splitPos.
  //
  size_t splitIndex = indexLeft;
  for (i = indexLeft; i < indexRight; i++) {

    if (nodes[i].bmax[splitAxis] < splitPos) {
      splitIndex++;
    }
  }

  return splitIndex;
}

template <typename T>
static void CalculateBoundingBox(double *bmin, double *bmax,
                                 std::vector<T> &nodes, size_t startIndex,
                                 size_t endIndex) {
  assert(startIndex != endIndex);

  bmin[0] = nodes[startIndex].bmin[0];
  bmin[1] = nodes[startIndex].bmin[1];
  bmin[2] = nodes[startIndex].bmin[2];

  bmax[0] = nodes[startIndex].bmax[0];
  bmax[1] = nodes[startIndex].bmax[1];
  bmax[2] = nodes[startIndex].bmax[2];

  for (size_t i = startIndex + 1; i < endIndex; ++i) {
    for (int j = 0; j < 3; ++j) {
      bmin[j] = std::min(bmin[j], nodes[i].bmin[j]);
      bmax[j] = std::max(bmax[j], nodes[i].bmax[j]);
    }
  }
}

template <typename T>
static void CalcSceneBBox(double *bmin, double *bmax, std::vector<T> &nodes) {
  assert(nodes.size() > 0);

  bmin[0] = nodes[0].bmin[0];
  bmin[1] = nodes[0].bmin[1];
  bmin[2] = nodes[0].bmin[2];
  bmax[0] = nodes[0].bmax[0];
  bmax[1] = nodes[0].bmax[1];
  bmax[2] = nodes[0].bmax[2];

  for (size_t i = 1; i < nodes.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      bmin[j] = std::min(bmin[j], nodes[i].bmin[j]);
      bmax[j] = std::max(bmax[j], nodes[i].bmax[j]);
    }
  }
}

static inline bool TestPointAABB(const double bmin[3], const double bmax[3],
                                 const double position[3]) {

  if ((position[0] >= bmin[0]) && (position[0] <= bmax[0]) &&
      (position[1] >= bmin[1]) && (position[1] <= bmax[1]) &&
      (position[2] >= bmin[2]) && (position[2] <= bmax[2])) {
    return true;
  }

  return false;
}

static inline int TestPointBVHNode(const std::vector<BVHNode> &nodesTree,
                                   const BVHNode *node,
                                   const double position[3]) {

  int hitLeft = 0;
  int hitRight = 0;

  if (node->child[0] > 0) {

    if (TestPointAABB(nodesTree[node->child[0]].bmin,
                      nodesTree[node->child[0]].bmax, position)) {
      hitLeft = 1;
    }
  }

  if (node->child[1] > 0) {

    if (TestPointAABB(nodesTree[node->child[1]].bmin,
                      nodesTree[node->child[1]].bmax, position)) {
      hitRight = 1;
    }
  }

  return (hitRight << 1) | hitLeft;
}

static int TestRayAABB(double &tMinOut, double &tMaxOut, const double bmin[3],
                       const double bmax[3], Ray &ray) {
  const double minX = ray.sign[0] ? bmax[0] : bmin[0]; // 1 op
  const double minY = ray.sign[1] ? bmax[1] : bmin[1]; // 1 op
  const double minZ = ray.sign[2] ? bmax[2] : bmin[2]; // 1 op

  const double maxX = ray.sign[0] ? bmin[0] : bmax[0]; // 1 op
  const double maxY = ray.sign[1] ? bmin[1] : bmax[1]; // 1 op
  const double maxZ = ray.sign[2] ? bmin[2] : bmax[2]; // 1 op

  // X
  const double tMinX = (minX - ray.org[0]) * ray.invdir[0]; // 2 op
  const double tMaxX = (maxX - ray.org[0]) * ray.invdir[0]; // 2 op

  // Y
  const double tMinY = (minY - ray.org[1]) * ray.invdir[1]; // 2 op
  const double tMaxY = (maxY - ray.org[1]) * ray.invdir[1]; // 2 op

  double tMin = (tMinX > tMinY) ? tMinX : tMinY; // 2 op
  double tMax = (tMaxX < tMaxY) ? tMaxX : tMaxY; // 2 op

  // Z
  const double tMinZ = (minZ - ray.org[2]) * ray.invdir[2]; // 2 op
  const double tMaxZ = (maxZ - ray.org[2]) * ray.invdir[2]; // 2 op

  tMin = (tMin > tMinZ) ? tMin : tMinZ; // 2 op
  tMax = (tMax < tMaxZ) ? tMax : tMaxZ; // 2 op

  if ((tMax > 0.0) && (tMin <= tMax)) { // 3 op
                                        // printf("ray.org = %f, %f, %f\n",
                                        //    ray.org[0],
                                        //    ray.org[1],
                                        //    ray.org[2]);
                                        // printf("ray.dir = %f, %f, %f\n",
                                        //    ray.dir[0],
                                        //    ray.dir[1],
                                        //    ray.dir[2]);
                                        // printf("bmin = %f, %f, %f\n",
                                        //    bmin[0],
                                        //    bmin[1],
                                        //    bmin[2]);
                                        // printf("bmax = %f, %f, %f\n",
                                        //    bmax[0],
                                        //    bmax[1],
                                        //    bmax[2]);
    tMinOut = tMin;
    tMaxOut = tMax;

    // printf("hit\n");

    return 1;
  }

  return 0;
}

static inline int TestRayNodeIntersection(double *tMinOut, // [2]
                                          double *tMaxOut, // [2]
                                          const std::vector<BVHNode> &nodesTree,
                                          const BVHNode *bboxNode, Ray &ray) {
  int hitLeft = 0;
  int hitRight = 0;

  if (bboxNode->child[0] > 0) {

#if 0 // DBG
    printf("bbox left = %f, %f, %f - %f, %f, %f\n",
           nodesTree[bboxNode->child[0]].bmin[0],
           nodesTree[bboxNode->child[0]].bmin[1],
           nodesTree[bboxNode->child[0]].bmin[2],
           nodesTree[bboxNode->child[0]].bmax[0],
           nodesTree[bboxNode->child[0]].bmax[1],
           nodesTree[bboxNode->child[0]].bmax[2]);

    printf("ray.org = %f, %f, %f\n", ray.org[0], ray.org[1], ray.org[2]);
    printf("ray.dir = %f, %f, %f\n", ray.dir[0], ray.dir[1], ray.dir[2]);
#endif

    if (TestRayAABB(tMinOut[0], tMaxOut[0], nodesTree[bboxNode->child[0]].bmin,
                    nodesTree[bboxNode->child[0]].bmax, ray)) {
      hitLeft = 1;
    }
  }

  if (bboxNode->child[1] > 0) {

#if 0
    printf("bbox right  = %f, %f, %f - %f, %f, %f\n",
           nodesTree[bboxNode->child[1]].bmin[0],
           nodesTree[bboxNode->child[1]].bmin[1],
           nodesTree[bboxNode->child[1]].bmin[2],
           nodesTree[bboxNode->child[1]].bmax[0],
           nodesTree[bboxNode->child[1]].bmax[1],
           nodesTree[bboxNode->child[1]].bmax[2]);

    printf("ray.org = %f, %f, %f\n", ray.org[0], ray.org[1], ray.org[2]);
    printf("ray.dir = %f, %f, %f\n", ray.dir[0], ray.dir[1], ray.dir[2]);
#endif

    if (TestRayAABB(tMinOut[1], tMaxOut[1], nodesTree[bboxNode->child[1]].bmin,
                    nodesTree[bboxNode->child[1]].bmax, ray)) {
      hitRight = 1;
    }
  }

  return (hitRight << 1) | hitLeft;
}

void BVHTree::ConstructTree(size_t indexRoot, size_t indexLeft,
                            size_t indexRight) {
  int axisList[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
  size_t n = indexRight - indexLeft;

  assert(n > 0);

  if (n == 1) {
    // Create leaf.

    m_nodesTree[indexRoot].isLeaf = true;
    m_nodesTree[indexRoot].child[0] = 0;
    m_nodesTree[indexRoot].child[1] = 0;
    m_nodesTree[indexRoot].nodeID = m_nodes[indexLeft].nodeID;

    return;
  }

  //
  // Spaital median.
  //
  double bsize[3];
  double bmin[3];

  bsize[0] = m_nodesTree[indexRoot].bmax[0] - m_nodesTree[indexRoot].bmin[0];
  bsize[1] = m_nodesTree[indexRoot].bmax[1] - m_nodesTree[indexRoot].bmin[1];
  bsize[2] = m_nodesTree[indexRoot].bmax[2] - m_nodesTree[indexRoot].bmin[2];

  bmin[0] = m_nodesTree[indexRoot].bmin[0];
  bmin[1] = m_nodesTree[indexRoot].bmin[1];
  bmin[2] = m_nodesTree[indexRoot].bmin[2];

  double maxLen = bsize[0];
  int axisTry = 0;

  if (maxLen < bsize[1]) {
    maxLen = bsize[1];
    axisTry = 1;
  }

  if (maxLen < bsize[2]) {
    maxLen = bsize[2];
    axisTry = 2;
  }

  size_t splitIndex = indexLeft;
  size_t nLeft = 0;
  size_t nRight = (indexRight - indexLeft);

  int splitAxis;
  double splitPos;

  for (int i = 0; i < 3; i++) {
    splitAxis = axisList[axisTry][i];
    splitPos = bmin[splitAxis] + bsize[splitAxis] * 0.5;

    splitIndex =
        PartitionNodes(m_nodes, indexLeft, indexRight, splitPos, splitAxis);

    // splitIndex = indexLeft = n / 2;
    if (splitIndex != indexLeft && splitIndex != indexRight) {
      // Found good split.
      nLeft = splitIndex - indexLeft;
      nRight = n - nLeft;
      break;
    }
  }

  if (splitIndex == indexLeft || splitIndex == indexRight) {
    // Cannot find good split. Force use object median.
    nLeft = n / 2;
    nRight = n - nLeft;
    splitIndex = indexLeft + (indexRight - indexLeft) / 2;
  }

  //
  // Subdivide
  //
  {

    m_nodesTree[indexRoot].axis = splitAxis;
    m_nodesTree[indexRoot].isLeaf = false;

    // left
    if (nLeft > 0) {
      BVHNode nodeLeft;

      double bminLeft[3], bmaxLeft[3];
      CalculateBoundingBox(bminLeft, bmaxLeft, m_nodes, indexLeft, splitIndex);

      // Slightly extend bbox to avoid numeric error.
      double eps = 1.0e-6;
      for (int i = 0; i < 3; i++) {
        double extent = bmaxLeft[i] - bminLeft[i];
        bminLeft[i] = bminLeft[i] - extent * eps;
        bmaxLeft[i] = bmaxLeft[i] + extent * eps;
      }

      nodeLeft.bmin[0] = bminLeft[0];
      nodeLeft.bmin[1] = bminLeft[1];
      nodeLeft.bmin[2] = bminLeft[2];
      nodeLeft.bmax[0] = bmaxLeft[0];
      nodeLeft.bmax[1] = bmaxLeft[1];
      nodeLeft.bmax[2] = bmaxLeft[2];

      size_t childLeftIndex = m_nodesTree.size();
      m_nodesTree[indexRoot].child[0] = childLeftIndex;

      m_nodesTree.push_back(nodeLeft);

      ConstructTree(childLeftIndex, indexLeft, splitIndex);

    } else {

      m_nodesTree[indexRoot].child[0] = 0;
    }

    // right
    if (nRight > 0) {
      BVHNode nodeRight;

      double bminRight[3], bmaxRight[3];
      CalculateBoundingBox(bminRight, bmaxRight, m_nodes, splitIndex,
                           indexRight);

      // Slightly extend bbox to avoid numeric error.
      double eps = 1.0e-6;
      for (int i = 0; i < 3; i++) {
        double extent = bmaxRight[i] - bminRight[i];
        bminRight[i] = bminRight[i] - extent * eps;
        bmaxRight[i] = bmaxRight[i] + extent * eps;
      }

      nodeRight.bmin[0] = bminRight[0];
      nodeRight.bmin[1] = bminRight[1];
      nodeRight.bmin[2] = bminRight[2];

      nodeRight.bmax[0] = bmaxRight[0];
      nodeRight.bmax[1] = bmaxRight[1];
      nodeRight.bmax[2] = bmaxRight[2];

      size_t childRightIndex = m_nodesTree.size();
      m_nodesTree[indexRoot].child[1] = childRightIndex;

      m_nodesTree.push_back(nodeRight);

      ConstructTree(childRightIndex, splitIndex, indexRight);

    } else {

      m_nodesTree[indexRoot].child[1] = 0;
    }
  }
}

//
//
//

void BVHTree::BuildTree() {
  double sceneBMin[3], sceneBMax[3];

  if (m_nodes.size() == 0) { // empty tree
  } else {

    CalcSceneBBox(sceneBMin, sceneBMax, m_nodes);

    // create root node.
    BVHNode rootNode;
    rootNode.bmin[0] = sceneBMin[0];
    rootNode.bmin[1] = sceneBMin[1];
    rootNode.bmin[2] = sceneBMin[2];

    rootNode.bmax[0] = sceneBMax[0];
    rootNode.bmax[1] = sceneBMax[1];
    rootNode.bmax[2] = sceneBMax[2];

    m_nodesTree.push_back(rootNode);

    ConstructTree(0, 0, m_nodes.size());
  }

  isBuiltTree = true;
}

bool BVHTree::Trace(std::vector<BVHNodeIntersection> &isects, double rayorg[3],
                    double raydir[3], double maxdist) const {
  uint64_t nodeStack[TREE_MAXDEPTH];
  uint64_t nodeStackIndex = 0;

  assert(isBuiltTree);

  if (m_nodes.size() == 0) { // empty tree
    return false;
  }

  // Precaluclate ray coeffs.
  Ray ray;
  ray.org[0] = rayorg[0];
  ray.org[1] = rayorg[1];
  ray.org[2] = rayorg[2];
  ray.dir[0] = raydir[0];
  ray.dir[1] = raydir[1];
  ray.dir[2] = raydir[2];
  SetupRay(ray);

  const BVHNode *root = &m_nodesTree[0];

  while (1) {

    if (root->isLeaf) {

      double tMin, tMax;
      // Use node's bbox(this is tighter than BBoxNode's bbox)
      bool hit = TestRayAABB(tMin, tMax, root->bmin, root->bmax, ray);
      if (hit && (tMin <= maxdist)) {

        BVHNodeIntersection isect;
        isect.nodeID = root->nodeID;
        isect.tMin = tMin;
        isect.tMax = tMax;
        isects.push_back(isect);
      }

      if (nodeStackIndex < 1)
        break;
      nodeStackIndex--;
      root = &m_nodesTree[nodeStack[nodeStackIndex]];

    } else {

      double tMin[2], tMax[2];
      int ret = TestRayNodeIntersection(tMin, tMax, m_nodesTree, root, ray);

      if (ret == 0) { // no hit

        if (nodeStackIndex < 1)
          break;
        nodeStackIndex--;
        root = &m_nodesTree[nodeStack[nodeStackIndex]];

      } else if (ret == 1) { // left

        root = &m_nodesTree[root->child[0]];

      } else if (ret == 2) { // right

        root = &m_nodesTree[root->child[1]];

      } else { // both

        int order = ray.sign[root->axis];

        nodeStack[nodeStackIndex] = root->child[1 - order];
        nodeStackIndex++;

        root = &m_nodesTree[root->child[order]];
      }
    }
  }

  //
  // Sort by tMin.
  //
  std::sort(isects.begin(), isects.end(), sorter_isect);

  return (isects.size() > 0);
}

bool BVHTree::Locate(std::vector<BVHNodeLocator> &locators,
                     const double position[3]) const {
  uint64_t nodeStack[TREE_MAXDEPTH];
  uint64_t nodeStackIndex = 0;

  assert(isBuiltTree);

  locators.clear();

  if (m_nodes.size() == 0) { // empty tree
    return false;
  }

  const BVHNode *root = &m_nodesTree[0];

  while (1) {

    if (root->isLeaf) {

      // Use node's bbox(this is tighter than BBoxNode's bbox)
      bool hit = TestPointAABB(root->bmin, root->bmax, position);
      if (hit) {

        BVHNodeLocator locator;
        locator.nodeID = root->nodeID;
        locators.push_back(locator);
      }

      if (nodeStackIndex < 1) {
        break;
      }
      nodeStackIndex--;
      root = &m_nodesTree[nodeStack[nodeStackIndex]];

    } else {

      int ret = TestPointBVHNode(m_nodesTree, root, position);

      if (ret == 0) { // no hit

        if (nodeStackIndex < 1)
          break;
        nodeStackIndex--;
        root = &m_nodesTree[nodeStack[nodeStackIndex]];

      } else if (ret == 1) { // left

        root = &m_nodesTree[root->child[0]];

      } else if (ret == 2) { // right

        root = &m_nodesTree[root->child[1]];

      } else { // both

        int order = 0;

        nodeStack[nodeStackIndex] = root->child[1 - order];
        nodeStackIndex++;

        root = &m_nodesTree[root->child[order]];
      }
    }
  }

  return (locators.size() > 0);
}

void BVHTree::BoundingBox(double bmin[3], double bmax[3]) const {
  if (m_nodesTree.size() < 1) {
    bmin[0] = bmin[1] = bmin[2] = 1.0e+30;
    bmax[0] = bmax[1] = bmax[2] = -1.0e+30;
  } else {
    const BVHNode *root = &m_nodesTree[0];
    bmin[0] = root->bmin[0];
    bmin[1] = root->bmin[1];
    bmin[2] = root->bmin[2];
    bmax[0] = root->bmax[0];
    bmax[1] = root->bmax[1];
    bmax[2] = root->bmax[2];
  }
}
