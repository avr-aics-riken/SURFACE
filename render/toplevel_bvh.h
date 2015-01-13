/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __RENDER_TOPLEVEL_BVH_H__
#define __RENDER_TOPLEVEL_BVH_H__

#include "render_common.h"
#include "matrix.h"

//
// Top level BVH is a BVH of bottom level BVH nodes.
// This class is vital to represent scene graph of input scene primitives.
//

#include <ctime>
#include <vector>

namespace lsgl {
namespace render {

struct ToplevelBVHNode {

  // local bbox
  double bmin[3];
  double bmax[3];

  // Child index offset
  int child[2]; ///< 0 or -1 = leaf.

  // Inner
  int axis; ///< Split axis

  // Leaf
  bool isLeaf;
  unsigned int nodeID;
};

struct ToplevelBVHData {
  // local bbox
  double bmin[3];
  double bmax[3];

  unsigned int nodeID; ///< Node ID
};

struct ToplevelBVHNodeIntersection {
  double tMin, tMax;   ///< Enter and exit distance.
  unsigned int nodeID; ///< BVH node id
  unsigned int rayID;
};

class ToplevelBVHTree {
public:
  ToplevelBVHTree() : isBuiltTree(false) {}
  ~ToplevelBVHTree() {}

  // Add node to the list.
  void AddNode(ToplevelBVHData &node) {
    m_nodes.push_back(node);
    m_orgnodes.push_back(node);
  }

  // Build toplevel BVH tree. Valid after adding nodes using AddNode().
  void BuildTree();

  // Trace ray into toplevel BVH tree and list up possible intersection list.
  bool Trace(std::vector<ToplevelBVHNodeIntersection> &isects, double rayorg[3],
             double raydir[3], double maxdist);

  std::vector<ToplevelBVHNode> &GetNodesTree() { return m_nodesTree; }

  void SetNodesTree(std::vector<ToplevelBVHNode> &nodesTree) {
    m_nodesTree = nodesTree;
  }

  // Get bounding box of this toplevel BVH tree.
  void BoundingBox(double bmin[3], double bmax[3]);

  // for testing. @fixme { remove. }
  std::vector<ToplevelBVHData> m_orgnodes;

private:
  void ConstructTree(size_t indexRoot, size_t indexLeft, size_t indexRight);

  // Tree reprenstation of ToplevelBVH
  std::vector<ToplevelBVHNode> m_nodesTree;

  // Input ToplevelBVH list.
  std::vector<ToplevelBVHData> m_nodes;

  bool isBuiltTree;
};

} // namespace render
} // namespace lsgl

#endif // __RENDER_TOPLEVEL_BVH_H__
