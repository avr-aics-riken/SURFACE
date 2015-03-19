/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_ACCEL_PARTICLE_H__
#define __LSGL_RENDER_ACCEL_PARTICLE_H__

#include <vector>

#include "render_common.h"
#include "render_prim_particle.h"
#include "render_ray.h"
#include "render_intersection.h"

namespace lsgl {
namespace render {

/// Represents particle BVH node.
class ParticleNode {
public:
  ParticleNode(){};
  ~ParticleNode(){};

  // 2-ary BVH
  real bmin[2][3];
  real bmax[2][3];

  // leaf
  //   data[0] = npoints
  //   data[1] = index
  //
  // branch
  //   data[0] = child[0]
  //   data[1] = child[1]
  unsigned int data[2];

  short flag; // 1 = leaf node, 0 = branch node
  short axis;
};

/// Particle build option.
struct ParticleBuildOptions {
  bool debugPrint;
  int minLeafPrimitives;
  int maxTreeDepth;

  // Set default value: Taabb = 0.2
  ParticleBuildOptions() : minLeafPrimitives(16), maxTreeDepth(256) {}
};

/// Particle build statistics.
struct ParticleBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  ParticleBuildStatistics()
      : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

/// Particle BVH acceleration class.
class ParticleAccel {
public:
  ParticleAccel(){};
  ~ParticleAccel(){};

  /// Build Particle for input particles.
  bool Build(const Particles *particles, const ParticleBuildOptions &options);

  /// Build Particle for input particles. Optimzied for ~2G points.
  bool Build32(const Particles *particles, const ParticleBuildOptions &options);

  /// Get statistics of built Particle tree. Valid after Build()
  ParticleBuildStatistics GetStatistics() const { return stats_; }

  /// Dump built Particle to the file.
  bool Dump(const char *filename);

  /// Traverse into Particle along ray and find closest hit point if found
  bool Traverse(Intersection &isect, Ray &ray) const;

  /// Get the bounding box of BVH. Valid after Build().
  void BoundingBox(double bmin[3], double bmax[3]);

private:
  /// Builds particle tree recursively.
  size_t BuildTree(const Particles *particles, real3 &bmin, real3 &bmax,
                   unsigned int leftIdx, unsigned int rightIdx, int depth);

#if 0
  /// Builds blob particle tree recursively.
  size_t BuildTreeBlob(const Particles *particles,
                       std::vector<unsigned int> &indices, unsigned int prevNum,
                       int depth);
#endif

  ParticleBuildOptions options_;
  std::vector<ParticleNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G paritcles.
  ParticleBuildStatistics stats_;
  const Particles *particles_;

  double bmin_[3], bmax_[3]; // Bounding box of whole tree
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_RENDER_ACCEL_PARTICLE_H__
