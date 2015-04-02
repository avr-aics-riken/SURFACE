/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cstring>

#include "gles_render_graph.h"
#include "../render/render_matrix.h"

using namespace lsgl;

static void TransformBoundingBox(double xbmin[3], // out
                                 double xbmax[3], // out
                                 double bmin[3], double bmax[3],
                                 double m[4][4]) {

  // create bounding vertex from (bmin, bmax)
  double b[8][3];

  b[0][0] = bmin[0];
  b[0][1] = bmin[1];
  b[0][2] = bmin[2];
  b[1][0] = bmax[0];
  b[1][1] = bmin[1];
  b[1][2] = bmin[2];
  b[2][0] = bmin[0];
  b[2][1] = bmax[1];
  b[2][2] = bmin[2];
  b[3][0] = bmax[0];
  b[3][1] = bmax[1];
  b[3][2] = bmin[2];

  b[4][0] = bmin[0];
  b[4][1] = bmin[1];
  b[4][2] = bmax[2];
  b[5][0] = bmax[0];
  b[5][1] = bmin[1];
  b[5][2] = bmax[2];
  b[6][0] = bmin[0];
  b[6][1] = bmax[1];
  b[6][2] = bmax[2];
  b[7][0] = bmax[0];
  b[7][1] = bmax[1];
  b[7][2] = bmax[2];

  double xb[8][3];
  for (int i = 0; i < 8; i++) {
    Matrixd::MultV(xb[i], m, b[i]);
  }

  xbmin[0] = xb[0][0];
  xbmin[1] = xb[0][1];
  xbmin[2] = xb[0][2];
  xbmax[0] = xb[0][0];
  xbmax[1] = xb[0][1];
  xbmax[2] = xb[0][2];

  for (int i = 1; i < 8; i++) {

    xbmin[0] = (std::min)(xb[i][0], xbmin[0]);
    xbmin[1] = (std::min)(xb[i][1], xbmin[1]);
    xbmin[2] = (std::min)(xb[i][2], xbmin[2]);

    xbmax[0] = (std::max)(xb[i][0], xbmax[0]);
    xbmax[1] = (std::max)(xb[i][1], xbmax[1]);
    xbmax[2] = (std::max)(xb[i][2], xbmax[2]);
  }
}

void RenderElement::Update() {
  // Compute world space bounding box.
  TransformBoundingBox(xbmin_, xbmax_, bmin_, bmax_, T_);

  //
  // Precompute tranformation matrices
  //

  // inverse of T
  memcpy(invT_, T_, sizeof(double) * 4 * 4);
  Matrixd::Inverse(invT_);

  memcpy(invT33_, T_, sizeof(double) * 4 * 4);

  // inverset of upper 3x3 transform.
  // clear translation
  invT33_[3][0] = 0.0;
  invT33_[3][1] = 0.0;
  invT33_[3][2] = 0.0;

  Matrixd::Inverse(invT33_);

  // inverse transpose of upper 3x3 transform.
  memcpy(invTransposeT33_, invT33_, sizeof(double) * 4 * 4);
  Matrixd::Transpose(invTransposeT33_);
}

bool RenderGraph::Build() {
  isBuiltGraph_ = true;

  delete tree_;
  tree_ = new render::BVHTree();

  // Add node info to toplevel bvh tree.
  for (size_t i = 0; i < renderElements_.size(); i++) {
    render::BVHData node;

    // Use world space bbox.
    node.bmin[0] = renderElements_[i].xbmin_[0];
    node.bmin[1] = renderElements_[i].xbmin_[1];
    node.bmin[2] = renderElements_[i].xbmin_[2];
    node.bmax[0] = renderElements_[i].xbmax_[0];
    node.bmax[1] = renderElements_[i].xbmax_[1];
    node.bmax[2] = renderElements_[i].xbmax_[2];
    node.nodeID = renderElements_[i].nodeID_;

    printf("[LSGL] DBG: node[%d] bmin (%f, %f, %f).\n", (int)i, node.bmin[0],
           node.bmin[1], node.bmin[2]);
    printf("[LSGL] DBG: node[%d] bmax (%f, %f, %f).\n", (int)i, node.bmax[0],
           node.bmax[1], node.bmax[2]);

    tree_->AddNode(node);
  }

  printf("[LSGL] DBG: added %d nodes.\n", (int)renderElements_.size());

  tree_->BuildTree();

  double bmin[3], bmax[3];

  tree_->BoundingBox(bmin, bmax);

  printf("[LSGL] DBG: bmin (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL] DBG: bmax (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);

  return false;
}

bool RenderGraph::Trace(Intersection &isect, Ray &ray) {
  std::vector<render::BVHNodeIntersection> isects;

  if (!tree_ || !isBuiltGraph_) {
    return false;
  }

  double maxdist = (std::numeric_limits<double>::max)();
  double rayorg[3];
  rayorg[0] = ray.origin()[0];
  rayorg[1] = ray.origin()[1];
  rayorg[2] = ray.origin()[2];

  double raydir[3];
  raydir[0] = ray.direction()[0];
  raydir[1] = ray.direction()[1];
  raydir[2] = ray.direction()[2];

  // Get list of potential hits by tracing toplevel BVH.
  bool mayHit = tree_->Trace(isects, rayorg, raydir, maxdist);

  if (mayHit) {

    double tMax = (std::numeric_limits<double>::max)();
    double tNearest = tMax;
    unsigned int hitNodeID;
    bool hasHit = false;
    const char *renderElement = NULL;

    // Find detailed intersection.
    // Note that intersection list are already sorted in its distance.
    for (size_t i = 0; i < isects.size(); i++) {

      // Early cull test.
      if (tNearest < isects[i].tMin) {
        continue;
      }

      bool hit = false;

      RenderElement *renderElement = &renderElements_[isects[i].nodeID];

      Intersection bottomIsect;
      bottomIsect.clear();
      bottomIsect.t = (std::numeric_limits<float>::max)();

      // Transform ray into local space(world -> local)
      double rayOrg[3];
      double rayDir[3];
      double localRayOrg[3];
      double localRayDir[3];

      rayOrg[0] = ray.origin()[0];
      rayOrg[1] = ray.origin()[1];
      rayOrg[2] = ray.origin()[2];

      rayDir[0] = ray.direction()[0];
      rayDir[1] = ray.direction()[1];
      rayDir[2] = ray.direction()[2];

      Matrixd::MultV(localRayOrg, renderElement->invT_, rayOrg);
      Matrixd::MultV(localRayDir, renderElement->invT33_, rayDir);

      Ray localRay(localRayOrg, localRayDir);
      localRay.double_sided = ray.double_sided;
      localRay.depth = ray.depth;

      // If previous ray hits this render element node, enable prev_prim_id
      // for self-intersection test
      if (ray.prev_node &&
          (ray.prev_node ==
           reinterpret_cast<const unsigned char *>(renderElement))) {
        localRay.prev_node = ray.prev_node;
        localRay.prev_prim_id = ray.prev_prim_id;
      }

      // @todo { optimzie }
      switch (renderElement->GetPrimType()) {
      case AccelBuilder::PRIMITIVE_TRIANGLES: {
        const AccelBuilder::MeshAccelerator *accel =
            reinterpret_cast<const AccelBuilder::MeshAccelerator *>(
                renderElement->GetAccel());
        if (accel) {
          hit = accel->Traverse(bottomIsect, localRay);
        }
      } break;
      case AccelBuilder::PRIMITIVE_POINTS: {
        const AccelBuilder::ParticleAccelerator *accel =
            reinterpret_cast<const AccelBuilder::ParticleAccelerator *>(
                renderElement->GetAccel());
        hit = accel->Traverse(bottomIsect, localRay);
      } break;
      case AccelBuilder::PRIMITIVE_LINES: {
        const AccelBuilder::LineAccelerator *accel =
            reinterpret_cast<const AccelBuilder::LineAccelerator *>(
                renderElement->GetAccel());
        hit = accel->Traverse(bottomIsect, localRay);
      } break;
      case AccelBuilder::PRIMITIVE_TETRAHEDRONS: {
        const AccelBuilder::TetraAccelerator *accel =
            reinterpret_cast<const AccelBuilder::TetraAccelerator *>(
                renderElement->GetAccel());
        hit = accel->Traverse(bottomIsect, localRay);
      } break;
      default:
        fprintf(stderr, "[LSGL] Unknown primitive type.\n");
        assert(0);
        break;
      }

      if (hit) {

        // First calulcate distance in world coordiante.
        real3 localP;
        localP[0] = localRayOrg[0] + bottomIsect.t * localRayDir[0];
        localP[1] = localRayOrg[1] + bottomIsect.t * localRayDir[1];
        localP[2] = localRayOrg[2] + bottomIsect.t * localRayDir[2];

        real3 worldP;
        Matrixd::MultV(worldP, renderElement->T_, localP);

        real3 po = worldP - ray.origin();

        float tWorld = (worldP - ray.origin()).length();

        // printf("lp = %f, %f, %f, lpLen = %f\n", localP[0], localP[1],
        // localP[2], localP.length());
        // printf("wp = %f, %f, %f, wpLen = %f\n", worldP[0], worldP[1],
        // worldP[2], worldP.length());
        // printf("po = %f, %f, %f, poLen = %f\n", po[0], po[1], po[2],
        // po.length());
        // printf("bt.t = %f, t = %f\n", bottomIsect.t, tWorld);

        if (tWorld < tNearest) {
          tNearest = tWorld;
          hasHit = true;
          isect = bottomIsect;
          isect.renderElement =
              reinterpret_cast<unsigned char *>(renderElement);

          // Transform intersection info(local -> world)
          isect.t = tWorld;
          Matrixd::MultV(isect.position, renderElement->T_, isect.position);
          Matrixd::MultV(isect.geometric, renderElement->invTransposeT33_,
                         isect.geometric);
          Matrixd::MultV(isect.normal, renderElement->invTransposeT33_,
                         isect.normal);
        }
      }
    }

    if (hasHit) {
      // isect is already filled.
      return true;
    }
  }

  return false;
}
