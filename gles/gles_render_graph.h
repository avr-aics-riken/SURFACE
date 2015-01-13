/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_GLES_RENDER_GRAPH_H__
#define __LSGL_GLES_RENDER_GRAPH_H__

#include <vector>
#include <cstring>

#include "gles_accel_builder.h"
#include "../render/toplevel_bvh.h"
#include "../render/intersection.h"

namespace lsgl {

//
// Render graph which manage internal render resources for raytracing engine.
//

class RenderElement {
public:
  RenderElement(AccelBuilder::PrimitiveType type, unsigned char *accel,
                double bmin[3], double bmax[3], double T[4][4],
                int drawStackIndex, const Program *program,
                const FragmentState &fragmentState,
                const ShadingState &shadingState)
      : primType_(type), accel_(accel), drawStackIndex_(drawStackIndex),
        program_(program), fragmentState_(fragmentState),
        shadingState_(shadingState), nodeID_(-1) {

    bmin_[0] = bmin[0];
    bmin_[1] = bmin[1];
    bmin_[2] = bmin[2];

    bmax_[0] = bmax[0];
    bmax_[1] = bmax[1];
    bmax_[2] = bmax[2];

    memcpy(T_, T, sizeof(double) * 4 * 4);
    Update();

    // nodeID_ will be assigned in the later phase
  }
  ~RenderElement() {}

  // Update internal transformation matrix.
  void Update();

  inline AccelBuilder::PrimitiveType GetPrimType() const { return primType_; }
  inline const unsigned char *GetAccel() const { return accel_; }
  inline const Program *GetProgram() const { return program_; }
  inline void GetBoundingBox(double bmin[3], double bmax[3]) const {
    bmin[0] = xbmin_[0];
    bmin[1] = xbmin_[1];
    bmin[2] = xbmin_[2];
    bmax[0] = xbmax_[0];
    bmax[1] = xbmax_[1];
    bmax[2] = xbmax_[2];
  }

  inline void GetLocalBoundingBox(double bmin[3], double bmax[3]) const {
    bmin[0] = bmin_[0];
    bmin[1] = bmin_[1];
    bmin[2] = bmin_[2];
    bmax[0] = bmax_[0];
    bmax[1] = bmax_[1];
    bmax[2] = bmax_[2];
  }

  inline ShadingState &GetShadingState() { return shadingState_; }

  inline FragmentState &GetFragmentState() { return fragmentState_; }

  inline int GetDrawStackIndex() const { return drawStackIndex_; }

  AccelBuilder::PrimitiveType primType_;
  const unsigned char *accel_; // Opaque pointer for accelerator.
  int nodeID_;
  int drawStackIndex_; // For vertex attributes
  const Program *program_;
  ShadingState shadingState_;
  FragmentState fragmentState_;

  // bounding box(local space)
  double bmin_[3];
  double bmax_[3];

  // bounding box after xform(world space)
  double xbmin_[3];
  double xbmax_[3];

  // Xform
  double T_[4][4];      /// Transformation matrix(local -> world)
  double invT_[4][4];   /// Inverse(T) (world -> local)
  double invT33_[4][4]; /// Inverse of T's upper 3x3 element (world -> local,
                        /// used for transforming direction vector)
  double invTransposeT33_[4][4]; /// Inverse transpose of T's upper 3x3 element
                                 /// (world -> local, used
                                 /// for transformin normal vector)
};

class RenderGraph {
public:
  RenderGraph() : isBuiltGraph_(false), tree_(NULL) {}
  ~RenderGraph() {
    if (isBuiltGraph_) {
      delete tree_;
    }
  }

  bool AddRenderElement(RenderElement &renderElement) {
    // Assign unique node ID. nodeID == array index.
    renderElement.nodeID_ = (int)renderElements_.size();
    renderElements_.push_back(renderElement);
    return true;
  }

  // Build spatial accelerator for ray tracing.
  // BuildAccel() should be called after adding render elements with
  // AddRenderElement().
  bool Build();

  // Trace ray into this render graph.
  // This functions is valid after BuildAccel() was called.
  bool Trace(Intersection &isect, Ray &ray);

  bool IsBuilt() { return isBuiltGraph_; }

  void Clear() {
    renderElements_.clear();
    isBuiltGraph_ = false;
  }

private:
  std::vector<RenderElement> renderElements_;
  render::ToplevelBVHTree *tree_;

  bool isBuiltGraph_;
};

} // namespace lsgl

#endif // __LSGL_GLES_RENDER_GRAPH_H__
