/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_GLES_RAYTRACING_ENGINE_H__
#define __LSGL_GLES_RAYTRACING_ENGINE_H__

#include <set>
#include <map>
#include <vector>
#include <string>

#include "gles_common.h"
#include "gles_accel_builder.h"
#include "gles_render_graph.h"

#include "../render/render_bvh_tree.h"

namespace lsgl {

class Context;
class RenderGraph;

/// Raytracing system
class RaytraceEngine {
public:
  RaytraceEngine();
  ~RaytraceEngine();

  /// Set camera where ray is traced from.
  void SetCamera(const GLfloat *eye, const GLfloat *target, const GLfloat *up,
                 GLfloat fov);

  /// Set panoramic environment camera.
  void SetStereoEnvCamera(const GLfloat *eye, const GLfloat *target,
                          const GLfloat *up, GLfloat zeroParallax,
                          GLfloat eyeSeparation);

  // Evented method
  //
  // OnPrepare ---> OnStart --+--> OnData --+--> OnEnd
  //                          |             |
  //                          +---- <- -----+
  //
  void OnPrepare(const Context *ctx);

  // @todo { remove }
  // bool OnStart(Framebuffer *fb,
  //             const std::vector<MeshBuilder::RenderMesh> &meshList);
  bool OnStart(Framebuffer *fb, RenderGraph *renderGraph);
  bool OnData();
  void OnEnd();

  // @todo { remove }
  // void SetFragmentState(FragmentState &fragmentState) {
  //  fragmentState_ = fragmentState;
  //}

  // void SetShadingState(ShadingState &shadingState) {
  //  shadingState_ = shadingState;
  //}

  bool Trace(Intersection &isect, Ray &r);

  static inline RaytraceEngine *GetRaytraceEngine() { return sRaytraceEngine; }

  /// Shade one fragment object for the intersection. Shaded color is stored in
  /// shadecol.
  void ShadeFragment(float shadecol[4], bool &discarded,
                     const Intersection &isect, const Ray &ray, int thread_id);

  RenderGraph *GetRenderGraph() const { return renderGraph_; }

  void SetPixelStep(int step) { pixelStep_ = step; }

  int GetPixelStep() const { return pixelStep_; }

  void SetProgressCallback(LSGLProgressCallback func, void *userdata) {
    progressCallbackFunc_ = func;
    callbackUserData_ = userdata;
  }

  void SetScreenParallelRendering(bool enabled) {
	screenParallelRendering_ = enabled;
  }

private:
  const Context *ctx_;
  Framebuffer *framebuffer_;

  Camera *camera_; // @todo { Support stereo camera? }

  RenderGraph *renderGraph_;

  timerutil renderTimer_;
  double accumTime_;
  double numRays_;

  static RaytraceEngine *sRaytraceEngine;

  int pixelStep_;
  LSGLProgressCallback progressCallbackFunc_;
  void *callbackUserData_;

  bool screenParallelRendering_;
};

} // namespace

#endif // __LSGL_GLES_RAYTRACING_ENGINE_H__
