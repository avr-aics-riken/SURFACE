/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>

#include "GLES2/gl2.h"
#include "gles_context.h"

using namespace lsgl;

void Context::glFlush() {
  TRACE_EVENT("()");
  // simply do nothing
}

void Context::glFinish() {
  TRACE_EVENT("()");
  // no need to do anything if our frame has not changed
  if (dirtyFrame_ == false)
    return;

  // no need to do anything if there's no shader
  Program *program = resourceManager_.GetProgram(state_.currentProgram);
  if (program == NULL) {
    return;
  }

  // build render graph
  assert(renderGraph_);
  renderGraph_->Build();

  // build sparse volume texture.
  for (size_t i = 0; i < sparseTextureList_.size(); i++) {
    Texture *tex = sparseTextureList_[i];
    tex->BuildSparseTexture();
  }

  // render setup!
  bool ret = engine_.OnStart(GetCurrentFrameBuffer(), renderGraph_);

  if (ret) {
    engine_.SetPixelStep(state_.pixelStep);

    // trace!
    engine_.OnData();
    engine_.OnEnd();
  }

  // Free internally generated buffers
  for (size_t i = 0; i < bufferFreeList_.size(); i++) {
    GLuint idx = bufferFreeList_[i];
    glDeleteBuffers(1, &idx);
  }
  bufferFreeList_.clear();

  // let the mesh builder know the frame ended, and clear our render list
  accelBuilder_.EndFrame();

  state_.currentDrawStackIndex = 0;

  // clear render graph.
  renderGraph_->Clear();

  // clear the dirty frame flag
  dirtyFrame_ = false;

  if (!ret) {
    return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
  }
}

void Context::glHint(GLenum target, GLenum mode) {
  // legacy
  TRACE_EVENT("(GLenum target = %d, GLenum mode = %d)", target, mode);
}
