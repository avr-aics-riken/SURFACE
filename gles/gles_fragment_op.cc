/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>

#include "GLES2/gl2.h"
#include "gles_context.h"

using namespace lsgl;

void Context::glColorMask(GLboolean red, GLboolean green, GLboolean blue,
                          GLboolean alpha) {
  TRACE_EVENT("(GLboolean red = %d, GLboolean green = %d, GLboolean blue = %d, "
              "GLboolean alpha = %d)",
              red, green, blue, alpha);
  state_.colorMaskRed = (bool)red;
  state_.colorMaskGreen = (bool)green;
  state_.colorMaskBlue = (bool)blue;
  state_.colorMaskAlpha = (bool)alpha;
}

void Context::glDepthFunc(GLenum func) {
  TRACE_EVENT("(GLenum func = %d)", func);
  state_.depthFunc = func;
}

void Context::glDepthMask(GLboolean flag) {
  TRACE_EVENT("(GLboolean flag = %d)", flag);
  state_.depthMask = flag;
}

void Context::glDepthRangef(GLclampf zNear, GLclampf zFar) {
  TRACE_EVENT("(GLclampf zNear = %f, GLclampf zFar = %f)", zNear, zFar);
  state_.nearDepth = zNear;
  state_.farDepth = zFar;
}

void Context::glBlendColor(GLclampf red, GLclampf green, GLclampf blue,
                           GLclampf alpha) {
  TRACE_EVENT("(GLclampf red = %f, GLclampf green = %f, GLclampf blue = %f, "
              "GLclampf alpha = %f)",
              red, green, blue, alpha);
  state_.blendColor.red = red;
  state_.blendColor.green = green;
  state_.blendColor.blue = blue;
  state_.blendColor.alpha = alpha;
}

void Context::glBlendEquation(GLenum mode) {
  state_.blendEquationRGB = mode;
  state_.blendEquationAlpha = mode;
}

void Context::glBlendEquationSeparate(GLenum modeRGB, GLenum modeAlpha) {
  state_.blendEquationRGB = modeRGB;
  state_.blendEquationAlpha = modeAlpha;
}

void Context::glBlendFunc(GLenum sfactor, GLenum dfactor) {
  state_.sourceBlendRGB = sfactor;
  state_.sourceBlendAlpha = sfactor;
  state_.destBlendRGB = dfactor;
  state_.destBlendAlpha = dfactor;
}

void Context::glBlendFuncSeparate(GLenum srcRGB, GLenum dstRGB, GLenum srcAlpha,
                                  GLenum dstAlpha) {
  state_.sourceBlendRGB = srcRGB;
  state_.sourceBlendAlpha = srcAlpha;
  state_.destBlendRGB = dstRGB;
  state_.destBlendAlpha = dstAlpha;
}

void Context::glScissor(GLint x, GLint y, GLsizei width, GLsizei height) {
  TRACE_EVENT(
      "(GLint x = %d, GLint y = %d, GLsizei width = %d, GLsizei height = %d)",
      x, y, width, height);
  state_.scissorX = x;
  state_.scissorY = y;
  state_.scissorWidth = width;
  state_.scissorHeight = height;
}

void Context::glStencilFunc(GLenum func, GLint ref, GLuint mask) {
  glStencilFuncSeparate(GL_FRONT_AND_BACK, func, ref, mask);
}

void Context::glStencilFuncSeparate(GLenum face, GLenum func, GLint ref,
                                    GLuint mask) {
  // stencil function front face setting
  if ((face == GL_FRONT) || (face == GL_FRONT_AND_BACK)) {
    state_.stencilFunc = func;
    state_.stencilRef = ref;
    state_.stencilMask = mask;
  }

  // stencil function back face setting
  if ((face == GL_BACK) || (face == GL_FRONT_AND_BACK)) {
    state_.stencilBackFunc = func;
    state_.stencilBackRef = ref;
    state_.stencilBackMask = mask;
  }
}

void Context::glStencilMask(GLuint mask) {
  glStencilMaskSeparate(GL_FRONT_AND_BACK, mask);
}

void Context::glStencilMaskSeparate(GLenum face, GLuint mask) {
  // stencil mask front face setting
  if ((face == GL_FRONT) || (face == GL_FRONT_AND_BACK)) {
    state_.stencilWritemask = mask;
  }

  // stencil mask back face setting
  if ((face == GL_BACK) || (face == GL_FRONT_AND_BACK)) {
    state_.stencilBackWritemask = mask;
  }
}

void Context::glStencilOp(GLenum fail, GLenum zfail, GLenum zpass) {
  glStencilOpSeparate(GL_FRONT_AND_BACK, fail, zfail, zpass);
}

void Context::glStencilOpSeparate(GLenum face, GLenum fail, GLenum zfail,
                                  GLenum zpass) {
  // stencil operation front face setting
  if ((face == GL_FRONT) || (face == GL_FRONT_AND_BACK)) {
    state_.stencilFail = fail;
    state_.stencilPassDepthFail = zfail;
    state_.stencilPassDepthPass = zpass;
  }

  // stencil operation back face setting
  if ((face == GL_BACK) || (face == GL_FRONT_AND_BACK)) {
    state_.stencilBackFail = fail;
    state_.stencilBackPassDepthFail = zfail;
    state_.stencilBackPassDepthPass = zpass;
  }
}

void Context::glSampleCoverage(GLclampf value, GLboolean invert) {
  TRACE_EVENT("(GLclampf value = %f, GLboolean invert %d)", value, invert);

  // ensure value is valid
  if ((value < 1.0f) || (value > kMaxSamplesPerPixel))
    return SetGLError(GL_INVALID_VALUE);

  state_.sampleCoverageValue = value;
  state_.sampleCoverageInvert = invert;
}
