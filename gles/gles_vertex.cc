/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cassert>

#include "gles_common.h"
#include "gles_context.h"

using namespace lsgl;

void Context::glVertexAttrib1f(GLuint indx, GLfloat x) {
  TRACE_EVENT("(GLuint indx = %d, GLfloat x = %f)", indx, x);

  assert(indx < kMaxVertexAttribs);

  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].size = 1;
  state_.vertexAttributes[k][indx].type = GL_FLOAT;
  state_.vertexAttributes[k][indx].ptr = NULL;
  state_.vertexAttributes[k][indx].uniform = true;
  state_.vertexAttributes[k][indx].u.fval1 = x;
}

void Context::glVertexAttrib2f(GLuint indx, GLfloat x, GLfloat y) {
  TRACE_EVENT("(GLuint indx = %d, GLfloat x = %f, GLfloat y = %f)", indx, x, y);

  assert(indx < kMaxVertexAttribs);

  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].size = 1;
  state_.vertexAttributes[k][indx].type = GL_FLOAT_VEC2;
  state_.vertexAttributes[k][indx].ptr = NULL;
  state_.vertexAttributes[k][indx].uniform = true;
  state_.vertexAttributes[k][indx].u.fval2[0] = x;
  state_.vertexAttributes[k][indx].u.fval2[1] = y;
}

void Context::glVertexAttrib3f(GLuint indx, GLfloat x, GLfloat y, GLfloat z) {
  TRACE_EVENT(
      "(GLuint indx = %d, GLfloat x = %f, GLfloat y = %f, GLfloat z = %f)",
      indx, x, y, z);

  assert(indx < kMaxVertexAttribs);

  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].size = 1;
  state_.vertexAttributes[k][indx].type = GL_FLOAT_VEC3;
  state_.vertexAttributes[k][indx].ptr = NULL;
  state_.vertexAttributes[k][indx].uniform = true;
  state_.vertexAttributes[k][indx].u.fval3[0] = x;
  state_.vertexAttributes[k][indx].u.fval3[1] = y;
  state_.vertexAttributes[k][indx].u.fval3[2] = z;
}

void Context::glVertexAttrib4f(GLuint indx, GLfloat x, GLfloat y, GLfloat z,
                               GLfloat w) {
  TRACE_EVENT("(GLuint indx = %d, GLfloat x = %f, GLfloat y = %f, GLfloat z = "
              "%f, GLfloat w = %f)",
              indx, x, y, z, w);

  assert(indx < kMaxVertexAttribs);

  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].size = 1;
  state_.vertexAttributes[k][indx].type = GL_FLOAT_VEC4;
  state_.vertexAttributes[k][indx].ptr = NULL;
  state_.vertexAttributes[k][indx].uniform = true;
  state_.vertexAttributes[k][indx].u.fval4[0] = x;
  state_.vertexAttributes[k][indx].u.fval4[1] = y;
  state_.vertexAttributes[k][indx].u.fval4[2] = z;
  state_.vertexAttributes[k][indx].u.fval4[3] = w;
}

void Context::glVertexAttrib1fv(GLuint indx, const GLfloat *values) {
  TRACE_EVENT("(GLuint indx = %d, const GLfloat* values = %p)", indx, values);
  assert(0 && "VertexAttrib1fv is not yet supported");
  assert(indx < kMaxVertexAttribs);
}

void Context::glVertexAttrib2fv(GLuint indx, const GLfloat *values) {
  TRACE_EVENT("(GLuint indx = %d, const GLfloat* values = %p)", indx, values);
  assert(0 && "VertexAttrib2fv is not yet supported");
  assert(indx < kMaxVertexAttribs);
}

void Context::glVertexAttrib3fv(GLuint indx, const GLfloat *values) {
  TRACE_EVENT("(GLuint indx = %d, const GLfloat* values = %p)", indx, values);
  assert(0 && "VertexAttrib3fv is not yet supported");
  assert(indx < kMaxVertexAttribs);
}

void Context::glVertexAttrib4fv(GLuint indx, const GLfloat *values) {
  TRACE_EVENT("(GLuint indx = %d, const GLfloat* values = %p)", indx, values);
  assert(0 && "VertexAttrib4fv is not yet supported");
  assert(indx < kMaxVertexAttribs);
}

void Context::glVertexAttribPointer(GLuint indx, GLint size, GLenum type,
                                    GLboolean normalized, GLsizei stride,
                                    const GLvoid *ptr) {
  TRACE_EVENT("(GLuint indx = %d, GLint size = %d, GLenum type = %d, GLboolean "
              "normalized = %d, GLsizei stride = %d, const GLvoid* ptr = %p)",
              indx, size, type, normalized, stride, ptr);
  assert(indx < kMaxVertexAttribs);

  int handle = -1;
  if (ptr == NULL) {
    // vertex buffer
    handle = state_.arrayBuffer;
  }

  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].size = size;
  state_.vertexAttributes[k][indx].type = type;
  state_.vertexAttributes[k][indx].normalized = normalized;
  state_.vertexAttributes[k][indx].stride = stride;
  state_.vertexAttributes[k][indx].ptr = ptr;
  state_.vertexAttributes[k][indx].uniform = false;
  state_.vertexAttributes[k][indx].handle = handle;
}

void Context::glEnableVertexAttribArray(GLuint indx) {
  TRACE_EVENT("(GLuint indx = %d)", indx);
  assert(indx < kMaxVertexAttribs);
  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].enabled = true;
}

void Context::glDisableVertexAttribArray(GLuint indx) {
  TRACE_EVENT("(GLuint indx = %d)", indx);
  assert(indx < kMaxVertexAttribs);
  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);
  state_.vertexAttributes[k][indx].enabled = false;
}
