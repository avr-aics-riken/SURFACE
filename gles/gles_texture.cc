/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "gles_context.h"

using namespace lsgl;

#include "../render/texture.h"

Texture *Context::HandleToTexture(GLenum target) {
  // get current selected texture index
  const GLuint idx = state_.activeSampler - GL_TEXTURE0;

  // get currently bound handle
  GLuint handle;
  switch (target) {
  case GL_TEXTURE_1D:
    handle = state_.texture1D[idx];
    break;

  case GL_TEXTURE_2D:
    handle = state_.texture2D[idx];
    break;

  case GL_TEXTURE_3D:
    handle = state_.texture3D[idx];
    break;

  case GL_TEXTURE_CUBE_MAP:
    handle = state_.textureCubeMap[idx];
    break;

  default:
    SetGLError(GL_INVALID_ENUM);
    return NULL;
  }

  // ensure bound handle is valid
  if (handle == 0) {
    SetGLError(GL_INVALID_OPERATION);
    return NULL;
  }

  // get texture pointer from handle
  Texture *tex = resourceManager_.GetTexture(handle);
  if (tex == NULL) {
    SetGLError(GL_INVALID_OPERATION);
    return NULL;
  }

  return tex;
}

void Context::glActiveTexture(GLenum texture) {
  TRACE_EVENT("(GLenum texture = %d)", texture);

  if ((texture < GL_TEXTURE0) ||
      (texture >= (GL_TEXTURE0 + kMaxTextureImageUnits))) {
    return SetGLError(GL_INVALID_ENUM);
  }

  state_.activeSampler = texture;
}

void Context::glGenerateMipmap(GLenum target) {
  TRACE_EVENT("(GLenum target = %d)", target);
  // mipmaps are not currently supported
  fprintf(stderr, "[LSGLES] glGenerateMipmap() is not currently supported.\n");
}

void Context::glGenTextures(GLsizei n, GLuint *textures) {
  TRACE_EVENT("(GLsizei n = %d, GLuint* textures = %p)", n, textures);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(textures && "textures are null");

  for (int i = 0; i < n; i++) {
    textures[i] = resourceManager_.CreateTexture();
  }
}

void Context::glBindTexture(GLenum target, GLenum texture) {
  TRACE_EVENT("(GLenum target = %d, GLenum texture = %d)", target, texture);

  // get current selected texture index
  const GLuint idx = state_.activeSampler - GL_TEXTURE0;

  // store new texture handle
  switch (target) {
  case GL_TEXTURE_1D:
    state_.texture1D[idx] = texture;
    break;

  case GL_TEXTURE_2D:
    state_.texture2D[idx] = texture;
    break;

  case GL_TEXTURE_3D:
    state_.texture3D[idx] = texture;
    break;

  case GL_TEXTURE_CUBE_MAP:
    state_.textureCubeMap[idx] = texture;
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glTexImage2D(GLenum target, GLint level, GLint internalformat,
                           GLsizei width, GLsizei height, GLint border,
                           GLenum format, GLenum type, const GLvoid *pixels) {
  TRACE_EVENT("(GLenum target = %d, GLint level = %d, GLint internalformat = "
              "%d, GLsizei width = %d, GLsizei height = %d, GLint border = %d, "
              "GLenum format = %d, GLenum type = %d, const GLvoid* pixels= %p)",
              target, level, internalformat, width, height, border, format,
              type, pixels);

  // lookup texture pointer
  Texture *tex = HandleToTexture(target);
  if (tex == NULL) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // only support level zero (no mipmaps) for now
  if (level != 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  int compos = -1;

  // only support 4 component or 1 component texture.
  if (internalformat == GL_RGBA) {
    // OK
    compos = 4;
  } else if (internalformat == GL_LUMINANCE) {
    // OK
    compos = 1;
  } else {
    return SetGLError(GL_INVALID_ENUM);
  }

  if (compos == -1) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // ensure width and height are valid
  if ((width < 0) || (height < 0)) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // borders are not supported yet
  if (border != 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  if (internalformat != format) {
    return SetGLError(GL_INVALID_ENUM);
  }

  if ((format == GL_RGBA) && (type == GL_UNSIGNED_BYTE)) {
    // OK
  } else if ((format == GL_RGBA) && (type == GL_FLOAT)) {
    // OK
  } else if ((format == GL_RGBA) && (type == GL_DOUBLE)) {
    // OK
  } else if ((format == GL_LUMINANCE) && (type == GL_UNSIGNED_BYTE)) {
    // OK
  } else if ((format == GL_LUMINANCE) && (type == GL_FLOAT)) {
    // OK
  } else if ((format == GL_LUMINANCE) && (type == GL_DOUBLE)) {
    // OK
  } else {
    return SetGLError(GL_INVALID_ENUM);
  }

  // upload the texture data. pixels may be null.
  tex->Data2D(pixels, width, height, compos, type);
}

void Context::glTexImage3D(GLenum target, GLint level, GLint internalformat,
                           GLsizei width, GLsizei height, GLsizei depth,
                           GLint border, GLenum format, GLenum type,
                           const GLvoid *pixels) {
  TRACE_EVENT("(GLenum target = %d, GLint level = %d, GLint internalformat = "
              "%d, GLsizei width = %d, GLsizei height = %d, GLsizei depth = "
              "%d, GLint border = %d, GLenum format = %d, GLenum type = %d, "
              "const GLvoid* pixels= %p)",
              target, level, internalformat, width, height, depth, border,
              format, type, pixels);

  // lookup texture pointer
  Texture *tex = HandleToTexture(target);
  if (tex == NULL) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // only support level zero (no mipmaps) for now
  if (level != 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // only 1, 3 or 4 component format are supported at this time.
  if ((internalformat == GL_LUMINANCE) || (internalformat == GL_RGB) || (internalformat == GL_RGBA)) {
    // OK
  } else {
    fprintf(stderr, "[LSGL] Unsupported internalformat for 3D texture.\n");
    assert(0);
    return SetGLError(GL_INVALID_ENUM);
  }

  // ensure width and height are valid
  if ((width < 0) || (height < 0) || (depth < 0)) {
    assert(0);
    return SetGLError(GL_INVALID_VALUE);
  }

  // borders are not supported yet
  if (border != 0) {
    assert(0);
    return SetGLError(GL_INVALID_VALUE);
  }

  // only support reading 1 or 3 component uchar or fp data for now
  if (((format == GL_LUMINANCE) && (type == GL_UNSIGNED_BYTE)) ||
      ((format == GL_RGB) && (type == GL_UNSIGNED_BYTE))) {
    // OK. uchar format
  } else if (((format == GL_LUMINANCE) && (type == GL_FLOAT)) ||
             ((format == GL_RGB) && (type == GL_FLOAT)) ||
             ((format == GL_RGBA) && (type == GL_FLOAT))) {
    // OK. float format
  } else if (((format == GL_LUMINANCE) && (type == GL_DOUBLE)) ||
             ((format == GL_RGB) && (type == GL_DOUBLE))) {
    // OK. double format
  } else {
    fprintf(stderr, "[LSGL] Unsupported format/type pair for 3D texture.\n");
    assert(0);
    return SetGLError(GL_INVALID_ENUM);
  }

  if (internalformat != format) {
    fprintf(stderr, "[LSGL] internalformat and format must be identical.\n");
    assert(0);
    return SetGLError(GL_INVALID_ENUM);
  }

  // upload the texture data. pixels may be null.
  int compos = -1;
  if (format == GL_LUMINANCE) {
    compos = 1;
  } else if (format == GL_RGB) {
    compos = 3;
  } else if (format == GL_RGBA) {
    compos = 4;
  }

  assert(compos != -1);
  tex->Data3D(pixels, width, height, depth, compos, type);
}

// Zero-copy version. User must retain texture memory until rendering finishes.
void Context::lsglTexImage3DPointer(GLenum target, GLint level,
                                    GLint internalformat, GLsizei width,
                                    GLsizei height, GLsizei depth, GLint border,
                                    GLenum format, GLenum type,
                                    const GLvoid *pixels) {
  TRACE_EVENT("(GLenum target = %d, GLint level = %d, GLint internalformat = "
              "%d, GLsizei width = %d, GLsizei height = %d, GLsizei depth = "
              "%d, GLint border = %d, GLenum format = %d, GLenum type = %d, "
              "const GLvoid* pixels= %p)",
              target, level, internalformat, width, height, depth, border,
              format, type, pixels);

  // lookup texture pointer
  Texture *tex = HandleToTexture(target);
  if (tex == NULL) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // only support level zero (no mipmaps) for now
  if (level != 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // only 1, 3 or 4 component format are supported at this time.
  if ((internalformat == GL_LUMINANCE) || (internalformat == GL_RGB) || (internalformat == GL_RGBA)) {
    // OK
  } else {
    fprintf(stderr, "[LSGL] Unsupported internalformat for 3D texture.\n");
    assert(0);
    return SetGLError(GL_INVALID_ENUM);
  }

  // ensure width and height are valid
  if ((width < 0) || (height < 0) || (depth < 0)) {
    assert(0);
    return SetGLError(GL_INVALID_VALUE);
  }

  // borders are not supported yet
  if (border != 0) {
    assert(0);
    return SetGLError(GL_INVALID_VALUE);
  }

  // only support reading 1 or 3 component uchar or fp data for now
  if (((format == GL_LUMINANCE) && (type == GL_UNSIGNED_BYTE)) ||
      ((format == GL_RGB) && (type == GL_UNSIGNED_BYTE))) {
    // OK. uchar format
  } else if (((format == GL_LUMINANCE) && (type == GL_FLOAT)) ||
             ((format == GL_RGB) && (type == GL_FLOAT)) ||
             ((format == GL_RGBA) && (type == GL_FLOAT))) {
    // OK. float format
  } else if (((format == GL_LUMINANCE) && (type == GL_DOUBLE)) ||
             ((format == GL_RGB) && (type == GL_DOUBLE))) {
    // OK. double format
  } else {
    fprintf(stderr, "[LSGL] Unsupported format/type pair for 3D texture.\n");
    assert(0);
    return SetGLError(GL_INVALID_ENUM);
  }

  if (internalformat != format) {
    fprintf(stderr, "[LSGL] internalformat and format must be identical.\n");
    assert(0);
    return SetGLError(GL_INVALID_ENUM);
  }

  // just retain a pointer of texture data.
  int compos = -1;
  if (format == GL_LUMINANCE) {
    compos = 1;
  } else if (format == GL_RGB) {
    compos = 3;
  } else if (format == GL_RGBA) {
    compos = 4;
  }
  tex->Retain3D(pixels, width, height, depth, compos, type);
}

void Context::glTexParameterf(GLenum target, GLenum pname, GLfloat param) {
  TRACE_EVENT("(GLenum target = %d, GLenum pname = %d, GLfloat param = %f\n",
              target, pname, param);
  assert(0);
}

void Context::glTexParameterfv(GLenum target, GLenum pname,
                               const GLfloat *params) {
  TRACE_EVENT(
      "(GLenum target = %d, GLenum pname = %d, const GLfloat* param = %p\n",
      target, pname, params);
  assert(0);
}

void Context::glTexParameteri(GLenum target, GLenum pname, GLint param) {
  TRACE_EVENT("(GLenum target = %d, GLenum pname = %d, GLint param = %d\n",
              target, pname, param);

  // lookup texture pointer
  Texture *tex = HandleToTexture(target);
  if (tex == NULL) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // @todo { Actually set texture sampler parameter. }

  switch (pname) {
  case GL_TEXTURE_MIN_FILTER:
    assert((param == GL_LINEAR) || (param == GL_NEAREST));
    break;
  case GL_TEXTURE_MAG_FILTER:
    assert((param == GL_LINEAR) || (param == GL_NEAREST));
    break;
  case GL_TEXTURE_WRAP_S:
    assert(param == GL_CLAMP_TO_EDGE);
    break;
  case GL_TEXTURE_WRAP_T:
    assert(param == GL_CLAMP_TO_EDGE);
    break;
  case GL_TEXTURE_WRAP_R:
    assert(param == GL_CLAMP_TO_EDGE);
    break;
  default:
    assert(0 && "Unsupported parameter.");
  }
}

void Context::glTexParameteriv(GLenum target, GLenum pname,
                               const GLint *params) {
  TRACE_EVENT(
      "(GLenum target = %d, GLenum pname = %d, const GLint* params = %p\n",
      target, pname, params);
  assert(0);
}

void Context::glTexSubImage2D(GLenum target, GLint level, GLint xoffset,
                              GLint yoffset, GLsizei width, GLsizei height,
                              GLenum format, GLenum type,
                              const GLvoid *pixels) {
  TRACE_EVENT("(GLenum target = %d, GLint level = %d, GLint xoffiset = %d, "
              "GLint yoffset = %d, GLsizei width = %d, GLsizei height = %d, "
              "GLenum format = %d, GLenum type = %d, const GLvoid* pixels= %p)",
              target, level, xoffset, yoffset, width, height, format, type,
              pixels);

  assert(0 && "glTexSubImage2D is not supported yet");
}

void Context::glDeleteTextures(GLsizei n, const GLuint *textures) {
  TRACE_EVENT("(GLsizei n = %d, const GLuint* textures = %p)", n, textures);
  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(textures && "textures are null");

  for (int i = 0; i < n; i++) {
    resourceManager_.DeleteTexture(textures[i]);
  }
}

GLboolean Context::glIsTexture(GLuint texture) {
  TRACE_EVENT("(GLuint texture = %d)", texture);
  return resourceManager_.IsValidTexture(texture);
}

void Context::glCompressedTexImage2D(GLenum target, GLint level,
                                     GLenum internalformat, GLsizei width,
                                     GLsizei height, GLint border,
                                     GLsizei imageSize, const GLvoid *data) {
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glCompressedTexSubImage2D(GLenum target, GLint level,
                                        GLint xoffset, GLint yoffset,
                                        GLsizei width, GLsizei height,
                                        GLenum format, GLsizei imageSize,
                                        const GLvoid *data) {
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glCopyTexImage2D(GLenum target, GLint level,
                               GLenum internalformat, GLint x, GLint y,
                               GLsizei width, GLsizei height, GLint border) {
  TRACE_EVENT("(GLenum target = %d, GLint level = %d, GLenum internalformat = "
              "%d, GLint x = %d, GLint y = %d, GLsizei width = %d, GLsizei "
              "height = %d, GLint border = %d)",
              target, level, internalformat, x, y, width, height, border);

  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glCopyTexSubImage2D(GLenum target, GLint level, GLint xoffset,
                                  GLint yoffset, GLint x, GLint y,
                                  GLsizei width, GLsizei height) {
  TRACE_EVENT("(GLenum target = %d, GLint level = %d, GLint xoffiset = %d, "
              "GLint yoffset = %d, GLint x = %d, GLint y = %d, GLsizei width = "
              "%d, GLsizei height = %d)",
              target, level, xoffset, yoffset, x, y, width, height);

  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetTexParameterfv(GLenum target, GLenum pname,
                                  GLfloat *params) {
  TRACE_EVENT("(GLenum target = %d, GLenum pname = %d, GLfloat* params = %p\n",
              target, pname, params);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetTexParameteriv(GLenum target, GLenum pname, GLint *params) {
  TRACE_EVENT("(GLenum target = %d, GLenum pname = %d, GLint* params = %p\n",
              target, pname, params);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}
