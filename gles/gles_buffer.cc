/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cstring>

#include "gles_context.h"
#include "gles_resource_manager.h"

using namespace lsgl;

Buffer *Context::HandleToBuffer(GLenum target) {
  // get currently bound handle
  GLuint handle;
  switch (target) {
  case GL_ELEMENT_ARRAY_BUFFER:
    handle = state_.elementArrayBuffer;
    break;

  case GL_ARRAY_BUFFER:
    handle = state_.arrayBuffer;
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

  // get buffer pointer from handle
  Buffer *buf = resourceManager_.GetBuffer(handle);
  if (buf == NULL) {
    SetGLError(GL_INVALID_OPERATION);
    return NULL;
  }

  return buf;
}

void Context::glGenBuffers(GLsizei n, GLuint *buffers) {
  TRACE_EVENT("(GLsizei n = %d, GLuint* buffers = %p)", n, buffers);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(buffers && "buffers are null");

  for (int i = 0; i < n; i++) {
    buffers[i] = resourceManager_.CreateBuffer();
  }
}

void Context::glGenFramebuffers(GLsizei n, GLuint *framebuffers) {
  TRACE_EVENT("(GLsizei n = %d, GLuint* framebuffers = %p)", n, framebuffers);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(framebuffers && "framebuffers are null");

  for (int i = 0; i < n; i++) {
    framebuffers[i] = resourceManager_.CreateFramebuffer();
  }
}

void Context::glGenRenderbuffers(GLsizei n, GLuint *renderbuffers) {
  TRACE_EVENT("(GLsizei n = %d, GLuint* renderbuffers = %p)", n, renderbuffers);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(renderbuffers && "renderbuffers are null");

  for (int i = 0; i < n; i++) {
    renderbuffers[i] = resourceManager_.CreateRenderbuffer();
  }
}

void Context::glGetBufferParameteriv(GLenum target, GLenum pname,
                                     GLint *params) {
  TRACE_EVENT("(GLenum target = %d, GLenum pname = %d, GLint* params = %p)",
              target, pname, params);

  // lookup buffer pointer
  Buffer *buf = HandleToBuffer(target);
  if (buf == NULL) {
    return;
  }

  // return value based on desired parameter
  switch (pname) {
  case GL_BUFFER_SIZE:
    *params = buf->GetSize();
    break;

  case GL_BUFFER_USAGE:
    *params = buf->GetUsage();
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glBindBuffer(GLenum target, GLuint buffer) {
  TRACE_EVENT("(GLenum target = %d, GLuint buffer = %d)", target, buffer);

  // store new buffer target handle
  switch (target) {
  case GL_ELEMENT_ARRAY_BUFFER:
    state_.elementArrayBuffer = buffer;
    break;

  case GL_ARRAY_BUFFER:
    state_.arrayBuffer = buffer;
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glBindFramebuffer(GLenum target, GLuint framebuffer) {
  TRACE_EVENT("(GLenum target = %d, GLuint framebuffer = %d)", target,
              framebuffer);

  if (target != GL_FRAMEBUFFER) {
    return SetGLError(GL_INVALID_ENUM);
  }

  state_.frameBuffer = framebuffer;
}

void Context::glBindRenderbuffer(GLenum target, GLuint renderbuffer) {
  TRACE_EVENT("(GLenum target = %d, GLuint renderbuffer = %d)", target,
              renderbuffer);

  if (target != GL_RENDERBUFFER) {
    return SetGLError(GL_INVALID_ENUM);
  }

  state_.renderBuffer = renderbuffer;
}

void Context::glDeleteBuffers(GLsizei n, const GLuint *buffers) {
  TRACE_EVENT("(GLsizei n = %d, const GLuint* buffers = %p)", n, buffers);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(buffers && "buffers are null");

  for (int i = 0; i < n; i++) {
    // invalidate the buffer in our mesh builder first, then delete it from the
    // resource manager
    Buffer *buf = resourceManager_.GetBuffer(buffers[i]);
    if (buf != NULL) {
      meshBuilder_.Invalidate(buf);
    }
    resourceManager_.DeleteBuffer(buffers[i]);
  }
}

void Context::glDeleteFramebuffers(GLsizei n, const GLuint *framebuffers) {
  TRACE_EVENT("(GLsizei n = %d, const GLuint* framebuffers = %p)", n,
              framebuffers);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(framebuffers && "framebuffers are null");

  for (int i = 0; i < n; i++) {
    resourceManager_.DeleteFramebuffer(framebuffers[i]);
  }
}

void Context::glDeleteRenderbuffers(GLsizei n, const GLuint *renderbuffers) {
  TRACE_EVENT("(GLsizei n = %d, const GLuint* renderbuffers = %p)", n,
              renderbuffers);

  if (n < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (n == 0) {
    return;
  }

  assert(renderbuffers && "renderbuffers are null");

  for (int i = 0; i < n; i++) {
    resourceManager_.DeleteRenderbuffer(renderbuffers[i]);
  }
}

GLboolean Context::glIsBuffer(GLuint buffer) {
  TRACE_EVENT("(GLuint buffer = %d)", buffer);

  return resourceManager_.IsValidBuffer(buffer);
}

GLboolean Context::glIsFramebuffer(GLuint framebuffer) {
  TRACE_EVENT("(GLuint framebuffer = %d)", framebuffer);

  return resourceManager_.IsValidFramebuffer(framebuffer);
}

GLboolean Context::glIsRenderbuffer(GLuint renderbuffer) {
  TRACE_EVENT("(GLuint renderbuffer = %d)", renderbuffer);

  return resourceManager_.IsValidRenderbuffer(renderbuffer);
}

void Context::glRenderbufferStorage(GLenum target, GLenum internalformat,
                                    GLsizei width, GLsizei height) {
  TRACE_EVENT("(GLeunm target = %d, GLenum internalformat = %d, GLsizei width "
              "= %d, GLsizei height = %d)",
              target, internalformat, width, height);

  // ensure target is set to render buffer
  if (target != GL_RENDERBUFFER) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // ensure width and height are valid
  if ((width < 0) || (height < 0)) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // lookup bound render buffer pointer
  Renderbuffer *rb = resourceManager_.GetRenderbuffer(state_.renderBuffer);
  if (rb == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // allocate render buffer memory
  bool res = false;
  switch (internalformat) {
  case GL_RGBA:
    res = rb->Allocate(width, height, 4, internalformat);
    break;

  case GL_RGBA8_OES:
    res = rb->Allocate(width, height, 4, internalformat);
    break;

  case GL_RGBA32F_EXT:
    res = rb->Allocate(width, height, 16, internalformat);
    break;

  // case GL_DEPTH_COMPONENT16:
  //    res = rb->Allocate(width, height, 2, internalformat);
  //    break;

  case GL_DEPTH_COMPONENT:
    // Treat as 32bit float depth.
    res = rb->Allocate(width, height, 4, internalformat);
    break;

  case GL_DEPTH_COMPONENT32_OES:
    res = rb->Allocate(width, height, 4, internalformat);
    break;

  case GL_STENCIL_INDEX8:
    res = rb->Allocate(width, height, 1, internalformat);
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }

  // check for out of memory condition
  if (res == false) {
    SetGLError(GL_OUT_OF_MEMORY);
  }
}

void Context::glFramebufferRenderbuffer(GLenum target, GLenum attachment,
                                        GLenum renderbuffertarget,
                                        GLuint renderbuffer) {
  TRACE_EVENT("(GLeunm target = %d, GLenum attachment = %d, GLenum "
              "renderbuffertarget = %d, GLuint renderbuffer = %d)",
              target, attachment, renderbuffertarget, renderbuffer);

  // ensure target is set to frame buffer and render buffer target is set to
  // render buffer
  if ((target != GL_FRAMEBUFFER) || (renderbuffertarget != GL_RENDERBUFFER)) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // lookup bound frame buffer pointer
  Framebuffer *fb = resourceManager_.GetFramebuffer(state_.frameBuffer);
  if (fb == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // if attaching, lookup the specified render buffer pointer
  Renderbuffer *rb = NULL;
  if (renderbuffer != 0) {
    rb = resourceManager_.GetRenderbuffer(renderbuffer);
    if (rb == NULL) {
      return SetGLError(GL_INVALID_OPERATION);
    }
  }

  // attach/detach this target to the frame buffer
  switch (attachment) {
  case GL_COLOR_ATTACHMENT0:
    fb->AttachColorBuffer(0, rb);
    break;

  case GL_DEPTH_ATTACHMENT:
    fb->AttachDepthBuffer(rb);
    break;

  case GL_STENCIL_ATTACHMENT:
    fb->AttachStencilBuffer(rb);
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glFramebufferTexture2D(GLenum target, GLenum attachment,
                                     GLenum textarget, GLuint texture,
                                     GLint level) {
  TRACE_EVENT("(GLeunm target = %d, GLenum attachment = %d, GLenum textarget = "
              "%d, GLuint teture = %d, GLint level = %d)",
              target, attachment, textarget, texture, level);

  assert(0 && "TODO");
}

void Context::glBufferData(GLenum target, GLsizeiptr size, const GLvoid *data,
                           GLenum usage) {
  TRACE_EVENT("(GLeunm target = %d, GLsizeiptr size = %ld, const GLvoid* data "
              "= %p, GLeunum usage = %d)",
              target, size, data, usage);

  // ensure size is valid (allow zero to pass as it is used to free the buffer
  // memory)
  if (size < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // lookup buffer pointer
  Buffer *buf = HandleToBuffer(target);
  if (buf == NULL) {
    return;
  }

  // store data into buffer
  buf->Data(size, data, usage);
}

void Context::glBufferSubData(GLenum target, GLintptr offset, GLsizeiptr size,
                              const GLvoid *data) {
  TRACE_EVENT("(GLeunm target = %d, GLintptr offset = %ld, GLsizeiptr size = "
              "%ld, const GLvoid* data = %p)",
              target, offset, size, data);

  // ensure size is valid
  if (size < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (size == 0) {
    return;
  }

  // lookup buffer pointer
  Buffer *buf = HandleToBuffer(target);
  if (buf == NULL) {
    return;
  }

  // store data into buffer
  buf->SubData(offset, size, data);
}

GLenum Context::glCheckFramebufferStatus(GLenum target) {
  TRACE_EVENT("(GLeunm target = %d)", target);

  // [4.4] If CheckFramebufferStatus generates an error, 0 is returned.

  // ensure the target is set to frame buffer
  if (target != GL_FRAMEBUFFER) {
    SetGLError(GL_INVALID_ENUM);
    return 0;
  }

  // only test if the currently bound framebuffer is not zero, as the window
  // system provided framebuffer is always complete
  if (state_.frameBuffer != 0) {
    // lookup the bound framebuffer pointer
    Framebuffer *fb = resourceManager_.GetFramebuffer(state_.frameBuffer);
    if (fb == NULL) {
      return 0; // @fixme { Set GL error? }
    }

    // fail if no render buffers are attached to this framebuffer
    if (fb->HasAttachment() == false) {
      return GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT;
    }

    // ensure all attached render buffers have the same width and height
    if (fb->HasValidSize() == false) {
      return GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS;
    }
  }

  // everything looks good!
  return GL_FRAMEBUFFER_COMPLETE;
}

void Context::glGetFramebufferAttachmentParameteriv(GLenum target,
                                                    GLenum attachment,
                                                    GLenum pname,
                                                    GLint *params) {
  TRACE_EVENT("(GLeunm target = %d, GLenum attachment = %d, GLenum pname = %d, "
              "GLint* params = %p)",
              target, attachment, pname, params);

  // ensure the target is set to frame buffer
  if (target != GL_FRAMEBUFFER) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // lookup the currently bound framebuffer
  Framebuffer *fb = GetCurrentFrameBuffer();
  if (fb == NULL) {
    return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
  }

  // lookup the desired renderbuffer
  Renderbuffer *rb;
  switch (attachment) {
  case GL_COLOR_ATTACHMENT0:
    rb = fb->GetColorBuffer(0);
    break;

  case GL_DEPTH_ATTACHMENT:
    rb = fb->GetDepthBuffer();
    break;

  case GL_STENCIL_ATTACHMENT:
    rb = fb->GetStencilBuffer();
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }

  // store the return parameter based on the pname value
  // TODO: texture bound to framebuffer support
  switch (pname) {
  case GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE:
    if (rb != NULL) {
      *params = GL_RENDERBUFFER;
    } else {
      *params = GL_NONE;
    }
    break;

  case GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME:
    if (rb != NULL) {
      *params = resourceManager_.GetRenderbufferHandle(rb);
    } else {
      return SetGLError(GL_INVALID_ENUM);
    }
    break;

  case GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL:
    return SetGLError(GL_INVALID_ENUM);

  case GL_FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE:
    return SetGLError(GL_INVALID_ENUM);

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glGetRenderbufferParameteriv(GLenum target, GLenum pname,
                                           GLint *params) {
  TRACE_EVENT("(GLeunm target = %d, GLenum pname = %d, GLint* params = %p)",
              target, pname, params);

  // ensure the target is set to renderbuffer
  if (target != GL_RENDERBUFFER) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // fail if the currently bound renderbuffer is the system-allocated one
  if (state_.renderBuffer == 0) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // lookup the currently bound renderbuffer
  Renderbuffer *rb = resourceManager_.GetRenderbuffer(state_.renderBuffer);
  if (rb == NULL) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // store the return parameter based on the pname value
  switch (pname) {
  case GL_RENDERBUFFER_WIDTH:
    *params = rb->GetWidth();
    break;

  case GL_RENDERBUFFER_HEIGHT:
    *params = rb->GetHeight();
    break;

  case GL_RENDERBUFFER_INTERNAL_FORMAT:
    *params = rb->GetFormat();
    break;

  case GL_RENDERBUFFER_RED_SIZE:
  case GL_RENDERBUFFER_GREEN_SIZE:
  case GL_RENDERBUFFER_BLUE_SIZE:
  case GL_RENDERBUFFER_ALPHA_SIZE:
  case GL_RENDERBUFFER_DEPTH_SIZE:
  case GL_RENDERBUFFER_STENCIL_SIZE:
    *params = rb->GetBitsPerPixel();
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height,
                           GLenum format, GLenum type, GLvoid *data) {
  TRACE_EVENT("(GLint x = %d, GLint y = %d, GLsizei width = %d, GLsizei height "
              "= %d, GLenum format = %d, GLenm type = %d, GLvoid* data = %p)",
              x, y, width, height, format, type, data);

  bool readDepth = false;

  if (data == NULL) {
    fprintf(stderr, "[LSGL] ERR: Null input for data.\n");
    return SetGLError(GL_INVALID_VALUE);
  }

  // ensure width and height are valid
  if ((width < 0) || (height < 0)) {
    return SetGLError(GL_INVALID_VALUE);
  }

  //
  // format/type check
  //
  if (((format == GL_DEPTH_COMPONENT) && (type == GL_FLOAT)) ||
      ((format == GL_DEPTH_COMPONENT32_OES) && (type == GL_FLOAT))) {
    // Depth read.
    readDepth = true;
  } else {
    // Colo read

    if ((format == GL_RGBA) && (type == GL_FLOAT)) {
      // OK. float HDR buffer read
    } else {
      // only support 32-bit RGBA reading from the color buffer for now
      if ((format != GL_RGBA) || (type != GL_UNSIGNED_BYTE)) {
        fprintf(stderr, "[LSGL] ERR: Invalid framebuffer format.\n");
        return SetGLError(GL_INVALID_ENUM);
      }
    }
  }

  // lookup the currently bound framebuffer
  Framebuffer *fb = GetCurrentFrameBuffer();
  if (fb == NULL) {
    fprintf(stderr, "[LSGL] ERR: Cannot find currently bound framebuffer.\n");
    return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
  }

  Renderbuffer *rb = NULL;
  if (readDepth) {
    rb = fb->GetDepthBuffer();
  } else {
    rb = fb->GetColorBuffer(0);
  }

  if (rb == NULL) {
    fprintf(stderr, "[LSGL] ERR: Cannot find color buffer.\n");
    return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
  }

  // ensure rendering has been completed first
  glFinish();

  if ((x == 0) && (y == 0) && (width == rb->GetWidth()) &&
      (height == rb->GetHeight())) {
    // check if we can do an optimized copy
    if ((rb->GetFormat() == GL_RGBA32F_EXT) && (format == GL_RGBA) &&
        (type == GL_FLOAT)) {
      memcpy(data, rb->GetBuffer(), rb->GetSize());
    } else if ((rb->GetFormat() == GL_RGBA) && (format == GL_RGBA) &&
               (type == GL_UNSIGNED_BYTE)) {
      memcpy(data, rb->GetBuffer(), rb->GetSize());
    } else if ((rb->GetFormat() == GL_RGBA8_OES) && (format == GL_RGBA) &&
               (type == GL_UNSIGNED_BYTE)) {
      memcpy(data, rb->GetBuffer(), rb->GetSize());
    } else if ((rb->GetFormat() == GL_DEPTH_COMPONENT32_OES) &&
               (format == GL_DEPTH_COMPONENT) && (type == GL_FLOAT)) {
      memcpy(data, rb->GetBuffer(), rb->GetSize());
    } else if ((rb->GetFormat() == GL_DEPTH_COMPONENT) &&
               (format == GL_DEPTH_COMPONENT) && (type == GL_FLOAT)) {
      memcpy(data, rb->GetBuffer(), rb->GetSize());
    } else {
      fprintf(stderr, "[LSGL] Unsupported internalformat & format/type pair. "
                      "internalformat = %d, format = %d, type = %d\n",
              rb->GetFormat(), format, type);
      assert(0 && "TODO");
    }

  } else {
    // @todo { Reading partial region }
    printf("rb->width %d, rb->height %d\n", rb->GetWidth(), rb->GetHeight());
    assert(0 && "TODO");
  }
}

void Context::glPixelStorei(GLenum pname, GLint param) {
  TRACE_EVENT("(GLenum pname = %d, GLint param = %d)", pname, param);

  // ensure alignment is valid
  if ((param != 1) && (param != 2) && (param != 4) && (param != 8))
    return SetGLError(GL_INVALID_VALUE);

  // store new alignment based on parameter name
  switch (pname) {
  case GL_PACK_ALIGNMENT:
    state_.packAlignment = param;
    break;

  case GL_UNPACK_ALIGNMENT:
    state_.unpackAlignment = param;
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::lsglBufferDataPointer(GLenum target, GLsizeiptr size,
                                    const GLvoid *data, GLenum usage) {
  TRACE_EVENT("(GLeunm target = %d, GLsizeiptr size = %ld, const GLvoid* data "
              "= %p, GLeunum usage = %d)",
              target, size, data, usage);

  // ensure size is valid (allow zero to pass as it is used to free the buffer
  // memory)
  if (size < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // lookup buffer pointer
  Buffer *buf = HandleToBuffer(target);
  if (buf == NULL) {
    return;
  }

  // store data into buffer
  buf->Retain(size, data, usage);
}

void Context::lsglInvalidateBuffer(GLenum target) {
  TRACE_EVENT("(GLenum target = %d)", target);

  Buffer *buf = HandleToBuffer(target);
  if (buf != NULL) {
    meshBuilder_.Invalidate(buf);
  }
}
