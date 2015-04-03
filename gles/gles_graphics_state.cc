/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cstring>
#include <cstdarg>

#include "GLES2/gl2.h"
#include "gles_context.h"

using namespace lsgl;

void Context::glClear(GLbitfield mask) {
  Renderbuffer *rb;

  // return an error if any bits are on in the mask that we don't recognize
  if ((mask & ~(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
                GL_STENCIL_BUFFER_BIT)) != 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // lookup the currently bound framebuffer
  Framebuffer *fb = GetCurrentFrameBuffer();
  if (fb == NULL) {
    return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
  }

  // clear color buffer if desired
  rb = fb->GetColorBuffer(0);
  if ((mask & GL_COLOR_BUFFER_BIT) && (rb != NULL)) {

    if ((rb->GetFormat() == GL_RGBA) || (rb->GetFormat() == GL_RGBA8_OES)) {

      // convert floating point color components to 8-bit depth
      const uint8_t red = (uint8_t)(state_.colorClearValue_.red * 255.0f);
      const uint8_t green = (uint8_t)(state_.colorClearValue_.green * 255.0f);
      const uint8_t blue = (uint8_t)(state_.colorClearValue_.blue * 255.0f);
      const uint8_t alpha = (uint8_t)(state_.colorClearValue_.alpha * 255.0f);

      // check if we can use an optimized clear
      if ((red == green) && (green == blue) && (blue == alpha)) {
        memset(rb->GetBuffer(), red, rb->GetSize());
      } else {

        // build combined 32-bit color value and clear the buffer with it
        // const uint32_t color = ((uint32_t)alpha << 24) | ((uint32_t)red <<
        // 16)
        //|
        //((uint32_t)green << 8) | blue;
        uint32_t color;
        uint8_t *cptr = reinterpret_cast<uint8_t *>(&color);
        cptr[0] = red;
        cptr[1] = green;
        cptr[2] = blue;
        cptr[3] = alpha;

        uint32_t *ptr = reinterpret_cast<uint32_t *>(rb->GetBuffer());
        GLuint sized4 = rb->GetSize() >> 2;

        while (sized4--) {
          *ptr++ = color;
        }
      }

    } else if (rb->GetFormat() == GL_RGBA32F_EXT) {

      float *ptr = reinterpret_cast<float *>(rb->GetBuffer());

      const float red = state_.colorClearValue_.red;
      const float green = state_.colorClearValue_.green;
      const float blue = state_.colorClearValue_.blue;
      const float alpha = state_.colorClearValue_.alpha;

      for (int i = 0; i < rb->GetWidth() * rb->GetHeight(); i++) {
        ptr[4 * i + 0] = red;
        ptr[4 * i + 1] = green;
        ptr[4 * i + 2] = blue;
        ptr[4 * i + 3] = alpha;
      }

    } else {

      fprintf(stderr, "[LSGL] glClear for unsupported buffer format.");
      assert(0);
    }
  }

  // clear depth buffer if desired
  rb = fb->GetDepthBuffer();
  if ((mask & GL_DEPTH_BUFFER_BIT) && (rb != NULL)) {
    // Assume 32bit float depth buffer
    float depth = state_.depthClearValue_;
    float *ptr = reinterpret_cast<float *>(rb->GetBuffer());
    GLuint sized4 = rb->GetSize() >> 2;

    while (sized4--) {
      *ptr++ = depth;
    }
  }

  // clear stencil buffer if desired
  rb = fb->GetStencilBuffer();
  if ((mask & GL_STENCIL_BUFFER_BIT) && (rb != NULL)) {
    memset(rb->GetBuffer(), state_.stencilClearValue_, rb->GetSize());
  }

  // mark our frame as dirty
  dirtyFrame_ = true;
}

void Context::glClearColor(GLclampf red, GLclampf green, GLclampf blue,
                           GLclampf alpha) {
  state_.colorClearValue_.red = red;
  state_.colorClearValue_.green = green;
  state_.colorClearValue_.blue = blue;
  state_.colorClearValue_.alpha = alpha;
}

void Context::glClearDepthf(GLclampf depth) { state_.depthClearValue_ = depth; }

void Context::glClearStencil(GLint s) { state_.stencilClearValue_ = s; }

void Context::SetState(GLenum cap, bool flag) {
  switch (cap) {
  case GL_BLEND:
    state_.blend = flag;
    break;

  case GL_CULL_FACE:
    state_.cullFace = flag;
    break;

  case GL_DEPTH_TEST:
    state_.depthTest = flag;
    break;

  case GL_DITHER:
    state_.dither = flag;
    break;

  case GL_POLYGON_OFFSET_FILL:
    state_.polygonOffsetFill = flag;
    break;

  case GL_SAMPLE_ALPHA_TO_COVERAGE:
    state_.sampleAlphaToCoverage = flag;
    break;

  case GL_SAMPLE_COVERAGE:
    state_.sampleCoverage = flag;
    break;

  case GL_SCISSOR_TEST:
    state_.scissorTest = flag;
    break;

  case GL_STENCIL_TEST:
    state_.stencilTest = flag;
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glEnable(GLenum cap) { SetState(cap, true); }

void Context::glDisable(GLenum cap) { SetState(cap, false); }

GLboolean Context::glIsEnabled(GLenum cap) {
  switch (cap) {
  case GL_BLEND:
    return state_.blend;

  case GL_CULL_FACE:
    return state_.cullFace;

  case GL_DEPTH_TEST:
    return state_.depthTest;

  case GL_DITHER:
    return state_.dither;

  case GL_POLYGON_OFFSET_FILL:
    return state_.polygonOffsetFill;

  case GL_SAMPLE_ALPHA_TO_COVERAGE:
    return state_.sampleAlphaToCoverage;

  case GL_SAMPLE_COVERAGE:
    return state_.sampleCoverage;

  case GL_SCISSOR_TEST:
    return state_.scissorTest;

  case GL_STENCIL_TEST:
    return state_.stencilTest;
  }

  SetGLError(GL_INVALID_ENUM);
  return false;
}

void Context::glLineWidth(GLfloat width) { state_.lineWidth = width; }

GLenum Context::glGetError(void) {
  GLenum res = state_.lastError;
  state_.lastError = GL_NO_ERROR;
  return res;
}

#define GLINT_TO_GLBOOLEAN(val) ((GLboolean)((val == 0) ? GL_FALSE : GL_TRUE))
#define GLFLOAT_TO_GLBOOLEAN(val)                                              \
  ((GLboolean)((val == 0.0) ? GL_FALSE : GL_TRUE))
#define GLBOOLEAN_TO_GLINT(val) ((GLint)val)
#define GLFLOAT_TO_GLINT(val) ((GLint)val)
#define GLFLOAT_TO_GLINT_SCALE(val)                                            \
  ((GLint)(val * (std::numeric_limits<int>::max)()))
#define GLBOOLEAN_TO_GLFLOAT(val) ((GLfloat)val)
#define GLINT_TO_GLFLOAT(val) ((GLfloat)val)

static void _lsglStoreb(GLboolean *dest, unsigned int count, ...) {
  // copy GLboolean to GLboolean
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = (GLboolean)va_arg(args, int);
  va_end(args);
}

static void _lsglStoreb(GLint *dest, unsigned int count, ...) {
  // convert GLboolean to GLint
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLBOOLEAN_TO_GLINT(va_arg(args, int));
  va_end(args);
}

static void _lsglStoreb(GLfloat *dest, unsigned int count, ...) {
  // convert GLboolean to GLfloat
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLBOOLEAN_TO_GLFLOAT(va_arg(args, int));
  va_end(args);
}

static void _lsglStorei(GLboolean *dest, unsigned int count, ...) {
  // convert GLint to GLboolean
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLINT_TO_GLBOOLEAN(va_arg(args, int));
  va_end(args);
}

static void _lsglStorei(GLint *dest, unsigned int count, ...) {
  // copy GLint to GLint
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = (GLint)va_arg(args, int);
  va_end(args);
}

static void _lsglStorei(GLfloat *dest, unsigned int count, ...) {
  // convert GLint to GLfloat
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLINT_TO_GLFLOAT(va_arg(args, int));
  va_end(args);
}

static void _lsglStoref(GLboolean *dest, unsigned int count, ...) {
  // convert GLfloat to GLboolean
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLFLOAT_TO_GLBOOLEAN(va_arg(args, double));
  va_end(args);
}

static void _lsglStoref(GLint *dest, unsigned int count, ...) {
  // convert GLfloat to GLint
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLFLOAT_TO_GLINT(va_arg(args, double));
  va_end(args);
}

static void _lsglStoref(GLfloat *dest, unsigned int count, ...) {
  // copy GLfloat to GLfloat
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = (GLfloat)va_arg(args, double);
  va_end(args);
}

static void _lsglStoreScalef(GLboolean *dest, unsigned int count, ...) {
  // convert GLfloat to GLboolean (not actually scaled, as it doesn't have any
  // affect on a boolean value)
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLFLOAT_TO_GLBOOLEAN(va_arg(args, double));
  va_end(args);
}

static void _lsglStoreScalef(GLint *dest, unsigned int count, ...) {
  // convert GLfloat to GLint, scaled
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = GLFLOAT_TO_GLINT_SCALE(va_arg(args, double));
  va_end(args);
}

static void _lsglStoreScalef(GLfloat *dest, unsigned int count, ...) {
  // copy GLfloat to GLfloat
  va_list args;
  va_start(args, count);
  while (count--)
    *dest++ = (GLfloat)va_arg(args, double);
  va_end(args);
}

template <typename T> void Context::GetValue(GLenum pname, T *params) {
  switch (pname) {
  case GL_ACTIVE_TEXTURE:
    _lsglStorei(params, 1, state_.activeSampler);
    break;

  case GL_ALIASED_LINE_WIDTH_RANGE:
    _lsglStoref(params, 2, kMinPointSize, kMaxPointSize);
    break;

  case GL_ALIASED_POINT_SIZE_RANGE:
    _lsglStoref(params, 2, kMinPointSize, kMaxPointSize);
    break;

  case GL_ALPHA_BITS:
  case GL_BLUE_BITS:
  case GL_GREEN_BITS:
  case GL_RED_BITS: {
    Framebuffer *fb = GetCurrentFrameBuffer();
    if (fb == NULL) {
      return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
    }

    Renderbuffer *rb = fb->GetColorBuffer(0);
    if (rb == NULL) {
      _lsglStorei(params, 1, 0);
    } else {
      _lsglStorei(params, 1, rb->GetBitsPerPixel());
    }
  } break;

  case GL_ARRAY_BUFFER_BINDING:
    _lsglStorei(params, 1, state_.arrayBuffer);
    break;

  case GL_BLEND:
    _lsglStoreb(params, 1, state_.blend);
    break;

  case GL_BLEND_COLOR:
    _lsglStoreScalef(params, 4, state_.blendColor.red, state_.blendColor.green,
                     state_.blendColor.blue, state_.blendColor.alpha);
    break;

  case GL_BLEND_DST_ALPHA:
    _lsglStorei(params, 1, state_.destBlendAlpha);
    break;

  case GL_BLEND_DST_RGB:
    _lsglStorei(params, 1, state_.destBlendRGB);
    break;

  case GL_BLEND_EQUATION_ALPHA:
    _lsglStorei(params, 1, state_.blendEquationAlpha);
    break;

  case GL_BLEND_EQUATION_RGB:
    _lsglStorei(params, 1, state_.blendEquationRGB);
    break;

  case GL_BLEND_SRC_ALPHA:
    _lsglStorei(params, 1, state_.sourceBlendAlpha);
    break;

  case GL_BLEND_SRC_RGB:
    _lsglStorei(params, 1, state_.sourceBlendRGB);
    break;

  case GL_COLOR_CLEAR_VALUE:
    _lsglStoreScalef(
        params, 4, state_.colorClearValue_.red, state_.colorClearValue_.green,
        state_.colorClearValue_.blue, state_.colorClearValue_.alpha);
    break;

  case GL_COLOR_WRITEMASK:
    _lsglStoreb(params, 4, state_.colorMaskRed, state_.colorMaskGreen,
                state_.colorMaskBlue, state_.colorMaskAlpha);
    break;

  case GL_COMPRESSED_TEXTURE_FORMATS:
    // compressed textures are not currently supported (we return zero for
    // GL_NUM_COMPRESSED_TEXTURE_FORMATS), so do nothing here
    break;

  case GL_CULL_FACE:
    _lsglStoreb(params, 1, state_.cullFace);
    break;

  case GL_CULL_FACE_MODE:
    _lsglStorei(params, 1, state_.cullMode);
    break;

  case GL_CURRENT_PROGRAM:
    _lsglStorei(params, 1, state_.currentProgram);
    break;

  case GL_DEPTH_BITS: {
    Framebuffer *fb = GetCurrentFrameBuffer();
    if (fb == NULL) {
      return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
    }

    Renderbuffer *rb = fb->GetDepthBuffer();
    if (rb == NULL) {
      _lsglStorei(params, 1, 0);
    } else {
      _lsglStorei(params, 1, rb->GetBitsPerPixel());
    }
  } break;

  case GL_DEPTH_CLEAR_VALUE:
    _lsglStoreScalef(params, 1, state_.depthClearValue_);
    break;

  case GL_DEPTH_FUNC:
    _lsglStorei(params, 1, state_.depthFunc);
    break;

  case GL_DEPTH_RANGE:
    _lsglStoreScalef(params, 2, state_.nearDepth, state_.farDepth);
    break;

  case GL_DEPTH_TEST:
    _lsglStorei(params, 1, state_.depthTest);
    break;

  case GL_DEPTH_WRITEMASK:
    _lsglStoreb(params, 1, state_.depthMask);
    break;

  case GL_DITHER:
    _lsglStoreb(params, 1, state_.dither);
    break;

  case GL_ELEMENT_ARRAY_BUFFER_BINDING:
    _lsglStorei(params, 1, state_.elementArrayBuffer);
    break;

  case GL_FRAMEBUFFER_BINDING:
    _lsglStorei(params, 1, state_.frameBuffer);
    break;

  case GL_FRONT_FACE:
    _lsglStorei(params, 1, state_.frontFace);
    break;

  case GL_GENERATE_MIPMAP_HINT:
    _lsglStorei(params, 1, state_.generateMipmapHint);
    break;

  case GL_IMPLEMENTATION_COLOR_READ_FORMAT:
    _lsglStorei(params, 1, GL_RGBA);
    break;

  case GL_IMPLEMENTATION_COLOR_READ_TYPE:
    _lsglStorei(params, 1, GL_UNSIGNED_BYTE);
    break;

  case GL_LINE_WIDTH:
    _lsglStoref(params, 1, state_.lineWidth);
    break;

  case GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS:
    _lsglStorei(params, 1, kMaxTextureImageUnits);
    break;

  case GL_MAX_CUBE_MAP_TEXTURE_SIZE:
    _lsglStorei(params, 1, kMaxTextureWidthHeight);
    break;

  case GL_MAX_FRAGMENT_UNIFORM_VECTORS:
    _lsglStorei(params, 1, kMaxFragmentUniformVectors);
    break;

  case GL_MAX_RENDERBUFFER_SIZE:
    _lsglStorei(params, 1, kMaxTextureWidthHeight);
    break;

  case GL_MAX_TEXTURE_IMAGE_UNITS:
    _lsglStorei(params, 1, kMaxTextureImageUnits);
    break;

  case GL_MAX_TEXTURE_SIZE:
    _lsglStorei(params, 1, kMaxTextureWidthHeight);
    break;

  case GL_MAX_VARYING_VECTORS:
    _lsglStorei(params, 1, kMaxVaryingVectors);
    break;

  case GL_MAX_VERTEX_ATTRIBS:
    _lsglStorei(params, 1, kMaxVertexAttribs);
    break;

  case GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS:
    _lsglStorei(params, 1, kMaxVertexTextureUnits);
    break;

  case GL_MAX_VERTEX_UNIFORM_VECTORS:
    _lsglStorei(params, 1, kMaxVertexUniformVectors);
    break;

  case GL_MAX_VIEWPORT_DIMS:
    _lsglStorei(params, 2, kMaxTextureWidthHeight, kMaxTextureWidthHeight);
    break;

  case GL_NUM_COMPRESSED_TEXTURE_FORMATS:
    _lsglStorei(params, 1, 0);
    break;

  case GL_NUM_SHADER_BINARY_FORMATS:
    _lsglStorei(params, 1, 0);
    break;

  case GL_PACK_ALIGNMENT:
    _lsglStorei(params, 1, state_.packAlignment);
    break;

  case GL_POLYGON_OFFSET_FACTOR:
    _lsglStoref(params, 1, state_.polygonOffsetFactor);
    break;

  case GL_POLYGON_OFFSET_FILL:
    _lsglStoreb(params, 1, state_.polygonOffsetFill);
    break;

  case GL_POLYGON_OFFSET_UNITS:
    _lsglStoref(params, 1, state_.polygonOffsetUnits);
    break;

  case GL_RENDERBUFFER_BINDING:
    _lsglStorei(params, 1, state_.renderBuffer);
    break;

  case GL_SAMPLE_ALPHA_TO_COVERAGE:
    _lsglStoreb(params, 1, state_.sampleAlphaToCoverage);
    break;

  case GL_SAMPLE_BUFFERS:
    _lsglStorei(params, 1, kSampleBuffersPerRenderbuffer);
    break;

  case GL_SAMPLE_COVERAGE:
    _lsglStoreb(params, 1, state_.sampleCoverage);
    break;

  case GL_SAMPLE_COVERAGE_INVERT:
    _lsglStoreb(params, 1, state_.sampleCoverageInvert);
    break;

  case GL_SAMPLE_COVERAGE_VALUE:
    _lsglStoref(params, 1, state_.sampleCoverageValue);
    break;

  case GL_SAMPLES:
    _lsglStorei(params, 1, kMaxSamplesPerPixel);
    break;

  case GL_SCISSOR_BOX:
    _lsglStorei(params, 4, state_.scissorX, state_.scissorY,
                state_.scissorWidth, state_.scissorHeight);
    break;

  case GL_SCISSOR_TEST:
    _lsglStoreb(params, 1, state_.scissorTest);
    break;

  case GL_SHADER_BINARY_FORMATS:
    // no binary shaders are currently supported (we return zero for
    // GL_NUM_SHADER_BINARY_FORMATS), so do nothing here
    break;

  case GL_SHADER_COMPILER:
    _lsglStoreb(params, 1, true);
    break;

  case GL_STENCIL_BACK_FAIL:
    _lsglStorei(params, 1, state_.stencilBackFail);
    break;

  case GL_STENCIL_BACK_FUNC:
    _lsglStorei(params, 1, state_.stencilBackFunc);
    break;

  case GL_STENCIL_BACK_PASS_DEPTH_FAIL:
    _lsglStorei(params, 1, state_.stencilBackPassDepthFail);
    break;

  case GL_STENCIL_BACK_PASS_DEPTH_PASS:
    _lsglStorei(params, 1, state_.stencilBackPassDepthPass);
    break;

  case GL_STENCIL_BACK_REF:
    _lsglStorei(params, 1, state_.stencilBackRef);
    break;

  case GL_STENCIL_BACK_VALUE_MASK:
    _lsglStorei(params, 1, state_.stencilBackMask);
    break;

  case GL_STENCIL_BACK_WRITEMASK:
    _lsglStorei(params, 1, state_.stencilBackWritemask);
    break;

  case GL_STENCIL_BITS: {
    Framebuffer *fb = GetCurrentFrameBuffer();
    if (fb == NULL) {
      return SetGLError(GL_INVALID_FRAMEBUFFER_OPERATION);
    }

    Renderbuffer *rb = fb->GetStencilBuffer();
    if (rb == NULL) {
      _lsglStorei(params, 1, 0);
    } else {
      _lsglStorei(params, 1, rb->GetBitsPerPixel());
    }
  } break;

  case GL_STENCIL_CLEAR_VALUE:
    _lsglStorei(params, 1, state_.stencilClearValue_);
    break;

  case GL_STENCIL_FAIL:
    _lsglStorei(params, 1, state_.stencilFail);
    break;

  case GL_STENCIL_FUNC:
    _lsglStorei(params, 1, state_.stencilFunc);
    break;

  case GL_STENCIL_PASS_DEPTH_FAIL:
    _lsglStorei(params, 1, state_.stencilPassDepthFail);
    break;

  case GL_STENCIL_PASS_DEPTH_PASS:
    _lsglStorei(params, 1, state_.stencilPassDepthPass);
    break;

  case GL_STENCIL_REF:
    _lsglStorei(params, 1, state_.stencilRef);
    break;

  case GL_STENCIL_TEST:
    _lsglStoreb(params, 1, state_.stencilTest);
    break;

  case GL_STENCIL_VALUE_MASK:
    _lsglStorei(params, 1, state_.stencilMask);
    break;

  case GL_STENCIL_WRITEMASK:
    _lsglStorei(params, 1, state_.stencilWritemask);
    break;

  case GL_SUBPIXEL_BITS:
    _lsglStorei(params, 1, 4);
    break;

  case GL_TEXTURE_BINDING_1D:
    _lsglStorei(params, 1,
                state_.texture1D[state_.activeSampler - GL_TEXTURE0]);
    break;

  case GL_TEXTURE_BINDING_2D:
    _lsglStorei(params, 1,
                state_.texture2D[state_.activeSampler - GL_TEXTURE0]);
    break;

  case GL_TEXTURE_BINDING_3D:
    _lsglStorei(params, 1,
                state_.texture3D[state_.activeSampler - GL_TEXTURE0]);
    break;

  case GL_TEXTURE_BINDING_CUBE_MAP:
    _lsglStorei(params, 1,
                state_.textureCubeMap[state_.activeSampler - GL_TEXTURE0]);
    break;

  case GL_UNPACK_ALIGNMENT:
    _lsglStorei(params, 1, state_.unpackAlignment);
    break;

  case GL_VIEWPORT:
    _lsglStorei(params, 4, state_.viewportX, state_.viewportY,
                state_.viewportWidth, state_.viewportHeight);
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glGetBooleanv(GLenum pname, GLboolean *params) {
  GetValue(pname, params);
}

void Context::glGetIntegerv(GLenum pname, GLint *params) {
  GetValue(pname, params);
}

void Context::glGetFloatv(GLenum pname, GLfloat *params) {
  GetValue(pname, params);
}

const GLubyte *Context::glGetString(GLenum name) {
  const char *result = NULL;

  switch (name) {
  case GL_VENDOR:
    result = "GLES LSGL Vendor";
    break;

  case GL_RENDERER:
    result = "GLES LSGL Renderer 1.0";
    break;

  case GL_VERSION:
    result = "2.0";
    break;

  case GL_SHADING_LANGUAGE_VERSION:
    result = "1.0";
    break;

  case GL_EXTENSIONS:
    result = "";
    break;

  default:
    SetGLError(GL_INVALID_ENUM);
    break;
  }

  return reinterpret_cast<const GLubyte *>(result);
}

void Context::glCullFace(GLenum mode) {
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glFrontFace(GLenum mode) {
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glPolygonOffset(GLfloat factor, GLfloat units) {
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}
