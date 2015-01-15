/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cstring>

#include "gles_context.h"

// must be included inside a namespace to avoid ambiguities
namespace lsgl {

static Context sContext;

Context::Context() : dirtyFrame_(false) {
  // initialize state to default values
  state_.arrayBuffer = 0;
  state_.activeSampler = GL_TEXTURE0;
  state_.blend = false;
  state_.blendColor.red = 0.0f;
  state_.blendColor.green = 0.0f;
  state_.blendColor.blue = 0.0f;
  state_.blendColor.alpha = 0.0f;
  state_.destBlendAlpha = GL_ZERO;
  state_.destBlendRGB = GL_ZERO;
  state_.blendEquationAlpha = GL_FUNC_ADD;
  state_.blendEquationRGB = GL_FUNC_ADD;
  state_.sourceBlendAlpha = GL_ONE;
  state_.sourceBlendRGB = GL_ONE;
  state_.colorClearValue_.red = 0.0f;
  state_.colorClearValue_.green = 0.0f;
  state_.colorClearValue_.blue = 0.0f;
  state_.colorClearValue_.alpha = 0.0f;
  state_.colorMaskRed = true;
  state_.colorMaskGreen = true;
  state_.colorMaskBlue = true;
  state_.colorMaskAlpha = true;
  state_.cullFace = false;
  state_.cullMode = GL_BACK;
  state_.currentProgram = 0;
  state_.depthClearValue_ = 1.0f;
  state_.depthFunc = GL_LESS;
  state_.nearDepth = 0.0f;
  state_.farDepth = 1.0f;
  state_.depthTest = false;
  state_.depthMask = true;
  state_.dither = true;
  state_.elementArrayBuffer = 0;
  state_.frameBuffer = 0;
  state_.frontFace = GL_CCW;
  state_.generateMipmapHint = GL_DONT_CARE;
  state_.lastError = GL_NO_ERROR;
  state_.lineWidth = 1.0f;
  state_.packAlignment = 4;
  state_.polygonOffsetFactor = 0.0f;
  state_.polygonOffsetFill = false;
  state_.polygonOffsetUnits = 0.0f;
  state_.renderBuffer = 0;
  state_.sampleAlphaToCoverage = false;
  state_.sampleCoverage = false;
  state_.sampleCoverageInvert = false;
  state_.sampleCoverageValue = 1.0f;
  state_.scissorX = 0;
  state_.scissorY = 0;
  state_.scissorWidth = kMaxTextureWidthHeight;
  state_.scissorHeight = kMaxTextureWidthHeight;
  state_.scissorTest = false;
  state_.stencilBackFail = GL_KEEP;
  state_.stencilBackFunc = GL_ALWAYS;
  state_.stencilBackPassDepthFail = GL_KEEP;
  state_.stencilBackPassDepthPass = GL_KEEP;
  state_.stencilBackRef = 0;
  state_.stencilBackMask = 0xFFFFFFFF;
  state_.stencilBackWritemask = 0xFFFFFFFF;
  state_.stencilClearValue_ = 0;
  state_.stencilFail = GL_KEEP;
  state_.stencilFunc = GL_ALWAYS;
  state_.stencilPassDepthFail = GL_KEEP;
  state_.stencilPassDepthPass = GL_KEEP;
  state_.stencilRef = 0;
  state_.stencilTest = false;
  state_.stencilMask = 0xFFFFFFFF;
  state_.stencilWritemask = 0xFFFFFFFF;
  state_.unpackAlignment = 4;
  state_.viewportX = 0;
  state_.viewportY = 0;
  state_.viewportWidth = kMaxTextureWidthHeight;
  state_.viewportHeight = kMaxTextureWidthHeight;

  for (GLuint i = 0; i < kMaxTextureImageUnits; i++) {
    state_.texture1D[i] = state_.texture2D[i] = state_.texture3D[i] =
        state_.textureCubeMap[i] = 0;
  }

  state_.currentDrawStackIndex = 0;
  for (int k = 0; k < kMaxDrawStack; k++) {
    state_.vertexAttributes[k].resize(kMaxVertexAttribs);
    for (GLuint i = 0; i < kMaxVertexAttribs; i++) {
      state_.vertexAttributes[k][i].enabled = false;
      state_.vertexAttributes[k][i].size = 0;
      state_.vertexAttributes[k][i].type = 0;
      state_.vertexAttributes[k][i].normalized = GL_FALSE;
      state_.vertexAttributes[k][i].stride = 0;
      state_.vertexAttributes[k][i].ptr = 0;
    }
  }

  // LSGL
  state_.pixelStep = 1;

  engine_.OnPrepare(this);

  renderGraph_ = new RenderGraph();
}

Context::~Context() { delete renderGraph_; }

Context &Context::GetCurrentContext() { return sContext; }

void Context::glViewport(GLint x, GLint y, GLsizei width, GLsizei height) {
  state_.viewportX = x;
  state_.viewportY = y;
  state_.viewportWidth = width;
  state_.viewportHeight = height;
}

void Context::lsglSetProgressCallback(LSGLProgressCallback func) {
  engine_.SetProgressCallback(func);
}

void Context::lsglSetPixelStep(GLint step) { state_.pixelStep = step; }

GLint Context::lsglGetPixelStep() { return state_.pixelStep; }

// @fixme { Deprecated }
void Context::lsglSetPointSize(GLfloat size) {
  // state_.pointSize = size;
  fprintf(stderr, "[LSGL] lsglSetPointSize is deprecated. use "
                  "\"lsgl_PointSize\" in Uniform insted.");
}

void Context::lsglSetPointSizev(GLsizei num, const GLfloat *size) {
  // state_.pointSizev.resize(num);
  // memcpy(&state_.pointSizev.at(0), size, sizeof(float) * num);
  fprintf(stderr, "[LSGL] lsglSetPointSizev is deprecated. use "
                  "\"lsgl_PointSize\" in VertexAttribute insted.");
}

void Context::lsglSetCamera(const GLfloat *eye, const GLfloat *target, const GLfloat *up,
                            GLfloat fov) {
  engine_.SetCamera(eye, target, up, fov);
}

void Context::lsglSetStereoEnvCamera(const GLfloat *eye, const GLfloat *target, const GLfloat *up,
                            GLfloat zeroParallax, GLfloat eyeSeparation) {
  engine_.SetStereoEnvCamera(eye, target, up, zeroParallax, eyeSeparation);
}

void Context::lsglSetShaderCompiler(const char *path, const char *options) {
// Set path to GLSL_COMPILER environment variable.
#ifdef _WIN32
  _putenv_s("GLSL_COMPILER", path);
#else
  setenv("GLSL_COMPILER", path, 1 /* overwrite */);
#endif
}

static unsigned char fclamp(float f) {
  int i = f * 255;
  if (i < 0)
    i = 0;
  if (i > 255)
    i = 255;

  return (unsigned char)i;
}

void Context::lsglEvalFragmentShader() {
  // Use current porgram
  Program *program = resourceManager_.GetProgram(state_.currentProgram);

  // Render quad.
  const int width = state_.viewportWidth;
  const int height = state_.viewportHeight;

  FragmentState fragmentState;
  ShadingState shadingState;
  program->PrepareEval(fragmentState, shadingState, state_.vertexAttributes[0],
                       (*this));

  // Assume RGBA 8bit.
  char *image = GetCurrentFrameBuffer()->GetColorBuffer(0)->GetBuffer();

  const std::vector<VertexAttribute> &vertexAttributes =
      state_.vertexAttributes[0];

  FragmentShader::CameraInfo cameraInfo; // dummy

  int thread_id = 0;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {

      GLfloat fragcol[4]; // out
      GLfloat fragCoord[4];

      // @fixme { setup z and w correctly. }
      fragCoord[0] = 0.5f + x;
      fragCoord[1] = 0.5f + y;
      fragCoord[2] = 0.0f;
      fragCoord[3] = 0.0f;

      real3 pos, normal, raydir;
      float px = 0.5f + x;
      float py = 0.5f + y;
      float raydepth = 0.0f;
      float rayattrib = 0.0f;
      int doubleSided = 1;
      program->GetFragmentShader(0)->Eval(
          fragcol, fragmentState, shadingState, vertexAttributes, fragCoord,
          pos, normal, raydir, raydepth, px, py, doubleSided, rayattrib, NULL,
          0, 0.0f, 0.0f, 0, 0, 0, cameraInfo, thread_id);

      image[4 * ((height - y - 1) * width + x) + 0] = fclamp(fragcol[0]);
      image[4 * ((height - y - 1) * width + x) + 1] = fclamp(fragcol[1]);
      image[4 * ((height - y - 1) * width + x) + 2] = fclamp(fragcol[2]);
      image[4 * ((height - y - 1) * width + x) + 3] = fclamp(fragcol[3]);
    }
  }
}

void Context::lsglEvalSingleFragmentShader() {

  // @todo { Deprecated. to be removed. }

  // Use current porgram
  Program *program = resourceManager_.GetProgram(state_.currentProgram);

  const std::vector<VertexAttribute> &vertexAttributes =
      state_.vertexAttributes[0];

  FragmentState fragmentState;
  ShadingState shadingState;
  program->PrepareEval(fragmentState, shadingState, vertexAttributes, (*this));

  FragmentShader::CameraInfo cameraInfo;

  GLfloat fragcol[4]; // out
  GLfloat fragCoord[4];
  fragCoord[0] = 0.5f;
  fragCoord[1] = 0.5f;
  fragCoord[2] = 0.0f;
  fragCoord[3] = 0.0f;

  real3 pos, normal, raydir;
  float px = 0.5f;
  float py = 0.5f;
  float raydepth = 0.0f;
  float rayattrib = 0.0f;

  int doubleSided = 1;
  int thread_id = 0;

  program->GetFragmentShader(0)
      ->Eval(fragcol, fragmentState, shadingState, vertexAttributes, fragCoord,
             pos, normal, raydir, raydepth, px, py, doubleSided, rayattrib,
             NULL, 0, 0.0f, 0.0f, 0, 0, 0, cameraInfo, thread_id);
}

Framebuffer *Context::GetCurrentFrameBuffer() {
  // lookup current framebuffer
  Framebuffer *fb;
  if (state_.frameBuffer != 0) {
    fb = resourceManager_.GetFramebuffer(state_.frameBuffer);
    if (fb == NULL) {
      return NULL;
    }
  } else {
    fb = &frameBuffer_;
  }

  // fail if framebuffer is not valid
  if (fb->HasValidSize() == false) {
    return NULL;
  }

  return fb;
}

void Context::SetGLError(GLenum err) {
  state_.lastError = err;
#ifdef LSGL_DEBUG_TRACE
  if (err != GL_NO_ERROR) {
    GLenum errCode = err;

    const char *str = "";

    switch (errCode) {
    case GL_INVALID_ENUM:
      str = "GL_INVALID_ENUM";
      break;
    case GL_INVALID_VALUE:
      str = "GL_INVALID_VALUE";
      break;
    case GL_INVALID_OPERATION:
      str = "GL_INVALID_OPERATION";
      break;
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      str = "GL_INVALID_FRAMEBUFFER_OPERATION";
      break;
    case GL_OUT_OF_MEMORY:
      str = "GL_INVALID_OPERATION";
      break;
    default:
      fprintf(stderr, "%d\n", err);
      str = "??? Unknown GL Error";
      break;
    }

    fprintf(stderr, "%s\n", str);
  }
#endif
}

} // namespace lsgl
