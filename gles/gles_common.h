/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2016 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_GLES_COMMON_H__
#define __LSGL_GLES_COMMON_H__

#include <set>
#include <vector>
#include <string>
#include <algorithm>

#include "../include/GLES2/gl2.h"
#include "../include/GLES2/gl2ext.h"
#include "../render/render_timerutil.h"
#include "../render/render_camera.h"
#include "../render/render_texture.h"
#include "../render/render_accel_volume.h"
#include "../render/stack_container.h"
#include "../glsl/glsl_runtime.h"

#ifdef LSGL_DEBUG_TRACE
#define TRACE_EVENT(message, ...)                                              \
  fprintf(stdout, "[LSGL] %s(%s:%d) " message "\n", __FUNCTION__, __FILE__,    \
          __LINE__, ##__VA_ARGS__)
#else
#define TRACE_EVENT(message, ...) (void(0))
#endif

#define LSGL_ASSERT(cond, message, ...) { \
  if (!(cond)) { \
    fprintf(stderr, "[LSGL] " #cond " : %s(%s:%d) " message "\n", __FUNCTION__, __FILE__,    \
     __LINE__, ##__VA_ARGS__); \
    abort(); \
  } \
}

#ifdef WIN32
#ifdef LSGL_EXPORT
#define LSGLES_EXPORT __declspec(dllexport)
#else
#define LSGLES_EXPORT __declspec(dllimport)
#endif
#else
#define LSGLES_EXPORT
#endif

extern "C" {

// Callback function from LSGL to user program
typedef bool (*LSGLProgressCallback)(int progress, int y, int height,
                                     void *userdata);
}

namespace lsgl {

using namespace lsgl::render;

// Value must be less than `glsl/glsl_runtine.h`
const int kMaxVertexAttribs = 8;
const int kMaxVertexUniformVectors = 128;
const int kMaxFragmentUniformVectors = 32;
const int kMaxVaryingVectors = 8;
const int kMaxTextureImageUnits = 8;
const int kMaxColorAttachments = 1;
const int kMaxTextureWidthHeight = 4096;
const int kMaxVertexTextureUnits = 0;
const int kSampleBuffersPerRenderbuffer = 1;
const int kMaxSamplesPerPixel = 8;
const int kMaxDrawStack = 16;
const float kMinPointSize = 1.0f;
const float kMaxPointSize = 1.0f;
const float kMinLineWidth = 1.0f;
const float kMaxLineWidth = 1.0f;

const int kVtxAttrPosition = 0;
const int kVtxAttrNormal = 1;
const int kVtxAttrTexCoord = 2;
const int kVtxAttrMaterial = 3;

class Context;

/// GLES color
struct Color {
  float red;
  float green;
  float blue;
  float alpha;
};

/// GLES vertex attribute
struct VertexAttribute {
  bool enabled;

  // For glVertexAttribPointer
  GLint size;
  GLenum type;
  GLboolean normalized;
  GLsizei stride;
  const GLvoid *ptr;

  // For vertex buffer
  int handle;

  // For glVertexAttrib
  bool uniform;
  union {
    GLfloat fval1;
    GLfloat fval2[2];
    GLfloat fval3[3];
    GLfloat fval4[4];
  } u;
};

/// GLES state
struct State {
  Color colorClearValue_;
  GLclampf depthClearValue_;
  GLint stencilClearValue_;

  bool cullFace;
  GLenum cullMode;
  GLenum frontFace;
  bool depthTest;
  GLenum depthFunc;
  GLclampf nearDepth;
  GLclampf farDepth;
  bool blend;
  GLenum sourceBlendRGB;
  GLenum destBlendRGB;
  GLenum sourceBlendAlpha;
  GLenum destBlendAlpha;
  GLenum blendEquationRGB;
  GLenum blendEquationAlpha;
  Color blendColor;
  bool stencilTest;
  GLenum stencilFunc;
  GLint stencilRef;
  GLuint stencilMask;
  GLenum stencilFail;
  GLenum stencilPassDepthFail;
  GLenum stencilPassDepthPass;
  GLuint stencilWritemask;
  GLenum stencilBackFunc;
  GLint stencilBackRef;
  GLuint stencilBackMask;
  GLenum stencilBackFail;
  GLenum stencilBackPassDepthFail;
  GLenum stencilBackPassDepthPass;
  GLuint stencilBackWritemask;
  bool polygonOffsetFill;
  GLfloat polygonOffsetFactor;
  GLfloat polygonOffsetUnits;
  bool sampleAlphaToCoverage;
  bool sampleCoverage;
  GLclampf sampleCoverageValue;
  bool sampleCoverageInvert;
  bool scissorTest;
  bool dither;

  GLfloat lineWidth;

  GLenum generateMipmapHint;
  GLenum fragmentShaderDerivativeHint;

  GLint viewportX;
  GLint viewportY;
  GLsizei viewportWidth;
  GLsizei viewportHeight;
  float zNear;
  float zFar;

  GLint scissorX;
  GLint scissorY;
  GLsizei scissorWidth;
  GLsizei scissorHeight;

  bool colorMaskRed;
  bool colorMaskGreen;
  bool colorMaskBlue;
  bool colorMaskAlpha;
  bool depthMask;

  GLenum activeSampler; // Active texture unit selector - GL_TEXTURE0
  GLuint arrayBuffer;
  GLuint elementArrayBuffer;
  GLuint frameBuffer;
  GLuint renderBuffer;

  GLuint readFramebuffer;
  GLuint drawFramebuffer;
  GLuint currentProgram;
  GLuint texture1D[kMaxTextureImageUnits];
  GLuint texture2D[kMaxTextureImageUnits];
  GLuint texture3D[kMaxTextureImageUnits];
  GLuint textureCubeMap[kMaxTextureImageUnits];

  std::vector<VertexAttribute> vertexAttributes[kMaxDrawStack];
  int currentDrawStackIndex;

  GLint unpackAlignment;
  GLint packAlignment;
  bool packReverseRowOrder;

  GLenum lastError;

  // LSGL
  GLint pixelStep;
};

// Used in sparse LoD texture.
struct Region {
  GLint level;          ///< LoD level
  GLint offset[3];      ///< Offset in 3D coordinate
  GLsizei extent[3];    ///< Extent in 3D coordinate
  GLsizei size[3];      ///< Cell(Voxel) size
  unsigned char *data;  ///< Voxel data.
  GLboolean commit;     ///< Page commit.
};

/// Base class for GLES vertex buffer.
class Buffer {
public:
  Buffer();
  virtual ~Buffer();

  void Data(size_t size, const void *data, GLenum usage);
  void SubData(unsigned int offset, size_t size, const void *data);

  void Retain(size_t size, const void *data, GLenum usage);

  inline const GLubyte *GetData() const {
    if (retained_) {
      return ptr_;
    } else {
      return &data_[0];
    }
  }
  inline GLsizeiptr GetSize() const { return size_; }
  inline GLenum GetUsage() const { return usage_; }

  // inline bool IsStatic() const { return (GetUsage() == GL_STATIC_DRAW) ||
  //(GetUsage() == GL_STATIC_READ) || (GetUsage() == GL_STATIC_COPY); }
  inline bool IsStatic() const { return (GetUsage() == GL_STATIC_DRAW); }

  bool IsRetained() const { return retained_; }

private:
  std::vector<GLubyte> data_;
  GLsizeiptr size_;
  GLenum usage_;
  const unsigned char *ptr_;
  bool retained_;
};

namespace local {

// Input range must be [0, 1)
static inline float remap(float x, const float *table, int n) {
  int idx = x * n;
  idx = (std::max)((std::min)(n - 1, idx), 0);

  return table[idx];
}
}

/// Base class for GLES texture data provider.
class Texture {
public:
  Texture();
  virtual ~Texture();

  // Make a internal copy of texture data.
  void Data2D(const void *data, GLuint width, GLuint height, int compos,
              GLenum type);
  void Data3D(const void *data, GLuint width, GLuint height, GLuint depth,
              int compos, GLenum type);

  // Just retain a pointer to the texture data(no internal copy happens).
  void Retain2D(const void *data, GLuint width, GLuint height, int compos,
                GLenum type);
  void Retain3D(const void *data, GLuint width, GLuint height, GLuint depth,
                int compos, GLenum type);

  void SubImage3D(GLint xoffset, GLint yoffset, GLint zoffset, GLsizei width,
                  GLsizei height, GLsizei depth, int compos, GLenum type,
                  const GLvoid *data);

  // Just retain a pointer to the texture data(no internal copy happens).
  void SubImage3DRetain(GLint level, GLint xoffset, GLint yoffset, GLint zoffset,
                        GLsizei width, GLsizei height, GLsizei depth,
                        GLsizei cellWidth, GLsizei cellHeight, GLsizei cellDepth,
                        int compos, GLenum type, const GLvoid *data);

  void SetMinFilter(bool filter) { minFiltering_ = filter; }

  void SetMagFilter(bool filter) { magFiltering_ = filter; }

  void SetWrapS(GLenum mode) {
    if (mode == GL_REPEAT) {
      if (texture2D_) texture2D_->setWrapS(render::Texture2D::WRAP_REPEAT);
      if (texture3D_) texture3D_->wrapS = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
    } else if (mode == GL_CLAMP_TO_EDGE) {
      if (texture2D_) texture2D_->setWrapS(render::Texture2D::WRAP_CLAMP_TO_EDGE);
      if (texture3D_) texture3D_->wrapS = LSGL_RENDER_TEXTURE3D_WRAP_CLAMP_TO_EDGE;
    } else {
      // @todo
    }
  }

  void SetWrapT(GLenum mode) {
    if (mode == GL_REPEAT) {
      if (texture2D_) texture2D_->setWrapT(render::Texture2D::WRAP_REPEAT);
      if (texture3D_) texture3D_->wrapT = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
    } else if (mode == GL_CLAMP_TO_EDGE) {
      if (texture2D_) texture2D_->setWrapT(render::Texture2D::WRAP_CLAMP_TO_EDGE);
      if (texture3D_) texture3D_->wrapT = LSGL_RENDER_TEXTURE3D_WRAP_CLAMP_TO_EDGE;
    } else {
      // @todo
    }
  }

  void SetWrapR(GLenum mode) {
    if (mode == GL_REPEAT) {
      if (texture3D_) texture3D_->wrapR = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
    } else if (mode == GL_CLAMP_TO_EDGE) {
      if (texture3D_) texture3D_->wrapR = LSGL_RENDER_TEXTURE3D_WRAP_CLAMP_TO_EDGE;
    } else {
      // @todo
    }
  }

  inline void Fetch(float *rgba, float u, float v) const {
    if (!texture2D_)
      return;
    texture2D_->fetch(rgba, u, v, minFiltering_, magFiltering_);
  }

  inline void Fetch(float *rgba, float u, float v, float r) const {

    if (isSparse_) {

      rgba[0] = rgba[1] = rgba[2] = 0.0f;
      rgba[3] = 0.0f; // Alpha zero = miss

      if (!sparseVolumeAccel_) {
        return;
      }

      if (sparseVolumeAccel_ && sparseVolume_) { // this shuld be always true.
        double value[4];
        double position[3];

        // Asume u, v, r is in [0, 1) range.
        position[0] = u * sparseVolume_->globalDim[0];
        position[1] = v * sparseVolume_->globalDim[1];
        position[2] = r * sparseVolume_->globalDim[2];
        sparseVolumeAccel_->Sample(value, position);
        rgba[0] = value[0];
        rgba[1] = value[1];
        rgba[2] = value[2];
        rgba[3] = value[3]; // @todo { set alpha to 1.0 might be better. }
      }

    } else {

      if (!texture3D_)
        return;

      float uu = u;
      float vv = v;
      float rr = r;

      // Apply coordinate remapping.
      if (doRemap_[0] && remapTable_[0].size() > 0) {
        uu =
            local::remap(uu, &remapTable_[0].at(0), (int)remapTable_[0].size());
      }
      if (doRemap_[1] && remapTable_[1].size() > 0) {
        vv =
            local::remap(vv, &remapTable_[1].at(0), (int)remapTable_[1].size());
      }
      if (doRemap_[2] && remapTable_[2].size() > 0) {
        rr =
            local::remap(rr, &remapTable_[2].at(0), (int)remapTable_[2].size());
      }

      if (texture3D_->data_type == LSGL_RENDER_TEXTURE3D_FORMAT_DOUBLE) {
        FilterTexture3DDouble(rgba, texture3D_, uu, vv, rr, texture3D_->wrapS, texture3D_->wrapT, texture3D_->wrapR);
      } else if (texture3D_->data_type == LSGL_RENDER_TEXTURE3D_FORMAT_FLOAT) {
        FilterTexture3DFloat(rgba, texture3D_, uu, vv, rr, texture3D_->wrapS, texture3D_->wrapT, texture3D_->wrapR);
      } else if (texture3D_->data_type == LSGL_RENDER_TEXTURE3D_FORMAT_BYTE) {
        FilterTexture3DByte(rgba, texture3D_, uu, vv, rr, texture3D_->wrapS, texture3D_->wrapT, texture3D_->wrapR);
      } else {
        LSGL_ASSERT(0, "Unknown 3D texture format");
      }
    }
  }

  Texture3D *GetTexture3D() const { return texture3D_; }

  bool RemoveRemapTable(GLenum coord) {
    if (coord == GL_COORDINATE_X) {
      doRemap_[0] = false;
      remapTable_[0].clear();
    } else if (coord == GL_COORDINATE_Y) {
      doRemap_[1] = false;
      remapTable_[1].clear();
    } else if (coord == GL_COORDINATE_Z) {
      doRemap_[2] = false;
      remapTable_[2].clear();
    } else {
      return false;
    }

    return true;
  }

  bool SetRemapTable(GLenum coord, GLsizei size, const GLfloat *coords);

  bool IsSparse() { return isSparse_; }

  // Builds internal state of sparse texture.
  void BuildSparseTexture();

  void SetSparsity(bool isSparse) { isSparse_ = isSparse; }

  std::vector<Region> &GetRegionList() { return regionList_; }

  int GetCompos() { return compos_; }

  GLenum GetType() { return type_; }

private:
  void Free();

  int compos_;
  GLenum type_;

  render::Texture2D *texture2D_;
  Texture3D *texture3D_;
  std::vector<GLubyte> data_;
  GLsizeiptr size_;
  bool minFiltering_;
  bool magFiltering_;
  bool retained_;
  bool doRemap_[3];                    // x, y and z
  std::vector<GLfloat> remapTable_[3]; // x, y and z

  // For spase texture.
  // @todo { Refactor. }
  bool isSparse_;
  std::vector<Region> regionList_;
  SparseVolumeAccel *sparseVolumeAccel_;
  SparseVolume *sparseVolume_;
};

/// GLES renderbuffer.
class LSGLES_EXPORT Renderbuffer {
public:
  Renderbuffer();
  virtual ~Renderbuffer();

  /// Allocate buffer for renderbuffer.
  bool Allocate(GLuint width, GLuint height, GLuint bytesPerPixel,
                GLenum format);

  inline GLuint GetWidth() const { return width_; }
  inline GLuint GetHeight() const { return height_; }
  inline GLuint GetBytesPerPixel() const { return bytesPerPixel_; }
  inline GLuint GetBitsPerPixel() const { return GetBytesPerPixel() * 8; }
  inline GLuint GetSize() const {
    return GetWidth() * GetHeight() * GetBytesPerPixel();
  }
  inline GLenum GetFormat() const { return format_; }
  inline char *GetBuffer() { return buffer_; }

  /// True if a renderbuffer has same size for specified renderbuffer.
  inline bool IsSameSize(const Renderbuffer *rb) const {
    return (width_ == rb->width_) && (height_ == rb->height_);
  }

private:
  void Free();

  GLuint width_;
  GLuint height_;
  GLuint bytesPerPixel_;
  GLenum format_;
  char *buffer_;
};

/// GLES framebuffer
class LSGLES_EXPORT Framebuffer {
public:
  Framebuffer();
  virtual ~Framebuffer();

  void AttachColorBuffer(GLuint index, Renderbuffer *rb);
  void AttachDepthBuffer(Renderbuffer *rb);
  void AttachStencilBuffer(Renderbuffer *rb);

  inline Renderbuffer *GetColorBuffer(GLuint index) {
    assert(index < kMaxColorAttachments);
    return colorBuffers_[index];
  }
  inline Renderbuffer *GetDepthBuffer() { return depthBuffer_; }
  inline Renderbuffer *GetStencilBuffer() { return stencilBuffer_; }

  inline bool HasAttachment() const { return attachments_; }
  inline bool HasValidSize() const { return validSize_; }
  inline GLuint GetWidth() const { return width_; }
  inline GLuint GetHeight() const { return height_; }

private:
  void UpdateAttachmentStatus();

  Renderbuffer *colorBuffers_[kMaxColorAttachments];
  Renderbuffer *depthBuffer_;
  Renderbuffer *stencilBuffer_;
  GLuint width_;
  GLuint height_;
  bool attachments_;
  bool validSize_;
};

/// Represents GLES shader uniform variable
struct Uniform {
  Uniform(GLenum _type, const char *_name, int _arraySize)
      : name(std::string(_name)), type(_type), arraySize(_arraySize) {}

  bool isArray() { return (arraySize > 1); }

  GLenum type;
  std::string name;
  int arraySize;

  // Uniform variable is implemented as opaque data
  std::vector<unsigned char> data;
};

struct UniformLocation {
  UniformLocation(const char *_name, unsigned int _element,
                  unsigned int _index) {
    name = std::string(_name);
    element = _element;
    index = _index;
  }

  std::string name;
  unsigned int element;
  unsigned int index;
};

/// Represents GLES shader varying variable
class Varying {
public:
  Varying(GLenum _type, int _count, const char *_name, int _arraySize)
      : name(std::string(_name)), type(_type), count(_count),
        arraySize(_arraySize) {}
  ~Varying() {}

  Varying(const Varying &rhs)
      : type(rhs.type), name(rhs.name), count(rhs.count),
        arraySize(rhs.arraySize), data(rhs.data) {}

  bool IsArray() const { return (arraySize > 1) ? true : false; }

  GLenum type;
  std::string name;
  int count;
  int arraySize;

  // Varuing variable is implemented as opaque data
  std::vector<unsigned char> data;
};

struct VaryingLocation {
  VaryingLocation(const char *_name, unsigned int _element,
                  unsigned int _index) {
    name = std::string(_name);
    element = _element;
    index = _index;
  }

  std::string name;
  unsigned int element;
  unsigned int index;
};

///
/// GLES Attrib <-> GLSL Varying connection table
///
struct VaryingConnection {
  int srcIndex;             // the index to vertex attrib.
  int dstIndex;             // the slot index to varying variable.
  const unsigned char *bufPtr; // the pointer to vertex attribute data.
  size_t bufSize; // Size of vertex attribute data(data size passed to glBufferData).
  int size;
  int stride;
};

/// ShadingState holds runtime shader information
struct ShadingState {
  StackVector<VaryingConnection, kMaxVaryingVectors> varyingConnections;
  StackVector<Varying, kMaxVaryingVectors> varyings; // buffer for varying variables

  void Copy(ShadingState &src) {
    varyingConnections->clear();
    // printf("sz = %d\n", (int)src.varyingConnections.size());
    for (size_t i = 0; i < src.varyingConnections->size(); i++) {
      assert(src.varyingConnections[i].dstIndex < 100);
      varyingConnections->push_back(src.varyingConnections[i]);
    }

    varyings->clear();
    // printf("varuings.sz = %d\n", (int)src.varyings.size());
    for (size_t i = 0; i < src.varyings->size(); i++) {
      varyings->push_back(src.varyings[i]);
    }
  }
};

/// Store intersection state for the shading.
struct IntersectionState {
  real3 position;
  real3 normal;
  real3 geometricNormal;
  real3 tangent;
  real3 binormal;
  real3 raydir;
  int raydepth;
  float px;
  float py;
  int doubleSided;
  float rayattrib;
  GLenum primitiveType;  // GL_TRIANGLES, GL_TETRAHEDRONS_EXT, etc
  const unsigned char *prev_node;
  unsigned int prev_prim_id;
  double prev_hit_t;
  real3 prev_hit_normal;
  float u;
  float v;
  unsigned int f0;
  unsigned f1;
  unsigned f2;
  unsigned f3;
  unsigned f4;
  unsigned f5;
  unsigned f6;
  unsigned f7;
  float d0;
  float d1;
  float d2;
  float d3;
  float d4;
  float d5;
  float d6;
  float d7;
};

/// Base class of vertex and fragment shader.
class LSGLES_EXPORT Shader {
  //
  // DLL shader method
  //
  typedef struct {
    void *shaderEvalFunc;
    void *shaderInfoFunc;
    void *shaderInitFunc;
  } Method;

public:
  Shader();
  virtual ~Shader();

  virtual GLenum GetType() const = 0;

  inline bool IsCompiled() const { return compiled_; }
  inline const std::string &GetSource() const { return source_; }
#ifdef ENABLE_LLVM
  inline ShaderRuntime *GetShaderRT() const { return &shaderRT_; }
#endif

  /// Loads precompiled dll shader.
  bool LoadShaderBinary(std::string &filename);

  void Source(GLsizei count, const GLchar **string, const GLint *length);
  bool Compile();

  /// Get uniform variable information from the loaded shader and
  /// set result to uniforms and uniformLocations.
  /// Valid after LoadShaderBinary()
  void BuildUniformInfo(std::vector<Uniform> &uniforms,
                        std::vector<UniformLocation> &uniformLocations);

  /// Get varying variable information from the loaded shader and
  /// set result to varyings and varyingLocations.
  /// Valid after LoadShaderBinary()
  void BuildVaryingInfo(std::vector<Varying> &varyings,
                        std::vector<VaryingLocation> &varyingLocations);

  bool Release();

  void SetAttached() { isAttached_ = true; }

  bool IsAttached() const { return isAttached_; }

protected:
  virtual bool DoCompile() = 0;

  std::string source_;

  void *handle_;         /// Shader DLL handle
  std::string filename_; /// Filename(for binary shader)
  Method method_;

private:
  bool compiled_;
  bool isBinary_;
  bool isAttached_;
};

/// Vertex shader
class VertexShader : public Shader {
public:
  VertexShader();
  virtual ~VertexShader();

  inline virtual GLenum GetType() const { return GL_VERTEX_SHADER; }

private:
  virtual bool DoCompile();
};

/// Fragment shader
class FragmentShader : public Shader {
public:
  typedef struct {
    float frame[3][3]; // (eye, lookat, up) in world coord.
    float fov;
  } CameraInfo;

  FragmentShader();
  virtual ~FragmentShader();

  inline virtual GLenum GetType() const { return GL_FRAGMENT_SHADER; }

  /// Prepare shader evaluation. Setup varying/uniform variable, etc.
  bool PrepareEval(FragmentState &fragmentState, // [out]
                   ShadingState &shadingState,   // [out]
                   const std::vector<Uniform> &uniforms,
                   const std::vector<Varying> &varyings,
                   const std::vector<VaryingLocation> &varyingLocations,
                   const std::vector<VertexAttribute> &vertexAttributes,
                   const Context &ctx) const;

  /// Evaluate fragment shader. Returns true if success, false if the fragment
  /// was discarded.
  bool Eval(GLfloat fragColor[4], FragmentState &fragmentState,
            ShadingState &shadingState,
            const std::vector<VertexAttribute> &vertexAttributes,
            const GLfloat fragCoord[4], const IntersectionState &isectState,
            const CameraInfo &cameraInfo, int threadID) const;

private:
  virtual bool DoCompile();

  // FragmentState state_;
};

/// GLES Program object
class Program {
public:
  Program();
  ~Program();

  Program(const Program &prg) { // Copy constructor.

    linked_ = prg.linked_;
    vertexShaders_ = prg.vertexShaders_;
    fragmentShaders_ = prg.fragmentShaders_;

    samplerPSs_ = prg.samplerPSs_;
    samplerVSs_ = prg.samplerVSs_;
    uniforms_ = prg.uniforms_;
    uniformLocations_ = prg.uniformLocations_;
    varyings_ = prg.varyings_;
    varyingLocations_ = prg.varyingLocations_;
  }

  inline bool IsLinked() const { return linked_; }
  inline GLuint GetVertexShaderCount() const {
    return (GLuint)vertexShaders_.size();
  }
  // inline GLuint GetGeometryShaderCount() const { return
  // geometryShaders_.size(); }
  inline GLuint GetFragmentShaderCount() const {
    return (GLuint)fragmentShaders_.size();
  }
  // inline GLuint GetAttachedCount() const { return GetVertexShaderCount() +
  // GetGeometryShaderCount() + GetFragmentShaderCount(); }
  inline GLuint GetAttachedCount() const {
    return GetVertexShaderCount() + GetFragmentShaderCount();
  }

  inline const VertexShader *GetVertexShader(GLuint idx) const {
    return vertexShaders_[idx];
  }
  // inline const GeometryShader *GetGeometryShader(GLuint idx) const { return
  // geometryShaders_[idx]; }
  inline const FragmentShader *GetFragmentShader(GLuint idx) const {
    return fragmentShaders_[idx];
  }

  inline FragmentShader *GetFragmentShader(GLuint idx) {
    return fragmentShaders_[idx];
  }

  bool AttachShader(Shader *shader);
  bool DetachShader(Shader *shader);

  /// Links vertex and fragment shaders.
  bool Link();

  bool SetUniform1iv(GLint location, GLsizei count, const GLint *v);

  bool SetUniform1fv(GLint location, GLsizei count, const GLfloat *v);
  bool SetUniform2fv(GLint location, GLsizei count, const GLfloat *v);
  bool SetUniform3fv(GLint location, GLsizei count, const GLfloat *v);
  bool SetUniform4fv(GLint location, GLsizei count, const GLfloat *v);

  bool SetUniformMatrix2fv(GLint location, GLsizei count, GLboolean transpose,
                           const GLfloat *v);
  bool SetUniformMatrix3fv(GLint location, GLsizei count, GLboolean transpose,
                           const GLfloat *v);
  bool SetUniformMatrix4fv(GLint location, GLsizei count, GLboolean transpose,
                           const GLfloat *v);

  /// Look ups uniform slot from variable name
  GLint GetUniformLocation(const std::string &name) const;

  /// Get uniform variable at specified location.
  Uniform const *GetUniform(GLint location) const;

  /// Look ups varying slot from variable name
  GLint GetVaryingLocation(const std::string &name) const;

  /// Get varying variable at specified location.
  Varying const *GetVarying(GLint location) const;

  /// Prepare shader evaluation. Setup uniform variable for example.
  bool PrepareEval(FragmentState &fragmentState, ShadingState &shadingState,
                   const std::vector<VertexAttribute> &vertexAttributes,
                   const Context &ctx) const;

private:
  bool IsAttached(Shader *shader);

  bool linked_;
  std::vector<VertexShader *> vertexShaders_;
  std::vector<FragmentShader *> fragmentShaders_;

  struct Sampler {
    Sampler();

    bool active;
    GLint textureUnit;
  };

  // Program state
  std::vector<Sampler> samplerPSs_;
  std::vector<Sampler> samplerVSs_;
  std::vector<Uniform> uniforms_;
  std::vector<UniformLocation> uniformLocations_;
  std::vector<Varying> varyings_;
  std::vector<VaryingLocation> varyingLocations_;
};

/// Utility class for data buffer
template <typename T> class DataBuffer {
public:
  DataBuffer() : count_(0) {}

  T *Add(unsigned int count);

  inline void Clear() { count_ = 0; }
  inline GLuint GetCount() const { return count_; }
  inline T *GetBase() { return (count_ == 0) ? NULL : &buffer_[0]; }
  inline GLuint GetByteSize() const {
    return (GLuint)(buffer_.size() * sizeof(T));
  }

private:
  std::vector<T> buffer_;
  unsigned int count_;
};

template <typename T> T *DataBuffer<T>::Add(unsigned int count) {
  // increase buffer size first as necessary
  buffer_.resize(count_ + count);

  // get pointer to buffer location and increment internal count
  T *ptr = &buffer_[count_];
  count_ += count;
  return ptr;
}

// random number utils
void initialize_random();
float randomreal(int thread_id);

} // namespace lsgl

#endif // __LSGL_GLES_COMMON_H__
