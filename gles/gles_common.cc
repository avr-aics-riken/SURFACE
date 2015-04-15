/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cassert>
#include <cstring>
#include <cerrno>
#include <cstdlib>
#include <algorithm>

#include <fstream>

#include "GLES2/gl2.h"
#include "GLES2/gl2ext.h"

#include "gles_context.h"

#ifdef _WIN32
//#include "MemoryModule.h" // Load DLL from memory.
#else
#include <dlfcn.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

// Interface with shader function
#include "../glsl/glsl_runtime.h"
#include "../render/tinymt64.h"

using namespace lsgl;

static tinymt64_t tinymt64[512]; // Up to 512 threads.

// Shader function pointer signature
typedef void (*ShaderEvalFuncPtr)(Fragment *__fragment, FragmentState *__state);
typedef int (*ShaderInfoFuncPtr)(FragmentConfig *config);

static void GetGLType(GLenum &compoundType, GLenum &baseType, int &count,
                      const char *typeStr, int n) {
  if (strcmp("float", typeStr) == 0) {
    compoundType = GL_FLOAT;
    baseType = GL_FLOAT;
    count = 1;
  } else if (strcmp("int", typeStr) == 0) {
    baseType = GL_INT;
    compoundType = GL_INT;
    count = 1;
  } else if (strcmp("vec", typeStr) == 0) {
    baseType = GL_FLOAT;
    if (n == 2) {
      compoundType = GL_FLOAT_VEC2;
      count = 2;
    } else if (n == 3) {
      compoundType = GL_FLOAT_VEC3;
      count = 3;
    } else if (n == 4) {
      compoundType = GL_FLOAT_VEC4;
      count = 4;
    } else {
      assert(0);
    }
  } else if (strcmp("ivec", typeStr) == 0) {
    baseType = GL_INT;
    if (n == 2) {
      compoundType = GL_INT_VEC2;
      count = 2;
    } else if (n == 3) {
      compoundType = GL_INT_VEC3;
      count = 3;
    } else if (n == 4) {
      compoundType = GL_INT_VEC4;
      count = 4;
    } else {
      assert(0);
    }
  } else if (strcmp("bvec", typeStr) == 0) {
    baseType = GL_BOOL;
    if (n == 2) {
      compoundType = GL_BOOL_VEC2;
      count = 2;
    } else if (n == 3) {
      compoundType = GL_BOOL_VEC3;
      count = 3;
    } else if (n == 4) {
      compoundType = GL_BOOL_VEC4;
      count = 4;
    }
  } else if (strcmp("mat", typeStr) == 0) {
    baseType = GL_FLOAT;
    if (n == 2) {
      compoundType = GL_FLOAT_MAT2;
      count = 4;
    } else if (n == 3) {
      compoundType = GL_FLOAT_MAT3;
      count = 9;
    } else if (n == 4) {
      compoundType = GL_FLOAT_MAT4;
      count = 16;
    } else {
      assert(0);
    }
  } else if (strcmp("sampler2D", typeStr) == 0) {
    baseType = GL_SAMPLER_2D;
    compoundType = GL_SAMPLER_2D;
    count = 1;
  } else if (strcmp("sampler3D", typeStr) == 0) {
    baseType = GL_SAMPLER_3D_OES;
    compoundType = GL_SAMPLER_3D_OES;
    count = 1;
  } else if (strcmp("samplerCube", typeStr) == 0) {
    baseType = GL_SAMPLER_CUBE;
    compoundType = GL_SAMPLER_CUBE;
    count = 1;
  } else {
    // @todo { matrix type }
    // Unknown type.
    assert(0);
  }
}

// sizeof(OpenGL type)
static int GetGLTypeSize(GLenum type) {
  int sz = 0;
  switch (type) {
  case GL_FLOAT:
    sz = sizeof(float);
    break;
  case GL_INT:
    sz = sizeof(int);
    break;
  case GL_UNSIGNED_INT:
    sz = sizeof(int);
    break;
  case GL_FLOAT_VEC2:
    sz = sizeof(float) * 2;
    break;
  case GL_FLOAT_VEC3:
    sz = sizeof(float) * 3;
    break;
  case GL_FLOAT_VEC4:
    sz = sizeof(float) * 4;
    break;
  case GL_INT_VEC2:
    sz = sizeof(int) * 2;
    break;
  case GL_INT_VEC3:
    sz = sizeof(int) * 3;
    break;
  case GL_INT_VEC4:
    sz = sizeof(int) * 4;
    break;
  case GL_FLOAT_MAT2:
    sz = sizeof(float) * 2 * 2;
    break;
  case GL_FLOAT_MAT3:
    sz = sizeof(float) * 3 * 3;
    break;
  case GL_FLOAT_MAT4:
    sz = sizeof(float) * 4 * 4;
    break;
  case GL_SAMPLER_2D:
    sz = sizeof(int);
    break;
  case GL_SAMPLER_3D:
    sz = sizeof(int);
    break;
  default:
    sz = 0;
    break;
  }

  // @todo { bool, matrix }

  assert(sz && "Unsupported type specified.");
  return sz;
}

///< Generate unique temp filename. User must free memory of returned pointer
/// when it become unused.
static char *GenerateUniqueTempFilename() {

#if defined(_WIN32)

#ifdef LSGL_ENABLE_MPI
  char prefix[1024];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  sprintf(prefix, "lsgl_shader_source_%08d", rank);
#else
  char prefix[] = "lsgl_shader_source";
#endif

  // The string returned by tempnam() is allocated using malloc(3) and hence
  // should be freed by free(3).
  char *basename = tempnam(NULL, prefix);
  if (!basename) {
    if (errno == EINVAL) {
      fprintf(stderr, "Bad template parameter.\n");
      return NULL;
    } else if (errno == EEXIST) {
      fprintf(stderr, "Out of unique filename.\n");
      return NULL;
    } else {
      fprintf(stderr, "Failed to create unique name.\n");
      return NULL;
    }
  }

  char *name = strdup(basename);

  free(basename);

#else

// char *name = tempnam(NULL, prefix);

#ifdef LSGL_ENABLE_MPI
  char basename[1024];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  sprintf(basename, "lsgl_shader_source_%08d_XXXXXX", rank);
#else
  char basename[] = "lsgl_shader_source_XXXXXX";
#endif

  int fd = mkstemp(basename); // prefix will be overwritten with actual filename
  if (fd == -1) {
    fprintf(stderr, "Failed to create unique filename.\n");
    return NULL;
  }
  close(fd);

  char *name = strdup(basename);

#endif

  return name;
}

namespace lsgl {

void initialize_random() {
#ifdef _OPENMP
  assert(omp_get_max_threads() < 512);

  for (int i = 0; i < omp_get_max_threads(); i++) {
    tinymt64[i].mat1 = 0xfa051f40;
    tinymt64[i].mat2 = 0xffd0fff4;
    tinymt64[i].tmat = 0x58d02ffeffbfffbcULL;
    tinymt64_init(&tinymt64[i], (uint64_t)i);
  }
#else
  tinymt64[0].mat1 = 0xfa051f40;
  tinymt64[0].mat2 = 0xffd0fff4;
  tinymt64[0].tmat = 0x58d02ffeffbfffbcULL;
  tinymt64_init(&tinymt64[0], (uint64_t)1234);
#endif
}

float randomreal(int thread_id) {
  return (float)tinymt64_generate_double(&tinymt64[thread_id]);
}

} // namespace lsgl

//
// Buffer
//
Buffer::Buffer()
    : size_(0), usage_(GL_STATIC_DRAW), retained_(false), ptr_(NULL) {
  data_.clear();
}

Buffer::~Buffer() {}

void Buffer::Data(size_t size, const void *data, GLenum usage) {
  if (size == 0) {
    // Explicitly release memory of 'data_'.
    std::vector<GLubyte>().swap(data_);
    return;
  }

  size_ = size;
  usage_ = usage;
  data_.resize(size);
  memcpy(&data_[0], data, size_);
  retained_ = false;
}

void Buffer::SubData(unsigned int offset, size_t size, const void *data) {
  assert(retained_ == false);
  assert(offset + size < data_.size());
  memcpy(&data_[offset], data, size);
}

void Buffer::Retain(size_t size, const void *data, GLenum usage) {
  if (size == 0) {
    ptr_ = NULL;
    return;
  }

  size_ = size;
  usage_ = usage;
  ptr_ = reinterpret_cast<const unsigned char *>(data);
  retained_ = true;
}

//
// Shader
//
Shader::Shader()
    : compiled_(false), isBinary_(false), isAttached_(false), handle_(NULL),
      filename_("") {
  method_.shaderEvalFunc = NULL;
  method_.shaderInfoFunc = NULL;
  method_.shaderInitFunc = NULL;
}

Shader::~Shader() {
#if defined(_WIN32)
  Release();
#else
  Release();
#endif
}

void Shader::Source(GLsizei count, const GLchar **string, const GLint *length) {
  // clear out any previous source
  source_.clear();

  // concatenate all string buffers together
  for (int t = 0; t < count; t++) {
    if ((length == NULL) || (length[t] < 0)) {
      source_.append(string[t]);
    } else {
      source_.append(string[t], length[t]);
    }
  }

  //    printf("shader source: [%s]\n", source_.c_str());
}

bool Shader::Compile() {
  compiled_ = DoCompile();
  return compiled_;
}

bool Shader::Release() {
  if (compiled_) {
    compiled_ = false;

#if defined(_WIN32)
    if (handle_) {
      BOOL ret = FreeLibrary(reinterpret_cast<HMODULE>(handle_));
      assert(ret == TRUE);
      return (ret == TRUE) ? true : false;
    }

    if (!filename_.empty()) {
      int ret = unlink(filename_.c_str());
      if (ret == -1) {
        fprintf(stderr, "[LSGL] Failed to delete file: %s\n",
                filename_.c_str());
      }

      filename_ = std::string("");
    }

#else
    if (handle_) {
      int ret = dlclose(handle_);
      assert(ret == 0);
      return (ret == 0) ? true : false;
    }

    if (!filename_.empty()) {
      int ret = unlink(filename_.c_str());
      if (ret == -1) {
        fprintf(stderr, "[LSGL] Failed to delete file: %s\n",
                filename_.c_str());
      }

      filename_ = std::string("");
    }
#endif
  }

  return true;
}

bool Shader::LoadShaderBinary(std::string &filename) {

  void *handle = NULL;

#ifdef _WIN32
  HMODULE module = LoadLibrary(filename.c_str());
  if (module == NULL) {
    fprintf(stderr, "[LSGL] Cannot find/open shader file: %s\n",
            filename.c_str());
    fprintf(stderr, "[LSGL] Err = %d\n", GetLastError());
    return false;
  }

  // Find shader body function
  method_.shaderEvalFunc = GetProcAddress(module, "shader");
  if (method_.shaderEvalFunc == NULL) {
    fprintf(stderr, "[LSGL] Cannot find shader body function from: %s\n",
            filename.c_str());
    FreeLibrary(module);
    return false;
  }

  // Find shader info function
  method_.shaderInfoFunc = GetProcAddress(module, "shader_info");
  if (method_.shaderInfoFunc == NULL) {
    fprintf(stderr, "[LSGL] Cannot find shader info function from: %s\n",
            filename.c_str());
    FreeLibrary(module);
    return false;
  }

  // printf("[LSGL] load shader binary OK\n");

  handle_ = reinterpret_cast<void *>(handle);
  filename_ = filename;

#else
  std::string filepath = filename;
  handle = dlopen(filepath.c_str(), RTLD_NOW);

  if (handle != NULL) {
    // Will be safe to delete .so file after dlopen().
    unlink(filename.c_str());
  } else {
    if ((filename.size() > 1) &&
        ((filename[0] != '/') || (filename[0] != '.'))) {
      // try to load from current path(this might have security risk?).
      filepath = std::string("./") + filename;
      handle = dlopen(filepath.c_str(), RTLD_NOW);

      // Will be safe to delete .so file after dlopen().
      unlink(filename.c_str());

      if (handle == NULL) {
        fprintf(stderr, "[LSGL] Cannot find/open shader file: %s\n",
                filepath.c_str());
      }
    } else {
      fprintf(stderr, "[LSGL] Cannot find/open shader file: %s\n",
              filename.c_str());
      return false;
    }
  }

  // Find shader body function
  method_.shaderEvalFunc = dlsym(handle, "shader");
  if (method_.shaderEvalFunc == NULL) {
    fprintf(stderr, "[LSGL] Cannot find shader body function from: %s\n",
            filename.c_str());
    dlclose(handle);
    return false;
  }

  // Find shader info function
  method_.shaderInfoFunc = dlsym(handle, "shader_info");
  if (method_.shaderInfoFunc == NULL) {
    fprintf(stderr, "[LSGL] Cannot find shader info function from: %s\n",
            filepath.c_str());
    dlclose(handle);
    return false;
  }

  // printf("[LSGL] load shader binary OK\n");

  // Store handle & filename for later use.
  handle_ = handle;
  filename_ = filepath;
#endif

  // Binary shader was already compiled.
  compiled_ = true;

  return true;
}

void Shader::BuildUniformInfo(std::vector<Uniform> &uniforms,
                              std::vector<UniformLocation> &uniformLocations) {
  assert(method_.shaderInfoFunc);

  uniforms.clear();
  uniformLocations.clear();

  ShaderInfoFuncPtr func =
      reinterpret_cast<ShaderInfoFuncPtr>(method_.shaderInfoFunc);

  FragmentConfig config;
  func(&config);

  // printf("num uniforms = %d\n", config.numUniforms);

  for (int i = 0; i < config.numUniforms; i++) {
    GLSLUniformInfo info = config.uniformInfos[i];
    // printf("  [%d] name = %s, ty = %s:%d, qual = %d\n", i, info.name,
    // info.type.name, info.type.n, info.qualifier);

    // @todo { Array type. }
    UniformLocation uniformLocation(info.name, 0, i);
    GLenum compoundTy, baseTy;
    int count;
    GetGLType(compoundTy, baseTy, count, info.type.name, info.type.n);
    Uniform uniform(compoundTy, info.name, 1);

    // The spec says...
    //   all active user-defined uniform variables belonging to program will be
    // initialized to 0
    // so fill with zero.

    int tySize = GetGLTypeSize(compoundTy);
    uniform.data.resize(tySize);
    memset(&uniform.data.at(0), 0, sizeof(tySize));

    uniforms.push_back(uniform);
    uniformLocations.push_back(uniformLocation);
  }
}

void Shader::BuildVaryingInfo(std::vector<Varying> &varyings,
                              std::vector<VaryingLocation> &varyingLocations) {
  assert(method_.shaderInfoFunc);

  varyings.clear();
  varyingLocations.clear();

  ShaderInfoFuncPtr func =
      reinterpret_cast<ShaderInfoFuncPtr>(method_.shaderInfoFunc);

  FragmentConfig config;
  func(&config);

  // printf("num varyings = %d\n", config.numVaryings);

  // location 0 is reserved for position
  VaryingLocation varyingLocation("***[pos:reserved]***", 0, (unsigned int)-1);
  varyingLocations.push_back(varyingLocation);

  for (int i = 0; i < config.numVaryings; i++) {
    GLSLVaryingInfo info = config.varyingInfos[i];
    // printf("  [%d] name = %s, ty = %s:%d, qual = %d\n", i, info.name,
    // info.type.name, info.type.n, info.qualifier);

    // @todo { Array type. }
    VaryingLocation varyingLocation(info.name, 0, i);
    GLenum compoundTy, baseTy;
    int count;
    GetGLType(compoundTy, baseTy, count, info.type.name, info.type.n);
    Varying varying(baseTy, count, info.name, 1);

    // Allocate memory and fill it with zero.
    int tySize = GetGLTypeSize(compoundTy);
    varying.data.resize(tySize);
    memset(&varying.data.at(0), 0, sizeof(tySize));

    varyings.push_back(varying);
    varyingLocations.push_back(varyingLocation);
  }

  // Add LSGL extension variables.
  {
    int index = config.numVaryings;
    assert(index < kMaxVertexAttribs);
    VaryingLocation pointSizeVaryingLocation("lsgl_PointSize", 0, index);
    Varying pointSizeVarying(GL_FLOAT, 1, "lsgl_PointSize", 1);

    varyings.push_back(pointSizeVarying);
    varyingLocations.push_back(pointSizeVaryingLocation);
  }

  {
    int index = config.numVaryings + 1;
    assert(index < kMaxVertexAttribs);
    VaryingLocation lineSizeVaryingLocation("lsgl_LineWidth", 0, index);
    Varying lineSizeVarying(GL_FLOAT, 1, "lsgl_LineWidth", 1);

    varyings.push_back(lineSizeVarying);
    varyingLocations.push_back(lineSizeVaryingLocation);
  }
}

VertexShader::VertexShader() {}

VertexShader::~VertexShader() {}

bool VertexShader::DoCompile() { return false; /* not suppored */ }

FragmentShader::FragmentShader() {}

FragmentShader::~FragmentShader() {}

bool FragmentShader::DoCompile() {

  // @fixme { Assume this function is not re-entrant. }
  static int sCounter = 0;

  const char *src = GetSource().c_str();

  //
  // 1. Dump shader code to a file.
  //
  char *tempFilenameStr = GenerateUniqueTempFilename();
  assert(tempFilenameStr);

  std::string tempFilename = std::string(tempFilenameStr);
  tempFilename += std::string(".frag"); // add suffix

  int status = unlink(tempFilenameStr);
  if (status == -1) {
    perror("unlink");
  }
  free(tempFilenameStr);

  FILE *fp = fopen(tempFilename.c_str(), "w");
  if (!fp) {
    fprintf(stderr, "Can't open file to write.");
    return false;
  }

  // printf("Write to file: %s\n", tempFilename.c_str());

  fwrite(src, GetSource().size(), 1, fp);
  fclose(fp);

  //
  // 2. Invoke offline compiler
  //
  const char *glslc = "glslc";
  char *glslc_path = getenv("GLSL_COMPILER");
  if (glslc_path) {
    glslc = glslc_path;
    printf("[LSGL] glslc = %s\n", glslc);
  }

  // Generate unique filename.
  char outputFilename[2048];
#ifdef LSGL_ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  sprintf(outputFilename, "shader_%04d_%08d.so", sCounter, rank);
#else
  sprintf(outputFilename, "shader_%04d.so", sCounter);
#endif

  sCounter++;

  char cmd[4096];
  sprintf(cmd, "%s -o %s %s", glslc, outputFilename, tempFilename.c_str());
  printf("[LSGL] cmd: %s\n", cmd);


  // See popen(3) manual for details why calling fflush(NULL) here
  fflush(NULL);

#if defined(_WIN32)
  FILE *pfp = _popen(cmd, "r");
#else
  FILE *pfp = popen(cmd, "r");
#endif

  if (!pfp) {
    fprintf(stderr, "[LSGL] Failed to open pipe.\n");
    perror("popen");

    status = unlink(tempFilename.c_str());
    if (status == -1) {
      perror("unlink");
    }

    return false;
  }

  char buf[4096];
  while (fgets(buf, 4095, pfp) != NULL) {
    printf("[glslc] %s", buf);
  }

  status = unlink(tempFilename.c_str());
  if (status == -1) {
    perror("unlink");
  }

#if defined(_WIN32)
  status = _pclose(pfp);
  if (status == -1) {
    fprintf(stderr, "[LSGL] Failed to close pipe.\n");
    return false;
  }
#else
  status = pclose(pfp);
  if (status == -1) {
    fprintf(stderr, "[LSGL] Failed to close pipe.\n");
    return false;
  }
#endif

  //
  // 3. Check if compiled shader exists.
  //
  std::string shaderBinaryFilename(outputFilename);
  std::ifstream ifile(shaderBinaryFilename.c_str());
  if (!ifile) {
    fprintf(stderr, "[LSGL] Failed to compile shader.\n");
    return false;
  }

  //
  // 4. load compiled shader binary(in file)
  //
  bool ret = LoadShaderBinary(shaderBinaryFilename);
  return ret;
}

bool FragmentShader::PrepareEval(
    FragmentState &fragmentState, ShadingState &shadingState,
    const std::vector<Uniform> &uniforms, const std::vector<Varying> &varyings,
    const std::vector<VaryingLocation> &varyingLocations,
    const std::vector<VertexAttribute> &vertexAttributes,
    const Context &ctx) const {
  FragmentState *state = new FragmentState();

  //
  // Fill uniforms
  //
  for (size_t i = 0; i < uniforms.size(); i++) {
    if (uniforms[i].data.empty()) {
      state->uniforms[i].data = NULL;
    } else {
      unsigned char *ptr = const_cast<unsigned char *>(&uniforms[i].data.at(0));
#ifdef LSGL_OPTIMIZE_GLSL

      state->textures[i] = 0; // Init as invalid

      if (uniforms[i].type == GL_SAMPLER_3D) {
        // Store Texture3D value diretly.
        int handle = *reinterpret_cast<int *>(ptr);
        if (ctx.resourceManager_.IsValidTexture(handle) == true) {
          const Texture *tex = ctx.resourceManager_.GetTexture(handle);
          state->textures[i] = reinterpret_cast<uintptr_t>(tex->GetTexture3D());
        }
      } else if (uniforms[i].type == GL_SAMPLER_2D ||
                 uniforms[i].type == GL_SAMPLER_CUBE) {

        // Store pointer value as index;
        int handle = *reinterpret_cast<int *>(ptr);
        if (ctx.resourceManager_.IsValidTexture(handle) == true) {
          const Texture *tex = ctx.resourceManager_.GetTexture(handle);
          state->textures[i] = reinterpret_cast<uintptr_t>(tex);
        }
#else
      if (uniforms[i].type == GL_SAMPLER_2D ||
          uniforms[i].type == GL_SAMPLER_3D ||
          uniforms[i].type == GL_SAMPLER_CUBE) {
        state->textures[i] = *(reinterpret_cast<int *>(ptr));
#endif
      } else {
        // Pass uniform data as opeque pointer.
        state->uniforms[i].data = ptr;
      }
    }
  }

  //
  // Bind slot of varyings
  //
  shadingState.varyingConnections.clear();

  //
  // dst index                           src index
  //
  // GLSL(VaryingLocations)              GLES(Attrib)
  // +--------------------------+        +-----------------------+
  // | 0(reserved for position) | <----> | 0(basically position) |
  // +--------------------------+        +-----------------------+
  // | 1                        | <-+    | 1                     |
  // +--------------------------+   |    +-----------------------+
  // | 2                        |   +--> | 2                     |
  // +--------------------------+        +-----------------------+
  // .                          .        .                       .
  // .                          .        .                       .
  // +--------------------------+        +-----------------------+
  // | n                        |        | n                     |
  // +--------------------------+        +-----------------------+

  // varyingLocations[0] is reserved for position.
  for (size_t i = 1; i < varyingLocations.size(); i++) {
    const VaryingLocation &varyingLocation = varyingLocations[i];
    // printf("[%ld] varying.idx = %d\n", i, varyingLocation.index);
    if (varyingLocation.index == (unsigned int)-1) {
      // ??? -1 is reserved for position.
      continue;
    }

    assert(i <= kMaxVertexAttribs);

    const VertexAttribute &vertexAttribute = vertexAttributes[i];
    if (!vertexAttribute.enabled) {
      // not bound.
      continue;
    }

    const Varying &varying = varyings[varyingLocation.index];

    // printf("  src = %ld, dst = %d\n", i, varyingLocation.index);
    // printf("  attrib : ty = %d, n = %d, ptr = %p\n", vertexAttribute.type,
    // vertexAttribute.size, vertexAttribute.ptr);
    // printf("  varying: ty = %d, n = %d\n", varying.type, varying.count);
    // printf("  va.stride %d\n", vertexAttribute.stride);
    // printf("  varying.data.size = %ld\n", varying.data.size());

    if ((vertexAttribute.type == varying.type) &&
        (vertexAttribute.size == varying.count)) {
      // printf("  ==> type signature OK!\n");
    }

    const unsigned char *ptr = NULL;
    if (vertexAttribute.ptr == NULL) {
      // Use vertex buffer

      // get buffer pointer from handle
      int handle = vertexAttribute.handle;
      assert((handle > 0) && "Unexpected behavior");

      const Buffer *buf = ctx.resourceManager_.GetBuffer(handle);
      assert(buf && "Unexpected behavior");

      ptr = buf->GetData();

    } else {
      assert(0 && "TODO: Vertex attribute should be specified by vertex buffer "
                  "in current implelementation.");
    }

    assert(ptr);

    VaryingConnection varyingConn;
    varyingConn.srcIndex = i;
    varyingConn.dstIndex = varyingLocation.index;
    varyingConn.ptr = ptr;
    varyingConn.size = varying.data.size();
    varyingConn.stride = vertexAttribute.stride;

    // Don't forget to clear() in the beginning of PrepareEval();
    shadingState.varyingConnections.push_back(varyingConn);
    shadingState.varyings.push_back(varying);
  }

  // if (state_) {
  //    delete state_;
  //}
  // state_ = reinterpret_cast<unsigned char*>(state); // Reference
  // state_ = *state; // copy
  fragmentState = *state;
  delete state;

  return true;
}

namespace {

// f' = (1 - u - v) f0 + u f1 + v f2
inline float LerpFloat(float f0, float f1, float f2, float u, float v) {
  return (1.0f - u - v) * f0 + u * f1 + v * f2;
}

} // namespace

bool FragmentShader::Eval(GLfloat fragColor[4], FragmentState &fragmentState,
                          ShadingState &shadingState,
                          const std::vector<VertexAttribute> &vertexAttributes,
                          const GLfloat fragCoord[4],
                          const IntersectionState &isectState,
                          const CameraInfo &cameraInfo, int threadID) const {
  if (!method_.shaderEvalFunc) {
    return false;
  }

  // Copy state for shader evaluation. this makes 'state' variable thread-safe
  // @todo { uniform is invariant. just copy varying variables. }
  // FragmentState state = *(reinterpret_cast<FragmentState*>(state_));
  FragmentState state;
  // memcpy(&state, &state_, sizeof(FragmentState));
  memcpy(&state, &fragmentState, sizeof(FragmentState));

  //
  // Allocate zero filled reserved area for unresolved varying variables.
  // this buffer shoud have enough size to store all GLSL types
  // (i.e, up to mat44)
  //
  GLfloat varyingBuffer[MAX_FRAGMENT_VARYING_VARIABLES][16];
  memset(varyingBuffer, 0,
         sizeof(GLfloat) * MAX_FRAGMENT_VARYING_VARIABLES * 16);

  //
  // Initialize varying to point zeroBuffer so that unresolved varying
  // varyings will have zero value in the shader.
  //
  for (int i = 0; i < MAX_FRAGMENT_VARYING_VARIABLES; i++) {
    state.varyings[i].data =
        reinterpret_cast<unsigned char *>(&varyingBuffer[i][0]);
  }

  ShaderEvalFuncPtr eval =
      reinterpret_cast<ShaderEvalFuncPtr>(method_.shaderEvalFunc);

  Fragment frag;
  frag.fragDiscarded = 0;

  frag.fragCoord[0] = fragCoord[0];
  frag.fragCoord[1] = fragCoord[1];
  frag.fragCoord[2] = fragCoord[2];
  frag.fragCoord[3] = fragCoord[3];

  frag.position[0] = isectState.position[0];
  frag.position[1] = isectState.position[1];
  frag.position[2] = isectState.position[2];
  frag.normal[0] = isectState.normal[0];
  frag.normal[1] = isectState.normal[1];
  frag.normal[2] = isectState.normal[2];
  frag.geometricNormal[0] = isectState.geometricNormal[0];
  frag.geometricNormal[1] = isectState.geometricNormal[1];
  frag.geometricNormal[2] = isectState.geometricNormal[2];
  frag.tangent[0] = isectState.tangent[0];
  frag.tangent[1] = isectState.tangent[1];
  frag.tangent[2] = isectState.tangent[2];
  frag.binormal[0] = isectState.binormal[0];
  frag.binormal[1] = isectState.binormal[1];
  frag.binormal[2] = isectState.binormal[2];
  frag.indir[0] = isectState.raydir[0];
  frag.indir[1] = isectState.raydir[1];
  frag.indir[2] = isectState.raydir[2];
  frag.barycentric[0] = isectState.u;
  frag.barycentric[1] = isectState.v;
  frag.px = isectState.px;
  frag.py = isectState.py;
  frag.raydepth = isectState.raydepth;
  frag.doubleSided = isectState.doubleSided;
  frag.rayattrib = isectState.rayattrib;
  frag.prev_node = isectState.prev_node;
  frag.prev_prim_id = isectState.prev_prim_id;
  frag.prev_hit_t = isectState.prev_hit_t;
  frag.prev_hit_normal[0] = isectState.prev_hit_normal[0];
  frag.prev_hit_normal[1] = isectState.prev_hit_normal[1];
  frag.prev_hit_normal[2] = isectState.prev_hit_normal[2];
  frag.threadID = threadID;

  frag.cameraFrame[0][0] = cameraInfo.frame[0][0];
  frag.cameraFrame[0][1] = cameraInfo.frame[0][1];
  frag.cameraFrame[0][2] = cameraInfo.frame[0][2];
  frag.cameraFrame[1][0] = cameraInfo.frame[1][0];
  frag.cameraFrame[1][1] = cameraInfo.frame[1][1];
  frag.cameraFrame[1][2] = cameraInfo.frame[1][2];
  frag.cameraFrame[2][0] = cameraInfo.frame[2][0];
  frag.cameraFrame[2][1] = cameraInfo.frame[2][1];
  frag.cameraFrame[2][2] = cameraInfo.frame[2][2];
  frag.cameraFov = cameraInfo.fov;

  // Setup function pointer
  frag.texture2D = texture2D;
  frag.texture3D = texture3D;
  frag.shadow = shadow;
  frag.trace = trace;
  frag.random = randomreal;

  unsigned int f0 = isectState.f0;
  unsigned int f1 = isectState.f1;
  unsigned int f2 = isectState.f2;

  // Fill varying variables.
  // printf("index = %d\n", index);
  // printf("f = %d, %d, %d\n", f0, f1, f2);
  // printf("conns = %d\n", shadingState.varyingConnections.size());
  for (int i = 0; i < shadingState.varyingConnections.size(); i++) {
    const VaryingConnection &varyingConn = shadingState.varyingConnections[i];
    const Varying &varying = shadingState.varyings[i];

    // @fixme { Currently only support for GL_FLOAT varyings. }
    // Assume GL_TRIANGLES
    int elems = varying.data.size() / sizeof(GL_FLOAT);
    const GLfloat *f0ptr = reinterpret_cast<const GLfloat *>(
        &varyingConn.ptr[f0 * varyingConn.stride]);
    const GLfloat *f1ptr = reinterpret_cast<const GLfloat *>(
        &varyingConn.ptr[f1 * varyingConn.stride]);
    const GLfloat *f2ptr = reinterpret_cast<const GLfloat *>(
        &varyingConn.ptr[f2 * varyingConn.stride]);
    // GLfloat* dst =
    // reinterpret_cast<GLfloat*>(&(shadingState.varyings[i].data.at(0)));

    // Use local(stack) buffer
    GLfloat *dst = reinterpret_cast<GLfloat *>(&varyingBuffer[i][0]);

    // printf("dst = %p, f0 = %d\n", dst, f0);

    if ((f0 == f1) && (f1 == f2)) {
      // No interpolation required.
      for (int k = 0; k < elems; k++) {
        float f = f0ptr[k];
        dst[k] = f; // store to varying storage
      }
    } else {
      for (int k = 0; k < elems; k++) {
        float f =
            LerpFloat(f0ptr[k], f1ptr[k], f2ptr[k], isectState.u, isectState.v);
        dst[k] = f; // store to varying storage
      }
    }

    // link to varying storage so that the shader can see (interpolated) varying
    // value.
    state.varyings[varyingConn.dstIndex].data =
        reinterpret_cast<unsigned char *>(dst);
  }

  // frag will be modified after function call.
  // Assume state will not be modified.
  // printf("shader func = %p\n", eval);
  // printf("  [BEFORE] frag.col = %f, %f, %f, %f\n",
  //    frag.fragColor[0],
  //    frag.fragColor[1],
  //    frag.fragColor[2],
  //    frag.fragColor[3]);

  eval(&frag, &state);

  if (frag.fragDiscarded) {
    return false;
  }

  // printf("  [AFTER ] frag.col = %f, %f, %f, %f\n",
  //   frag.fragColor[0],
  //    frag.fragColor[1],
  //    frag.fragColor[2],
  //    frag.fragColor[3]);

  fragColor[0] = frag.fragColor[0];
  fragColor[1] = frag.fragColor[1];
  fragColor[2] = frag.fragColor[2];
  fragColor[3] = frag.fragColor[3];

  return true;
}

//
// Texture(X)
//
Texture::Texture()
    : texture2D_(NULL), texture3D_(NULL), retained_(false),
      sparseVolumeAccel_(NULL), sparseVolume_(NULL), isSparse_(false),
      minFiltering_(true), magFiltering_(true) {
  doRemap_[0] = false;
  doRemap_[1] = false;
  doRemap_[2] = false;
}

Texture::~Texture() { Free(); }

void Texture::Data2D(const void *data, GLuint width, GLuint height, int compos,
                     GLenum type) {
  // first free any current texture
  Free();

  int tsize = 1; // UCHAR
  if (type == GL_FLOAT) {
    tsize = sizeof(float);
  }

  // calculate texture size in bytes
  size_ = width * height * compos * tsize;

  // allocate and copy the texture data
  data_.reserve(size_);
  if (data) {
    memcpy(&data_[0], data, size_);
  } else {
    // Zero fill.
    memset(&data_[0], 0, size_);
  }

  render::Texture2D::Format format = render::Texture2D::FORMAT_BYTE;
  if (type == GL_FLOAT) {
    format = render::Texture2D::FORMAT_FLOAT32;
  } else if (type == GL_DOUBLE) {
    format = render::Texture2D::FORMAT_FLOAT64;
  }

  // create a new texture object
  texture2D_ = new render::Texture2D(&data_[0], width, height, compos, format);
  texture3D_ = NULL;
}

void Texture::Data3D(const void *data, GLuint width, GLuint height,
                     GLuint depth, int compos, GLenum type) {
  // first free any current texture
  Free();

  // calculate texture size in bytes
  int data_type = LSGL_RENDER_TEXTURE3D_FORMAT_INVALID;
  if (type == GL_DOUBLE) {
    size_ = width * height * depth * sizeof(double) * compos;
    data_type = LSGL_RENDER_TEXTURE3D_FORMAT_DOUBLE;
  } else if (type == GL_FLOAT) {
    size_ = width * height * depth * sizeof(float) * compos;
    data_type = LSGL_RENDER_TEXTURE3D_FORMAT_FLOAT;
  } else if (type == GL_UNSIGNED_BYTE) {
    size_ = width * height * depth * sizeof(unsigned char) * compos;
    data_type = LSGL_RENDER_TEXTURE3D_FORMAT_BYTE;
  } else {
    assert(0 && "Unknown format.");
  }

  // allocate and copy the texture data
  data_.reserve(size_);
  if (data) {
    memcpy(&data_[0], data, size_);
  } else {
    // Zero fill.
    memset(&data_[0], 0, size_);
  }

  // create a new texture object
  texture2D_ = NULL;
  texture3D_ = (Texture3D *)malloc(sizeof(Texture3D));
  assert(texture3D_);

  texture3D_->image = reinterpret_cast<unsigned char *>(&data_[0]);
  texture3D_->width = width;
  texture3D_->height = height;
  texture3D_->depth = depth;
  texture3D_->components = compos;
  texture3D_->data_type = data_type;
  texture3D_->wrapS = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
  texture3D_->wrapT = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
  texture3D_->wrapR = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
}

void Texture::Retain3D(const void *data, GLuint width, GLuint height,
                       GLuint depth, int compos, GLenum type) {
  // first free any current texture
  Free();

  // calculate texture size in bytes
  int data_type = LSGL_RENDER_TEXTURE3D_FORMAT_INVALID;
  if (type == GL_DOUBLE) {
    size_ = width * height * depth * sizeof(double) * compos;
    data_type = LSGL_RENDER_TEXTURE3D_FORMAT_DOUBLE;
  } else if (type == GL_FLOAT) {
    size_ = width * height * depth * sizeof(float) * compos;
    data_type = LSGL_RENDER_TEXTURE3D_FORMAT_FLOAT;
  } else if (type == GL_UNSIGNED_BYTE) {
    size_ = width * height * depth * sizeof(unsigned char) * compos;
    data_type = LSGL_RENDER_TEXTURE3D_FORMAT_BYTE;
  } else {
    assert(0 && "Unknown format.");
  }

  // create a new texture object. Texture data is sent to internal object
  // directly(zero-copy).
  texture2D_ = NULL;
  texture3D_ = (Texture3D *)malloc(sizeof(Texture3D));
  assert(texture3D_);

  texture3D_->image =
      reinterpret_cast<unsigned char *>(const_cast<void *>(data));
  texture3D_->width = width;
  texture3D_->height = height;
  texture3D_->depth = depth;
  texture3D_->components = compos;
  texture3D_->data_type = data_type;
  texture3D_->wrapS = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
  texture3D_->wrapT = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
  texture3D_->wrapR = LSGL_RENDER_TEXTURE3D_WRAP_REPEAT;
  retained_ = true;
}

void Texture::Free() {

  // release the texture object
  delete texture2D_;
  texture2D_ = NULL;

  delete texture3D_;
  texture3D_ = NULL;

  // then free our internal copy of the texture data
  // note: it is also safe to call swap() for retained texture.
  std::vector<GLubyte>().swap(data_);
}

void Texture::SubImage3D(GLint xoffset, GLint yoffset, GLint zoffset,
                         GLsizei width, GLsizei height, GLsizei depth,
                         int compos, GLenum type, const GLvoid *data) {

  if (isSparse_) {

    if (regionList_.size() == 1) {
      compos_ = compos;
      type_ = type;
    } else {
      if (compos_ != compos) {
        return;
      }

      if (type_ != type) {
        return;
      }
    }

    int dataSize = -1;
    if (type_ == GL_UNSIGNED_BYTE) {
      dataSize = 1;
    } else if (type_ == GL_FLOAT) {
      dataSize = 4;
    } else if (type_ == GL_DOUBLE) {
      dataSize = 8;
    } else {
      assert(0 && "Unsupported data type.");
    }

    // Find region.
    // @todo { optimize region search. }
    for (size_t i = 0; i < regionList_.size(); i++) {
      if ((xoffset == regionList_[i].offset[0]) &&
          (yoffset == regionList_[i].offset[1]) &&
          (zoffset == regionList_[i].offset[2]) &&
          (width == regionList_[i].extent[0]) &&
          (height == regionList_[i].extent[1]) &&
          (depth == regionList_[i].extent[2])) {
        // Gotcha!
        delete[] regionList_[i].data;
        regionList_[i].data =
            new unsigned char[width * height * depth * dataSize * compos];
        memcpy(reinterpret_cast<void *>(regionList_[i].data), data,
               width * height * depth * dataSize * compos);
      }
    }
  } else {
    // SubImage3D for non-sparse texture is not supported at this time.
    return;
  }
}

bool Texture::SetRemapTable(GLenum coord, GLsizei size, const GLfloat *coords) {
  int idx = -1;
  if (coord == GL_COORDINATE_X) {
    idx = 0;
  } else if (coord == GL_COORDINATE_Y) {
    idx = 1;
  } else if (coord == GL_COORDINATE_Z) {
    idx = 2;
  } else {
    // ???
    return false;
  }

  remapTable_[idx].resize(size);
  for (size_t i = 0; i < size; i++) {
    remapTable_[idx][i] = coords[i];
  }
  doRemap_[idx] = true;

  return true;
}

void Texture::SubImage3DRetain(GLint xoffset, GLint yoffset, GLint zoffset,
                               GLsizei width, GLsizei height, GLsizei depth,
                               int compos, GLenum type, const GLvoid *data) {

  if (isSparse_) {

    if (regionList_.size() == 1) {
      compos_ = compos;
      type_ = type;
    } else {
      if (compos_ != compos) {
        return;
      }

      if (type_ != type) {
        return;
      }
    }

    // Find region.
    // @todo { optimize region search. }
    for (size_t i = 0; i < regionList_.size(); i++) {
      if ((xoffset == regionList_[i].offset[0]) &&
          (yoffset == regionList_[i].offset[1]) &&
          (zoffset == regionList_[i].offset[2]) &&
          (width == regionList_[i].extent[0]) &&
          (height == regionList_[i].extent[1]) &&
          (depth == regionList_[i].extent[2])) {
        // Just save an poineter(no local copy)
        regionList_[i].data =
            reinterpret_cast<unsigned char *>(const_cast<GLvoid *>(data));
      }
    }
  } else {
    // SubImage3DRetain for non-sparse texture is not supported at this time.
    return;
  }
}

void Texture::BuildSparseTexture() {
  if (!isSparse_) {
    return;
  }

  if (sparseVolumeAccel_) {
    return;
  }

  if (sparseVolume_) {
    return;
  }

  sparseVolume_ = new SparseVolume();
  sparseVolume_->components = 1; // @fixme

  if (type_ == GL_UNSIGNED_BYTE) {
    sparseVolume_->format = LSGL_RENDER_TEXTURE3D_FORMAT_BYTE;
  } else if (type_ == GL_FLOAT) {
    sparseVolume_->format = LSGL_RENDER_TEXTURE3D_FORMAT_FLOAT;
  } else if (type_ == GL_DOUBLE) {
    sparseVolume_->format = LSGL_RENDER_TEXTURE3D_FORMAT_DOUBLE;
  } else { // ???
    assert(0);
    sparseVolume_->format = 0;
  }

  // Find largest extent and make it global dim.
  int dim[3] = {0, 0, 0};
  for (size_t i = 0; i < regionList_.size(); i++) {
    dim[0] =
        (std::max)(dim[0], regionList_[i].offset[0] + regionList_[i].extent[0]);
    dim[1] =
        (std::max)(dim[1], regionList_[i].offset[1] + regionList_[i].extent[1]);
    dim[2] =
        (std::max)(dim[2], regionList_[i].offset[2] + regionList_[i].extent[2]);
  }

  sparseVolume_->globalDim[0] = dim[0];
  sparseVolume_->globalDim[1] = dim[1];
  sparseVolume_->globalDim[2] = dim[2];

  printf("[SparseTexture] Global dim: %d x %d x %d\n", dim[0], dim[1], dim[2]);

  for (size_t i = 0; i < regionList_.size(); i++) {
    VolumeBlock block;
    block.offset[0] = regionList_[i].offset[0];
    block.offset[1] = regionList_[i].offset[1];
    block.offset[2] = regionList_[i].offset[2];

    block.extent[0] = regionList_[i].extent[0];
    block.extent[1] = regionList_[i].extent[1];
    block.extent[2] = regionList_[i].extent[2];

    block.id = i;
    block.data = regionList_[i].data;

    sparseVolume_->blocks.push_back(block);
  }

  sparseVolumeAccel_ = new SparseVolumeAccel();
  sparseVolumeAccel_->Build(sparseVolume_);

  double bmin[3], bmax[3];
  sparseVolumeAccel_->BoundingBox(bmin, bmax);
  printf("[SparseTexture] bbox: (%f, %f, %f) - (%f, %f, %f)\n", bmin[0],
         bmin[1], bmin[2], bmax[0], bmax[1], bmax[2]);
}

//
// Renderbuffer
//
Renderbuffer::Renderbuffer() : buffer_(NULL) { Free(); }

Renderbuffer::~Renderbuffer() { Free(); }

bool Renderbuffer::Allocate(GLuint width, GLuint height, GLuint bytesPerPixel,
                            GLenum format) {
  // free any previous buffer first
  Free();

  // store parameters
  width_ = width;
  height_ = height;
  format_ = format;
  bytesPerPixel_ = bytesPerPixel;

  // attempt to allocate memory for a new buffer
  buffer_ = new char[width_ * height_ * bytesPerPixel_];
  return (buffer_ != NULL);
}

void Renderbuffer::Free() {
  // free buffer
  delete[] buffer_;
  buffer_ = NULL;

  // clear parameters
  width_ = 0;
  height_ = 0;
  bytesPerPixel_ = 0;
  format_ = 0;
}

//
// Framebuffer
//
Framebuffer::Framebuffer()
    : depthBuffer_(NULL), stencilBuffer_(NULL), width_(0), height_(0),
      attachments_(false), validSize_(false) {
  for (GLuint i = 0; i < kMaxColorAttachments; i++) {
    colorBuffers_[i] = NULL;
  }
}

Framebuffer::~Framebuffer() {}

void Framebuffer::AttachColorBuffer(GLuint index, Renderbuffer *rb) {
  assert(index < kMaxColorAttachments && "color buffer index out of range");
  colorBuffers_[index] = rb;
  UpdateAttachmentStatus();
}

void Framebuffer::AttachDepthBuffer(Renderbuffer *rb) {
  depthBuffer_ = rb;
  UpdateAttachmentStatus();
}

void Framebuffer::AttachStencilBuffer(Renderbuffer *rb) {
  stencilBuffer_ = rb;
  UpdateAttachmentStatus();
}

void Framebuffer::UpdateAttachmentStatus() {
  // initially assume attachments, but with an invalid size
  attachments_ = true;
  validSize_ = false;

  // iterate through all attached color buffers first
  const Renderbuffer *first = NULL, *rb;
  for (GLuint i = 0; i < kMaxColorAttachments; i++) {
    rb = GetColorBuffer(i);
    if (rb != NULL) {
      // store as first buffer found if not set yet
      if (first == NULL) {
        first = rb;
      }

      // fail if size of this renderbuffer is different from the first
      // renderbuffer found
      if (rb->IsSameSize(first) == false) {
        return;
      }
    }
  }

  // check depth buffer next after color buffers
  rb = GetDepthBuffer();
  if (rb != NULL) {
    // store as first buffer found if not set yet
    if (first == NULL) {
      first = rb;
    }

    // fail if the size is different
    if (rb->IsSameSize(first) == false) {
      return;
    }
  }

  // then finally check the stencil buffer
  rb = GetStencilBuffer();
  if (rb != NULL) {
    // store as first buffer found if not set yet
    if (first == NULL) {
      first = rb;
    }

    // fail if the size is different
    if (rb->IsSameSize(first) == false) {
      return;
    }
  }

  // after iterating through all attachment pointers, if the first renderbuffer
  // is null, then we didn't actually find any attachments
  if (first == NULL) {
    attachments_ = false;
  } else {
    // there is at least one attachment, and the sizes are valid across all
    // attachments, so we have a valid size
    width_ = first->GetWidth();
    height_ = first->GetHeight();
    validSize_ = true;
  }
}
