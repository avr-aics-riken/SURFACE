/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include "GLES2/gl2.h"

#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#include <cassert>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cerrno>

#if defined(_WIN32)
#include <io.h> // _mktemp
#endif

#include "gles_context.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

using namespace lsgl;

namespace {

/// Generate unique filename. User must free memory of returned pointer when it
/// become unused.
char *GenerateUniqueFilename() {

#if defined(_WIN32)
  char basename[] = "lsgl_shader_XXXXXX";
#else

#ifdef LSGL_ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char basename[1024];
  sprintf(basename, "lsgl_shader_%08d_XXXXXX", rank);
#else
  char basename[] = "lsgl_shader_XXXXXX";
#endif
#endif

#if defined(_WIN32)

  // NOTE: Windows can't exec DLL in Temp file.
  // char *name = strdup("lsgl_shader.dll"); //_tempnam(NULL, prefix);

  assert(strlen(basename) < 1023);

  char buf[1024];
  strncpy(buf, basename, strlen(basename));
  buf[strlen(basename)] = '\0';

  size_t len = strlen(buf) + 1;
  int err =
      _mktemp_s(buf, len); // prefix will be overwritten with actual filename
  if (err != 0) {
    fprintf(stderr, "Failed to create unique filename.\n");
    return NULL;
  }

  printf("DBG: Unique name: %s", buf);

  char *name = strdup(buf);

#else
  // The string returned by tempnam() is allocated using malloc(3) and hence
  // should be freed by free(3).
  // char *name = tempnam(NULL, prefix);
  // if (!name) {
  //  if (errno == EINVAL) {
  //    fprintf(stderr, "Bad template parameter.\n");
  //    return NULL;
  //  } else if (errno == EEXIST) {
  //    fprintf(stderr, "Out of unique filename.\n");
  //    return NULL;
  //  } else {
  //    fprintf(stderr, "Failed to create unique name.\n");
  //    return NULL;
  //  }
  //}

  int fd = mkstemp(basename); // prefix will be overwritten with actual filename
  if (fd == -1) {
    fprintf(stderr, "Failed to create unique filename.\n");
    return NULL;
  }
  close(fd);
  int ret = unlink(basename);
  if (ret == -1) {
    fprintf(stderr, "[LSGL] Failed to delete file: %s\n", basename);
  }

  char *name = strdup(basename);

#endif

  return name;
}

} // namespace

//
// --
//

void Context::glAttachShader(GLuint program, GLuint shader) {
  TRACE_EVENT("(GLuint program = %d, GLuint shader = %d)", program, shader);

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shader);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // attach the shader to the program
  if (prg->AttachShader(shd) == false) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // Set attach frag to the shader.
  shd->SetAttached();
}

void Context::glGetAttachedShaders(GLuint program, GLsizei maxcount,
                                   GLsizei *count, GLuint *shaders) {
  TRACE_EVENT("(GLuint program = %d, GLsizei maxcount = %d, GLsizei* count = "
              "%p, GLuint* shaders = %p)",
              program, maxcount, count, shaders);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glDetachShader(GLuint program, GLuint shader) {
  TRACE_EVENT("(GLuint program = %d, GLuint shader = %d)", program, shader);

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shader);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // detach the shader to the program
  if (prg->DetachShader(shd) == false) {
    return SetGLError(GL_INVALID_OPERATION);
  }
}

void Context::glCompileShader(GLuint shader) {
  TRACE_EVENT("(GLuint shader = %d)", shader);

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shader);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // compile the shader
  shd->Compile();
}

GLuint Context::glCreateProgram(void) {
  TRACE_EVENT("()");
  return resourceManager_.CreateProgram();
}

GLuint Context::glCreateShader(GLenum type) {
  TRACE_EVENT("(GLenum type = %d)", type);

  // check for valid type
  if ((type != GL_VERTEX_SHADER) && (type != GL_FRAGMENT_SHADER)) {
    SetGLError(GL_INVALID_ENUM);
    return 0;
  }

  return resourceManager_.CreateShader(type);
}

void Context::glDeleteProgram(GLuint program) {
  TRACE_EVENT("(GLuint program = %d)", program);
  return resourceManager_.DeleteProgram(program);
}

void Context::glDeleteShader(GLuint shader) {
  TRACE_EVENT("(GLuint shader = %d)", shader);

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shader);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  if (shd->IsAttached()) {
    // The spec says,
    //   If a shader object to be deleted is attached to a program object,
    //   it will be flagged for deletion, but it will not be deleted
    //   until it is no longer attached to any program object, for any
    //   rendering context
    //
    // http://www.opengl.org/sdk/docs/man/xhtml/glDeleteShader.xml

    // @fixme { Mark deletion to this shader object. }

    return;

  } else {

    // Release internal resource.
    shd->Release();

    return resourceManager_.DeleteShader(shader);
  }
}

GLboolean Context::glIsProgram(GLuint program) {
  TRACE_EVENT("(GLuint program = %d)", program);
  return resourceManager_.IsValidProgram(program);
}

GLboolean Context::glIsShader(GLuint shader) {
  TRACE_EVENT("(GLuint shader = %d)", shader);
  return resourceManager_.IsValidShader(shader);
}

void Context::glLinkProgram(GLuint program) {
  TRACE_EVENT("(GLuint program = %d)", program);

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // link the program
  prg->Link();
}

void Context::glShaderBinary(GLsizei n, const GLuint *shaders,
                             GLenum binaryformat, const GLvoid *binary,
                             GLsizei length) {
  TRACE_EVENT(
      "(GLsizei n = %d, const GLuint* shaders = %p, GLenum binaryformat = %d, "
      "const GLvoid* binary = %p, GLsizei length = %d)",
      n, shaders, binaryformat, binary, length);

  if (shaders == NULL) {
    return SetGLError(GL_INVALID_VALUE);
  }

  if (n <= 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  if (length <= 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // @fixme { Only handle first element of shaders }
  assert(n == 1);

  // binaryformat is just ignored because the spec says,
  //   OpenGL ES defines no specific binary formats, but does provide
  //   a mechanism to obtain symbolic constants for such formats provided
  //   by extensions

  // Dump shader binary to a file, then read it as DLL.
  // Since content of 'binary' is actually DLL in LSGL, this save-to-file
  // operation is a bit redundant.
  char *filename = GenerateUniqueFilename();
  assert(filename);

  printf("[DBG] tempfilename = %s\n", filename);

  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    free(filename);
    fprintf(stderr, "[LSGL] Failed to create file to save shader binary\n");
    return SetGLError(GL_INVALID_OPERATION);
  }

  size_t count = fwrite(binary, length, 1, fp);
  assert(count == 1);

  fclose(fp);

#if defined(_WIN32)
  // On windows, it seems dll must have .dll extension to get LoadLibrary
  // success.
  // Save .dll file also.
  // If we write .dll file directy, GenerateUniqueFilename() will generate same
  // filename...
  // @todo { remove this redundant file write. }

  {
    std::string dllname = std::string(filename) + ".dll";

    FILE *fp = fopen(dllname.c_str(), "wb");
    assert(fp);

    size_t count = fwrite(binary, length, 1, fp);
    assert(count == 1);

    fclose(fp);
  }
#endif

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shaders[0]);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

#if defined(_WIN32)
  std::string filepath = std::string(filename) + ".dll";
#else
  std::string filepath = std::string(filename);
#endif
  // free(filename);

  bool ret = shd->LoadShaderBinary(filepath);
  if (!ret) {
    return SetGLError(GL_INVALID_OPERATION);
  }
}

void Context::glShaderSource(GLuint shader, GLsizei count,
                             const GLchar **string, const GLint *length) {
  TRACE_EVENT("(GLuint shader = %d, GLsizei count = %d, const GLchar** string "
              "= %p, const GLint* length = %p)",
              shader, count, string, length);

  // fail if count is less than zero
  if (count < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shader);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_VALUE);
    // Potentially GL_INVALUE_OPERATION.
    // The spec says: GL_INVALID_OPERATION is generated if shader is not a
    // shader object.
  }

  // set the source for this shader
  shd->Source(count, string, length);
}

void Context::glUseProgram(GLuint program) {
  TRACE_EVENT("(GLuint program = %d)", program);

  // GLESv2 spec says,
  // If UseProgram is called with program set to zero, then the current
  // rendering state refers to an invalid program object, and the results of
  // vertex and fragment shader execution due to any DrawArrays or DrawElements
  // commands are undefined. However, this is not an error.
  if (program == 0) {
    // Undefined, but OK.
    state_.currentProgram = 0;
    return;
  }

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // ensure program linked properly
  if (prg->IsLinked() == false) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  state_.currentProgram = program;
}

void Context::glValidateProgram(GLuint program) {
  TRACE_EVENT("(GLuint program = %d)", program);
  assert(0 && "glValidateProgram is not yet implemented.");
}

void Context::glReleaseShaderCompiler(void) {
  TRACE_EVENT("()");
  // ignore for now
}

void Context::glGetShaderiv(GLuint shader, GLenum pname, GLint *params) {
  TRACE_EVENT("(GLuint shader = %d, GLenum name = %d, GLuint* params = %p)",
              shader, pname, params);

  // get shader pointer from handle
  Shader *shd = resourceManager_.GetShader(shader);
  if (shd == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  switch (pname) {
  case GL_SHADER_TYPE:
    *params = shd->GetType();
    break;

  case GL_DELETE_STATUS:
    // TODO
    *params = GL_FALSE;
    break;

  case GL_COMPILE_STATUS:
    *params = shd->IsCompiled();
    break;

  case GL_INFO_LOG_LENGTH:
    // TODO
    *params = 0;
    break;

  case GL_SHADER_SOURCE_LENGTH:
    // add one to include the null terminator
    *params = shd->GetSource().size() + 1;
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

void Context::glGetShaderPrecisionFormat(GLenum shadertype,
                                         GLenum precisiontype, GLint *range,
                                         GLint *precision) {
  TRACE_EVENT("(GLenum shadertype = %d, GLenum precisiontype = %d, GLint* "
              "range = %p, GLint* precision = %p)",
              shadertype, precisiontype, range, precision);

  assert(0 && "TODO");
  SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetProgramiv(GLuint program, GLenum pname, GLint *params) {
  TRACE_EVENT("(GLuint program = %d, GLenum name = %d, GLuint* params = %p)",
              program, pname, params);

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  switch (pname) {
  case GL_DELETE_STATUS:
    // TODO
    *params = GL_FALSE;
    break;

  case GL_LINK_STATUS:
    *params = prg->IsLinked();
    break;

  case GL_VALIDATE_STATUS:
    // TODO
    *params = 0;
    break;

  case GL_INFO_LOG_LENGTH:
    // TODO
    *params = 0;
    break;

  case GL_ATTACHED_SHADERS:
    *params = prg->GetAttachedCount();
    break;

  case GL_ACTIVE_ATTRIBUTES:
    // TODO
    *params = 0;
    break;

  case GL_ACTIVE_ATTRIBUTE_MAX_LENGTH:
    // TODO
    *params = 0;
    break;

  case GL_ACTIVE_UNIFORMS:
    // TODO
    *params = 0;
    break;

  case GL_ACTIVE_UNIFORM_MAX_LENGTH:
    // TODO
    *params = 0;
    break;

  default:
    return SetGLError(GL_INVALID_ENUM);
  }
}

int Context::glGetAttribLocation(GLuint program, const GLchar *name) {
  TRACE_EVENT("(GLuint program = %d, const GLchar* name = %s)", program, name);

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    SetGLError(GL_INVALID_OPERATION);
    return -1;
  }

  // ensure program has been successfully linked
  if (prg->IsLinked() == false) {
    SetGLError(GL_INVALID_OPERATION);
    return -1;
  }

  // @fixme { used fixed attribute locations for now }
  if (!strcmp(name, "position")) {
    return kVtxAttrPosition;
  } else if (!strcmp(name, "normal")) {
    return kVtxAttrNormal;
  } else if (!strcmp(name, "texCoord")) {
    return kVtxAttrTexCoord;
  } else if (!strcmp(name, "material")) {
    return kVtxAttrMaterial;
  }

  // Fall back to varying varible in the shader.
  // @fixme { Fall back to vertex attribute. }
  std::string s = std::string(name);
  GLint idx = prg->GetVaryingLocation(s);

  if (idx == -1) {
    // unknown attribute
    fprintf(stderr, "[LSGL] Cannot find a location for %s.\n", name);
    return -1;
  }

  return static_cast<GLuint>(idx);
}

void Context::glUniform1f(GLint location, GLfloat x) {
  TRACE_EVENT("(GLuint location = %d, GLfloat x = %f)", location, x);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  int n = 1;
  prg->SetUniform1fv(location, n, &x);
}

void Context::glUniform1fv(GLint location, GLsizei count, const GLfloat *v) {
  TRACE_EVENT(
      "(GLint location = %d, GLsiszei count = %d, const GLfloat* v = %p)",
      location, count, v);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform1i(GLint location, GLint x) {
  TRACE_EVENT("(GLint location = %d, GLint x = %d)", location, x);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  Uniform const *uniform = prg->GetUniform(location);
  if (uniform == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  GLint val;

  // If the uniform at the location is sampler type, do slot -> handle remap
  if (uniform->type == GL_SAMPLER_2D) {
    val = state_.texture2D[x];
  } else if (uniform->type == GL_SAMPLER_3D) {
    val = state_.texture3D[x];
  } else if (uniform->type == GL_SAMPLER_CUBE) {
    val = state_.textureCubeMap[x];
  } else {
    val = x;
  }

  int n = 1;
  prg->SetUniform1iv(location, n, &val);
}

void Context::glUniform1iv(GLint location, GLsizei count, const GLint *v) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, const GLint* v = %p)",
              location, count, v);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform2f(GLint location, GLfloat x, GLfloat y) {
  TRACE_EVENT("(GLuint location = %d, GLfloat x = %f, GLfloat y = %f)",
              location, x, y);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  GLfloat val[2] = {x, y};
  int n = 1;

  bool ret = prg->SetUniform2fv(location, n, val);
  assert(ret);
}

void Context::glUniform2fv(GLint location, GLsizei count, const GLfloat *v) {
  TRACE_EVENT(
      "(GLuint location = %d, GLsizei count = %d, const GLfloat* v = %p)",
      location, count, v);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  bool ret = prg->SetUniform2fv(location, count, v);
  assert(ret);

  return;
}

void Context::glUniform2i(GLint location, GLint x, GLint y) {
  TRACE_EVENT("(GLint location = %d, GLint x = %d, GLint y = %d)", location, x,
              y);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform2iv(GLint location, GLsizei count, const GLint *v) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, const GLint* v = %p)",
              location, count, v);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform3f(GLint location, GLfloat x, GLfloat y, GLfloat z) {
  TRACE_EVENT(
      "(GLuint location = %d, GLfloat x = %f, GLfloat y = %f, GLfloat = %f)",
      location, x, y, z);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  GLfloat val[3] = {x, y, z};
  int n = 1;

  bool ret = prg->SetUniform3fv(location, n, val);
  assert(ret);
}

void Context::glUniform3fv(GLint location, GLsizei count, const GLfloat *v) {
  TRACE_EVENT(
      "(GLuint location = %d, GLsizei count = %d, const GLfloat* v = %p)",
      location, count, v);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  bool ret = prg->SetUniform3fv(location, count, v);
  if (!ret) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  return;
}

void Context::glUniform3i(GLint location, GLint x, GLint y, GLint z) {
  TRACE_EVENT("(GLint location = %d, GLint x = %d, GLint y = %d, GLint z = %d)",
              location, x, y, z);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform3iv(GLint location, GLsizei count, const GLint *v) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, const GLint* v = %p)",
              location, count, v);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform4f(GLint location, GLfloat x, GLfloat y, GLfloat z,
                          GLfloat w) {
  TRACE_EVENT("(GLuint location = %d, GLfloat x = %f, GLfloat y = %f, GLfloat "
              "= %f, GLfloat w = %f)",
              location, x, y, z, w);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  GLfloat val[4] = {x, y, z, w};
  int n = 1;
  bool ret = prg->SetUniform4fv(location, n, val);

  return;
}

void Context::glUniform4fv(GLint location, GLsizei count, const GLfloat *v) {
  TRACE_EVENT(
      "(GLuint location = %d, GLsizei count = %d, const GLfloat* v = %p)",
      location, count, v);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  bool ret = prg->SetUniform4fv(location, count, v);
  assert(ret);
}

void Context::glUniform4i(GLint location, GLint x, GLint y, GLint z, GLint w) {
  TRACE_EVENT("(GLuint location = %d, GLint x = %d, GLint y = %d, GLint z = "
              "%d, GLint w = %d)",
              location, x, y, z, w);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniform4iv(GLint location, GLsizei count, const GLint *v) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, const GLint* v = %p)",
              location, count, v);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glUniformMatrix2fv(GLint location, GLsizei count,
                                 GLboolean transpose, const GLfloat *value) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, GLboolean transpose "
              "= %d, const GLfloat* value = %p)",
              location, count, transpose, value);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  bool ret = prg->SetUniformMatrix2fv(location, count, transpose, value);
  assert(ret);
}

void Context::glUniformMatrix3fv(GLint location, GLsizei count,
                                 GLboolean transpose, const GLfloat *value) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, GLboolean transpose "
              "= %d, const GLfloat* value = %p)",
              location, count, transpose, value);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  bool ret = prg->SetUniformMatrix3fv(location, count, transpose, value);
  assert(ret);
}

void Context::glUniformMatrix4fv(GLint location, GLsizei count,
                                 GLboolean transpose, const GLfloat *value) {
  TRACE_EVENT("(GLuint location = %d, GLsizei count = %d, GLboolean transpose "
              "= %d, const GLfloat* value = %p)",
              location, count, transpose, value);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  bool ret = prg->SetUniformMatrix4fv(location, count, transpose, value);
  assert(ret);
}

void Context::glGetShaderInfoLog(GLuint shader, GLsizei bufsize,
                                 GLsizei *length, GLchar *infolog) {
  TRACE_EVENT("(GLuint shader = %d, GLsizei bufsize = %d, GLsizei* length = "
              "%p, GLchar* infolog = %p)",
              shader, bufsize, length, infolog);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetProgramInfoLog(GLuint program, GLsizei bufsize,
                                  GLsizei *length, GLchar *infolog) {
  TRACE_EVENT("(GLuint program = %d, GLsizei bufsize = %d, GLsizei* length = "
              "%p, GLchar* infolog = %p)",
              program, bufsize, length, infolog);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetVertexAttribfv(GLuint index, GLenum pname, GLfloat *params) {
  TRACE_EVENT("(GLuint index = %d, GLenum pname = %d, GLfloat* params = %p)",
              index, pname, params);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetVertexAttribiv(GLuint index, GLenum pname, GLint *params) {
  TRACE_EVENT("(GLuint index = %d, GLenum pname = %d, GLint* params = %p)",
              index, pname, params);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetVertexAttribPointerv(GLuint index, GLenum pname,
                                        GLvoid **pointer) {
  TRACE_EVENT("(GLuint index = %d, GLenum pname = %d, GLvoid** pointer = %p)",
              index, pname, pointer);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetUniformfv(GLuint program, GLint location, GLfloat *params) {
  TRACE_EVENT(
      "(GLuint program = %d, GLint location = %d, GLfloat* params = %p)",
      program, location, params);

  if (location < 0) {
    return SetGLError(GL_INVALID_VALUE);
  }

  if (params == NULL) {
    return SetGLError(GL_INVALID_VALUE);
  }

  // get program pointer from current program handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  Uniform const *uniform = prg->GetUniform(location);
  if (uniform == NULL) {
    return SetGLError(GL_INVALID_VALUE);
  }

  switch (uniform->type) {

  case GL_FLOAT: {
    float val = *(reinterpret_cast<float const *>(&uniform->data.at(0)));
    params[0] = val;
  } break;
  case GL_FLOAT_VEC2: {
    float const *ptr = reinterpret_cast<float const *>(&uniform->data.at(0));
    memcpy(params, ptr, sizeof(float) * 2);
  } break;
  case GL_FLOAT_VEC3: {
    float const *ptr = reinterpret_cast<float const *>(&uniform->data.at(0));
    memcpy(params, ptr, sizeof(float) * 3);
  } break;
  case GL_FLOAT_VEC4: {
    float const *ptr = reinterpret_cast<float const *>(&uniform->data.at(0));
    memcpy(params, ptr, sizeof(float) * 4);
  } break;
  case GL_FLOAT_MAT4: {
    float const *ptr = reinterpret_cast<float const *>(&uniform->data.at(0));
    memcpy(params, ptr, sizeof(float) * 4 * 4);
  } break;
  default:
    // Unknown or unsupported data type.
    assert(0 && "TODO");
    break;
  }

  return;
}

void Context::glGetUniformiv(GLuint program, GLint location, GLint *params) {
  TRACE_EVENT("(GLuint program = %d, GLint location = %d, GLint* params = %p)",
              program, location, params);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetActiveAttrib(GLuint program, GLuint index, GLsizei bufsize,
                                GLsizei *length, GLint *size, GLenum *type,
                                GLchar *name) {
  TRACE_EVENT(
      "(GLuint program = %d, GLuint index = %d, GLsizei bufsize = %d, GLsizei* "
      "length = %p, GLint* size = %p, GLenum* type = %p, GLchar* name = %p)",
      program, index, bufsize, length, size, type, name);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetActiveUniform(GLuint program, GLuint index, GLsizei bufsize,
                                 GLsizei *length, GLint *size, GLenum *type,
                                 GLchar *name) {
  TRACE_EVENT(
      "(GLuint program = %d, GLuint index = %d, GLsizei bufsize = %d, GLsizei* "
      "length = %p, GLint* size = %p, GLenum* type = %p, GLchar* name = %p)",
      program, index, bufsize, length, size, type, name);
  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glBindAttribLocation(GLuint program, GLuint index,
                                   const GLchar *name) {
  TRACE_EVENT(
      "(GLuint program = %d, GLuint index = %d, const GLchar* name = %p)",
      program, index, name);

  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

void Context::glGetShaderSource(GLuint shader, GLsizei bufsize, GLsizei *length,
                                GLchar *source) {
  TRACE_EVENT("(GLuint shader = %d, GLsizei bufsize = %d, GLsizei* length = "
              "%p, GLchar* source = %p)",
              shader, bufsize, length, source);

  assert(0 && "TODO");
  return SetGLError(GL_INVALID_OPERATION);
}

int Context::glGetUniformLocation(GLuint program, const GLchar *name) {
  TRACE_EVENT("(GLuint program = %d, const GLchar* name = %s)", program, name);

  // return -1 if builtin variables
  if (strstr(name, "gl_") == name) {
    return -1;
  }

  // get program pointer from handle
  Program *prg = resourceManager_.GetProgram(program);
  if (prg == NULL) {
    SetGLError(GL_INVALID_OPERATION);
    return -1;
  }

  std::string str(name);
  return prg->GetUniformLocation(str);
}
