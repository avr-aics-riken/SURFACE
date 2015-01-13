/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "GLES2/gl2.h"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include "gles_common.h"
#include "../render/matrix.h"

using namespace lsgl;

namespace {

// Build lsgl predefined vertex uniform variables
void BuildPredefinedLSGLVertexUniforms(
    std::vector<Uniform> &vertUniforms,
    std::vector<UniformLocation> &vertUniformLocations, int indexOffset) {

  // index may exceeds kMaxVertexUniformVectors, but would be OK.
  UniformLocation locWorld("lsgl_World", /*elem=*/0, indexOffset);
  UniformLocation locView("lsgl_View", /*elem=*/0, indexOffset + 1);
  UniformLocation locProj("lsgl_Proj", /*elem=*/0, indexOffset + 2);
  UniformLocation locPointSize("lsgl_PointSize", /*elem=*/0, indexOffset + 3);

  // SHOULD not break this ordering
  vertUniformLocations.push_back(locWorld);
  vertUniformLocations.push_back(locView);
  vertUniformLocations.push_back(locProj);
  vertUniformLocations.push_back(locPointSize);

  Uniform world(GL_FLOAT_MAT4, "lsgl_World", 1);
  Uniform view(GL_FLOAT_MAT4, "lsgl_View", 1);
  Uniform proj(GL_FLOAT_MAT4, "lsgl_Proj", 1);
  Uniform psize(GL_FLOAT, "lsgl_PointSize", 1);

  // fill data with identity matrix.
  float ident[4][4];
  Matrixf::Identity(ident);

  int sz = sizeof(float) * 4 * 4; // 4x4 float matrix.
  world.data.resize(sz);
  memcpy(&world.data.at(0), ident, sz);
  view.data.resize(sz);
  memcpy(&view.data.at(0), ident, sz);
  proj.data.resize(sz);
  memcpy(&proj.data.at(0), ident, sz);

  // fill with 1.0
  float initPointSize = 1.0f;
  psize.data.resize(sizeof(float));
  memcpy(&psize.data.at(0), &initPointSize, sizeof(float));

  // SHOULD not break this ordering
  vertUniforms.push_back(world);
  vertUniforms.push_back(view);
  vertUniforms.push_back(proj);
  vertUniforms.push_back(psize);
}

} // namespace

//
// Program
//
Program::Program() : linked_(false) {}

Program::~Program() {}

bool Program::AttachShader(Shader *shd) {
  // fail if shader is already attached
  if (IsAttached(shd) == true) {
    return false;
  }

  // add shader to appropriate list based on it's type
  switch (shd->GetType()) {
  case GL_VERTEX_SHADER: {
    VertexShader *vtx = static_cast<VertexShader *>(shd);
    vertexShaders_.push_back(vtx);
  } break;

  // case GL_GEOMETRY_SHADER:
  //    {
  //        GeometryShader *geo = static_cast<GeometryShader *>(shd);
  //        geometryShaders_.push_back(geo);
  //    }
  //    break;

  case GL_FRAGMENT_SHADER: {
    FragmentShader *frag = static_cast<FragmentShader *>(shd);
    fragmentShaders_.push_back(frag);
  } break;
  }

  return true;
}

bool Program::DetachShader(Shader *shd) {
  int t;

  // fail if shader is already detached
  if (IsAttached(shd) == false) {
    return false;
  }

  // remove shader from appropriate list based on it's type
  switch (shd->GetType()) {
  case GL_VERTEX_SHADER:
    for (t = 0; t < vertexShaders_.size(); t++) {
      if (vertexShaders_[t] == shd) {
        vertexShaders_.erase(vertexShaders_.begin() + t);
        return true;
      }
    }
    break;

  // case GL_GEOMETRY_SHADER:
  //    for (t = 0; t < geometryShaders_.size(); t++) {
  //        if (geometryShaders_[t] == shd) {
  //            geometryShaders_.erase(geometryShaders_.begin() + t);
  //            return true;
  //        }
  //    }
  //    break;

  case GL_FRAGMENT_SHADER:
    for (t = 0; t < fragmentShaders_.size(); t++) {
      if (fragmentShaders_[t] == shd) {
        fragmentShaders_.erase(fragmentShaders_.begin() + t);
        return true;
      }
    }
    break;
  }

  return false;
}

bool Program::Link() {
  // reset linked state
  linked_ = false;

  // for now, just ensure that there is at least one attached shader, and that
  // all attached shaders compiled properly
  if (GetAttachedCount() < 1) {
    return false;
  }

  for (size_t t = 0; t < vertexShaders_.size(); t++) {
    if (vertexShaders_[t]->IsCompiled() == false) {
      return false;
    }
  }

  for (size_t t = 0; t < fragmentShaders_.size(); t++) {
    if (fragmentShaders_[t]->IsCompiled() == false) {
      fprintf(stderr, "FragmentShader::IsCompiled = false\n");
      return false;
    }
  }

  // Do post link operation.
  //   - Bind uniform location.
  // @todo { support vertex shader }

  uniforms_.clear();
  uniformLocations_.clear();
  varyings_.clear();
  varyingLocations_.clear();

  assert(fragmentShaders_.size() == 1);
  for (size_t t = 0; t < fragmentShaders_.size(); t++) {

    FragmentShader *shader = fragmentShaders_[t];

    std::vector<Uniform> fragUniforms;
    std::vector<UniformLocation> fragUniformLocations;
    shader->BuildUniformInfo(fragUniforms, fragUniformLocations);

    // Append

    uniforms_.insert(uniforms_.end(), fragUniforms.begin(), fragUniforms.end());
    uniformLocations_.insert(uniformLocations_.end(),
                             fragUniformLocations.begin(),
                             fragUniformLocations.end());

    // @todo { Full support for vertex attribute }

    std::vector<Varying> fragVaryings;
    std::vector<VaryingLocation> fragVaryingLocations;
    shader->BuildVaryingInfo(fragVaryings, fragVaryingLocations);

    // Append
    varyings_.insert(varyings_.end(), fragVaryings.begin(), fragVaryings.end());
    varyingLocations_.insert(varyingLocations_.end(),
                             fragVaryingLocations.begin(),
                             fragVaryingLocations.end());
  }

  // DBG
  printf("Linked: # of uniforms = %ld\n", uniforms_.size());
  printf("Linked: # of varyings = %ld\n", varyings_.size());

  // Add predefined uniform
  std::vector<Uniform> vertUniforms;
  std::vector<UniformLocation> vertUniformLocations;
  BuildPredefinedLSGLVertexUniforms(vertUniforms, vertUniformLocations,
                                    uniforms_.size());

  uniforms_.insert(uniforms_.end(), vertUniforms.begin(), vertUniforms.end());
  uniformLocations_.insert(uniformLocations_.end(),
                           vertUniformLocations.begin(),
                           vertUniformLocations.end());

  linked_ = true;
  return true;
}

bool Program::PrepareEval(FragmentState &fragmentState,
                          ShadingState &shadingState,
                          const std::vector<VertexAttribute> &vertexAttributes,
                          const Context &ctx) const {
  assert(GetFragmentShaderCount() == 1);

  // @fixme { Support vertex shader. }

  // Prepare fragment shader
  {
    bool ret = GetFragmentShader(0)->PrepareEval(
        fragmentState, shadingState, uniforms_, varyings_, varyingLocations_,
        vertexAttributes, ctx);
    assert(ret);
  }

  return true;
}

bool Program::IsAttached(Shader *shd) {
  int t;

  // check all shader lists to see if this shader is already attached
  for (t = 0; t < vertexShaders_.size(); t++) {
    if (vertexShaders_[t] == shd) {
      return true;
    }
  }

  // for (t = 0; t < geometryShaders_.size(); t++) {
  //    if (geometryShaders_[t] == shd) {
  //        return true;
  //    }
  //}

  for (t = 0; t < fragmentShaders_.size(); t++) {
    if (fragmentShaders_[t] == shd) {
      return true;
    }
  }

  return false;
}

GLint Program::GetUniformLocation(const std::string &name) const {
  unsigned int subscript = 0; // array num
  std::string var = name;

  // Strip any trailing array operator and retrieve the subscript
  size_t open = var.find_last_of('[');
  size_t close = var.find_last_of(']');
  if (open != std::string::npos && close == var.length() - 1) {
    subscript = atoi(var.substr(open + 1).c_str());
    var.erase(open);
  }

  // Simple linear search to find a location by name.
  unsigned int n = uniformLocations_.size();
  // printf("dbg: n loc = %d\n", n);
  for (unsigned int location = 0; location < n; location++) {
    // printf("dbg: %d| name = %s, elem = %d, subscript = %d\n", location,
    //       uniformLocations_[location].name.c_str(),
    //       uniformLocations_[location].element, subscript);

    if (uniformLocations_[location].name == name &&
        uniformLocations_[location].element == subscript) {
      return location;
    }
  }

  return -1;
}

Uniform const *Program::GetUniform(GLint location) const {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return NULL;
  }

  // Look up uniform variable info
  Uniform const *uniform = &uniforms_[uniformLocations_[location].index];

  return uniform;
}

bool Program::SetUniform1iv(GLint location, GLsizei count, const GLint *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_INT || uniform->type == GL_SAMPLER_2D ||
      uniform->type == GL_SAMPLER_3D || uniform->type == GL_SAMPLER_CUBE) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));

    assert(uniform->data.size() == sizeof(GLint) * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLint),
           v, sizeof(GLint) * n);

  } else if (uniform->type == GL_BOOL) {
    assert(0 && "Bool type is todo.");
  } else {
    // Unsupported type.
    assert(0 && "Unsupported type.");
    return false;
  }

  return true;
}

bool Program::SetUniform1fv(GLint location, GLsizei count, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf("MATCH:\n");
    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));
    // printf(" data = %f\n", v[0]);

    assert(uniform->data.size() == sizeof(GLfloat) * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat),
           v, sizeof(GLfloat) * n);

  } else {
    // Type mismatch.
    return false;
  }

  return true;
}

bool Program::SetUniform2fv(GLint location, GLsizei count, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    printf("[LSGL] ERR: location out of range. loc = %d\n", location);
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT_VEC2) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      printf("[LSGL] ERR: Type signature mismatch\n");
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));

    assert(uniform->data.size() == sizeof(GLfloat) * 2 * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat) * 2,
           v, sizeof(GLfloat) * 2 * n);

  } else {
    // Type mismatch.
    printf("[LSGL] ERR: Type signature mismatch. ty = %d\n", uniform->type);
    return false;
  }

  return true;
}

bool Program::SetUniform3fv(GLint location, GLsizei count, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT_VEC3) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));

    assert(uniform->data.size() == sizeof(GLfloat) * 3 * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat) * 3,
           v, sizeof(GLfloat) * 3 * n);

  } else {
    // Type mismatch.
    return false;
  }

  return true;
}

bool Program::SetUniform4fv(GLint location, GLsizei count, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT_VEC4) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));

    assert(uniform->data.size() == sizeof(GLfloat) * 4 * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat) * 4,
           v, sizeof(GLfloat) * 4 * n);

  } else {
    // Type mismatch.
    return false;
  }

  return true;
}

bool Program::SetUniformMatrix2fv(GLint location, GLsizei count,
                                  GLboolean transpose, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT_MAT2) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf("MATCH:\n");
    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));
    // printf(" data = %f\n", v[0]);

    assert(uniform->data.size() == sizeof(GLfloat) * 2 * 2 * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat) * 2 * 2,
           v, sizeof(GLfloat) * 2 * 2 * n);

  } else {
    // Type mismatch.
    return false;
  }

  return true;
}

bool Program::SetUniformMatrix3fv(GLint location, GLsizei count,
                                  GLboolean transpose, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT_MAT3) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf("MATCH:\n");
    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));
    // printf(" data = %f\n", v[0]);

    assert(uniform->data.size() == sizeof(GLfloat) * 3 * 3 * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat) * 3 * 3,
           v, sizeof(GLfloat) * 3 * 3 * n);

  } else {
    // Type mismatch.
    return false;
  }

  return true;
}

bool Program::SetUniformMatrix4fv(GLint location, GLsizei count,
                                  GLboolean transpose, const GLfloat *v) {
  if ((location < 0) || (location >= uniforms_.size())) {
    // Out of range.
    return false;
  }

  // Look up uniform variable info
  Uniform *uniform = &uniforms_[uniformLocations_[location].index];

  if (uniform->type == GL_FLOAT_MAT4) {
    int arraySize = uniform->arraySize;
    if ((arraySize == 1) && (count > 1)) {
      // Signature mismatch.
      return false;
    }

    int n =
        std::min(arraySize - (int)uniformLocations_[location].element, count);

    // printf("MATCH:\n");
    // printf(" n = %d\n", n);
    // printf(" elem = %d\n", uniformLocations_[location].element);
    // printf(" data.size() = %d\n", (int)(uniform->data.size()));
    // printf(" data = %f\n", v[0]);

    assert(uniform->data.size() == sizeof(GLfloat) * 4 * 4 * n);

    memcpy(&uniform->data.at(0) +
               uniformLocations_[location].element * sizeof(GLfloat) * 4 * 4,
           v, sizeof(GLfloat) * 4 * 4 * n);

  } else {
    // Type mismatch.
    return false;
  }

  return true;
}

GLint Program::GetVaryingLocation(const std::string &name) const {
  unsigned int subscript = 0; // array num
  std::string var = name;

  // Strip any trailing array operator and retrieve the subscript
  size_t open = var.find_last_of('[');
  size_t close = var.find_last_of(']');
  if (open != std::string::npos && close == var.length() - 1) {
    subscript = atoi(var.substr(open + 1).c_str());
    var.erase(open);
  }

  // Simple linear search to find a location by name.
  unsigned int n = varyingLocations_.size();
  // printf("n loc = %d\n", n);
  for (unsigned int location = 0; location < n; location++) {
    // printf("%d| name = %s, elem = %d, subscript = %d\n", location,
    //       varyingLocations_[location].name.c_str(),
    //       varyingLocations_[location].element, subscript);

    if (varyingLocations_[location].name == name &&
        varyingLocations_[location].element == subscript) {
      return location;
    }
  }

  return -1;
}

Varying const *Program::GetVarying(GLint location) const {
  if ((location < 0) || (location > varyings_.size())) {
    // Out of range.
    assert(0);
    return NULL;
  }

  // Look up varying variable info
  Varying const *varying = &varyings_[varyingLocations_[location].index];

  return varying;
}
