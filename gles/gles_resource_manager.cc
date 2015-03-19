/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "gles_resource_manager.h"

using namespace lsgl;

GLuint ResourceManager::CreateBuffer() {
  GLuint handle = bufferHandleAllocator_.Allocate();
  bufferMap_[handle] = new Buffer;
  return handle;
}

GLuint ResourceManager::CreateTexture() {
  GLuint handle = textureHandleAllocator_.Allocate();
  textureMap_[handle] = new Texture;
  return handle;
}

GLuint ResourceManager::CreateShader(GLenum type) {
  GLuint handle = shaderHandleAllocator_.Allocate();

  switch (type) {
  case GL_VERTEX_SHADER:
    shaderMap_[handle] = new VertexShader;
    break;

  // case GL_GEOMETRY_SHADER:
  //    shaderMap_[handle] = new GeometryShader;
  //    break;

  case GL_FRAGMENT_SHADER:
    shaderMap_[handle] = new FragmentShader;
    break;

  default:
    assert(0 && "unsupported shader type");
    break;
  }

  return handle;
}

GLuint ResourceManager::CreateProgram() {
  GLuint handle = programHandleAllocator_.Allocate();
  programMap_[handle] = new Program;
  return handle;
}

GLuint ResourceManager::CreateFramebuffer() {
  GLuint handle = framebufferHandleAllocator_.Allocate();
  framebufferMap_[handle] = new Framebuffer;
  return handle;
}

GLuint ResourceManager::CreateRenderbuffer() {
  GLuint handle = renderbufferHandleAllocator_.Allocate();
  renderbufferMap_[handle] = new Renderbuffer;
  return handle;
}

const Buffer *ResourceManager::GetBuffer(GLuint buffer) const {
  BufferMap::const_iterator it = bufferMap_.find(buffer);
  if (it == bufferMap_.end()) {
    fprintf(stderr, "ERR: cannot find buffer object for handle: %d\n", buffer);
    return NULL;
  }
  return it->second;
}

const Shader *ResourceManager::GetShader(GLuint shader) const {
  ShaderMap::const_iterator it = shaderMap_.find(shader);
  if (it == shaderMap_.end()) {
    fprintf(stderr, "ERR: cannot find shader object for handle: %d\n", shader);
    return NULL;
  }
  return it->second;
}

const Program *ResourceManager::GetProgram(GLuint program) const {
  ProgramMap::const_iterator it = programMap_.find(program);
  if (it == programMap_.end()) {
    fprintf(stderr, "ERR: cannot find program object for handle: %d\n",
            program);
    return NULL;
  }
  return it->second;
}

const Texture *ResourceManager::GetTexture(GLuint texture) const {
  TextureMap::const_iterator it = textureMap_.find(texture);
  if (it == textureMap_.end()) {
    fprintf(stderr, "ERR: cannot find texture object for handle: %d\n",
            texture);
    return NULL;
  }
  return it->second;
}

const Framebuffer *ResourceManager::GetFramebuffer(GLuint framebuffer) const {
  FramebufferMap::const_iterator it = framebufferMap_.find(framebuffer);
  if (it == framebufferMap_.end()) {
    fprintf(stderr, "ERR: cannot find framebuffer object for handle: %d\n",
            framebuffer);
    return NULL;
  }
  return it->second;
}

const Renderbuffer *
ResourceManager::GetRenderbuffer(GLuint renderbuffer) const {
  RenderbufferMap::const_iterator it = renderbufferMap_.find(renderbuffer);
  if (it == renderbufferMap_.end()) {
    fprintf(stderr, "ERR: cannot find renderbuffer object for handle: %d\n",
            renderbuffer);
    return NULL;
  }
  return it->second;
}

bool ResourceManager::IsValidBuffer(GLuint buffer) const {
  return bufferMap_.find(buffer) != bufferMap_.end();
}

bool ResourceManager::IsValidShader(GLuint shader) const {
  return shaderMap_.find(shader) != shaderMap_.end();
}

bool ResourceManager::IsValidProgram(GLuint program) const {
  return programMap_.find(program) != programMap_.end();
}

bool ResourceManager::IsValidTexture(GLuint texture) const {
  return textureMap_.find(texture) != textureMap_.end();
}

bool ResourceManager::IsValidFramebuffer(GLuint framebuffer) const {
  return framebufferMap_.find(framebuffer) != framebufferMap_.end();
}

bool ResourceManager::IsValidRenderbuffer(GLuint renderbuffer) const {
  return renderbufferMap_.find(renderbuffer) != renderbufferMap_.end();
}

GLuint ResourceManager::GetBufferHandle(Buffer *buffer) {
  BufferMap::iterator it = bufferMap_.begin();
  while (it != bufferMap_.end()) {
    if (it->second == buffer) {
      return it->first;
    }
  }
  return 0;
}

GLuint ResourceManager::GetShaderHandle(Shader *shader) {
  ShaderMap::iterator it = shaderMap_.begin();
  while (it != shaderMap_.end()) {
    if (it->second == shader) {
      return it->first;
    }
  }
  return 0;
}

GLuint ResourceManager::GetProgramHandle(Program *program) {
  ProgramMap::iterator it = programMap_.begin();
  while (it != programMap_.end()) {
    if (it->second == program) {
      return it->first;
    }
  }
  return 0;
}

GLuint ResourceManager::GetTextureHandle(Texture *texture) {
  TextureMap::iterator it = textureMap_.begin();
  while (it != textureMap_.end()) {
    if (it->second == texture) {
      return it->first;
    }
  }
  return 0;
}

GLuint ResourceManager::GetFramebufferHandle(Framebuffer *framebuffer) {
  FramebufferMap::iterator it = framebufferMap_.begin();
  while (it != framebufferMap_.end()) {
    if (it->second == framebuffer) {
      return it->first;
    }
  }
  return 0;
}

GLuint ResourceManager::GetRenderbufferHandle(Renderbuffer *renderbuffer) {
  RenderbufferMap::iterator it = renderbufferMap_.begin();
  while (it != renderbufferMap_.end()) {
    if (it->second == renderbuffer) {
      return it->first;
    }
  }
  return 0;
}

void ResourceManager::DeleteBuffer(GLuint buffer) {
  Buffer *buf = GetBuffer(buffer);
  if (buf == NULL) {
    fprintf(stderr,
            "ERR: cannot delete unknown buffer object with handle: %d\n",
            buffer);
    return;
  }

  delete buf;
  bufferHandleAllocator_.Release(buffer);
}

void ResourceManager::DeleteShader(GLuint shader) {
  Shader *shd = GetShader(shader);
  if (shd == NULL) {
    fprintf(stderr,
            "ERR: cannot delete unknown shader object with handle: %d\n",
            shader);
    return;
  }

  delete shd;
  shaderHandleAllocator_.Release(shader);
}

void ResourceManager::DeleteProgram(GLuint program) {
  Program *prg = GetProgram(program);
  if (prg == NULL) {
    fprintf(stderr,
            "ERR: cannot delete unknown program object with handle: %d\n",
            program);
    return;
  }

  delete prg;
  programHandleAllocator_.Release(program);
}

void ResourceManager::DeleteTexture(GLuint texture) {
  Texture *tex = GetTexture(texture);
  if (tex == NULL) {
    fprintf(stderr,
            "ERR: cannot delete unknown texture object with handle: %d\n",
            texture);
    return;
  }

  delete tex;
  textureHandleAllocator_.Release(texture);
}

void ResourceManager::DeleteFramebuffer(GLuint framebuffer) {
  Framebuffer *buf = GetFramebuffer(framebuffer);
  if (buf == NULL) {
    fprintf(stderr,
            "ERR: cannot delete unknown framebuffer object with handle: %d\n",
            framebuffer);
    return;
  }

  delete buf;
  framebufferHandleAllocator_.Release(framebuffer);
}

void ResourceManager::DeleteRenderbuffer(GLuint renderbuffer) {
  Renderbuffer *buf = GetRenderbuffer(renderbuffer);
  if (buf == NULL) {
    fprintf(stderr,
            "ERR: cannot delete unknown renderbuffer object with handle: %d\n",
            renderbuffer);
    return;
  }

  delete buf;
  renderbufferHandleAllocator_.Release(renderbuffer);
}
