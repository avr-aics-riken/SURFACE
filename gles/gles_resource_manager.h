/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_GLES_RESOURCE_MANAGER_H__
#define __LSGL_GLES_RESOURCE_MANAGER_H__

#include <vector>
#include <map>
#include <cassert>

#include "gles_common.h"
#include "gles_handle_allocator.h"

#include "GLES2/gl2.h"

namespace lsgl {

class ResourceManager {
public:
  ResourceManager() {};
  ~ResourceManager() {};

  GLuint CreateBuffer();
  GLuint CreateShader(GLenum type);
  GLuint CreateProgram();
  GLuint CreateTexture();
  GLuint CreateFramebuffer();
  GLuint CreateRenderbuffer();

  void DeleteBuffer(GLuint buffer);
  void DeleteShader(GLuint shader);
  void DeleteProgram(GLuint program);
  void DeleteTexture(GLuint texture);
  void DeleteFramebuffer(GLuint framebuffer);
  void DeleteRenderbuffer(GLuint renderbuffer);

  const Buffer *GetBuffer(GLuint buffer) const;
  const Shader *GetShader(GLuint shader) const;
  const Program *GetProgram(GLuint program) const;
  const Texture *GetTexture(GLuint texture) const;
  const Framebuffer *GetFramebuffer(GLuint framebuffer) const;
  const Renderbuffer *GetRenderbuffer(GLuint renderbuffer) const;

  // non-const calls simply call the const versions, then cast back
  inline Buffer *GetBuffer(GLuint buffer) {
    return const_cast<Buffer *>(
        const_cast<const ResourceManager *>(this)->GetBuffer(buffer));
  }
  inline Shader *GetShader(GLuint shader) {
    return const_cast<Shader *>(
        const_cast<const ResourceManager *>(this)->GetShader(shader));
  }
  inline Program *GetProgram(GLuint program) {
    return const_cast<Program *>(
        const_cast<const ResourceManager *>(this)->GetProgram(program));
  }
  inline Texture *GetTexture(GLuint texture) {
    return const_cast<Texture *>(
        const_cast<const ResourceManager *>(this)->GetTexture(texture));
  }
  inline Framebuffer *GetFramebuffer(GLuint framebuffer) {
    return const_cast<Framebuffer *>(
        const_cast<const ResourceManager *>(this)->GetFramebuffer(framebuffer));
  }
  inline Renderbuffer *GetRenderbuffer(GLuint renderbuffer) {
    return const_cast<Renderbuffer *>(const_cast<const ResourceManager *>(this)
                                          ->GetRenderbuffer(renderbuffer));
  }

  bool IsValidBuffer(GLuint buffer) const;
  bool IsValidShader(GLuint shader) const;
  bool IsValidProgram(GLuint program) const;
  bool IsValidTexture(GLuint texture) const;
  bool IsValidFramebuffer(GLuint framebuffer) const;
  bool IsValidRenderbuffer(GLuint renderbuffer) const;

  GLuint GetBufferHandle(Buffer *buffer);
  GLuint GetShaderHandle(Shader *shader);
  GLuint GetProgramHandle(Program *program);
  GLuint GetTextureHandle(Texture *texture);
  GLuint GetFramebufferHandle(Framebuffer *framebuffer);
  GLuint GetRenderbufferHandle(Renderbuffer *renderbuffer);

private:
  HandleAllocator bufferHandleAllocator_;
  HandleAllocator shaderHandleAllocator_;
  HandleAllocator programHandleAllocator_;
  HandleAllocator textureHandleAllocator_;
  HandleAllocator framebufferHandleAllocator_;
  HandleAllocator renderbufferHandleAllocator_;

  typedef std::map<GLuint, Buffer *> BufferMap;
  typedef std::map<GLuint, Shader *> ShaderMap;
  typedef std::map<GLuint, Program *> ProgramMap;
  typedef std::map<GLuint, Texture *> TextureMap;
  typedef std::map<GLuint, Framebuffer *> FramebufferMap;
  typedef std::map<GLuint, Renderbuffer *> RenderbufferMap;

  BufferMap bufferMap_;
  ShaderMap shaderMap_;
  ProgramMap programMap_;
  TextureMap textureMap_;
  FramebufferMap framebufferMap_;
  RenderbufferMap renderbufferMap_;
};
}

#endif // __LSGL_GLES_RESOURCE_MANAGER_H__
