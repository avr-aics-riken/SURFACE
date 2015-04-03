/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

//#define GL_APICALL  // prevent dllimport
#include "GLES2/gl2.h"
#include "gles_context.h"

//
// LSGL extenstions.
//
extern "C" {
GL_APICALL void GL_APIENTRY lsglSetPointSize(GLfloat size);

GL_APICALL void GL_APIENTRY lsglSetPointSizev(GLsizei num, const GLfloat *size);

GL_APICALL void GL_APIENTRY
    lsglSetCamera(GLfloat *eye, GLfloat *target, GLfloat *up, GLfloat fov);

GL_APICALL void GL_APIENTRY lsglInvalidateBuffer(GLenum target);

GL_APICALL void GL_APIENTRY lsglEvalFragmentShader();

GL_APICALL void GL_APIENTRY lsglEvalSingleFragmentShader();
};

using namespace lsgl;
// extern "C" {

void GL_APIENTRY glActiveTexture(GLenum texture) {
  Context::GetCurrentContext().glActiveTexture(texture);
}

void GL_APIENTRY glAttachShader(GLuint program, GLuint shader) {
  Context::GetCurrentContext().glAttachShader(program, shader);
}

void GL_APIENTRY
    glBindAttribLocation(GLuint program, GLuint index, const GLchar *name) {
  Context::GetCurrentContext().glBindAttribLocation(program, index, name);
}

void GL_APIENTRY glBindBuffer(GLenum target, GLuint buffer) {
  Context::GetCurrentContext().glBindBuffer(target, buffer);
}

void GL_APIENTRY glBindFramebuffer(GLenum target, GLuint framebuffer) {
  Context::GetCurrentContext().glBindFramebuffer(target, framebuffer);
}

void GL_APIENTRY glBindRenderbuffer(GLenum target, GLuint renderbuffer) {
  Context::GetCurrentContext().glBindRenderbuffer(target, renderbuffer);
}

void GL_APIENTRY glBindTexture(GLenum target, GLuint texture) {
  Context::GetCurrentContext().glBindTexture(target, texture);
}

void GL_APIENTRY
    glBlendColor(GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha) {
  Context::GetCurrentContext().glBlendColor(red, green, blue, alpha);
}

void GL_APIENTRY glBlendEquation(GLenum mode) {
  Context::GetCurrentContext().glBlendEquation(mode);
}

void GL_APIENTRY glBlendEquationSeparate(GLenum modeRGB, GLenum modeAlpha) {
  Context::GetCurrentContext().glBlendEquationSeparate(modeRGB, modeAlpha);
}

void GL_APIENTRY glBlendFunc(GLenum sfactor, GLenum dfactor) {
  Context::GetCurrentContext().glBlendFunc(sfactor, dfactor);
}

void GL_APIENTRY glBlendFuncSeparate(GLenum srcRGB, GLenum dstRGB,
                                     GLenum srcAlpha, GLenum dstAlpha) {
  Context::GetCurrentContext().glBlendFuncSeparate(srcRGB, dstRGB, srcAlpha,
                                                   dstAlpha);
}

void GL_APIENTRY glBufferData(GLenum target, GLsizeiptr size,
                              const GLvoid *data, GLenum usage) {
  Context::GetCurrentContext().glBufferData(target, size, data, usage);
}

void GL_APIENTRY glBufferSubData(GLenum target, GLintptr offset,
                                 GLsizeiptr size, const GLvoid *data) {
  Context::GetCurrentContext().glBufferSubData(target, offset, size, data);
}

GLenum GL_APIENTRY glCheckFramebufferStatus(GLenum target) {
  return Context::GetCurrentContext().glCheckFramebufferStatus(target);
}

void GL_APIENTRY glClear(GLbitfield mask) {
  Context::GetCurrentContext().glClear(mask);
}

void GL_APIENTRY
    glClearColor(GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha) {
  Context::GetCurrentContext().glClearColor(red, green, blue, alpha);
}

void GL_APIENTRY glClearDepthf(GLclampf depth) {
  Context::GetCurrentContext().glClearDepthf(depth);
}

void GL_APIENTRY glClearStencil(GLint s) {
  Context::GetCurrentContext().glClearStencil(s);
}

void GL_APIENTRY glColorMask(GLboolean red, GLboolean green, GLboolean blue,
                             GLboolean alpha) {
  Context::GetCurrentContext().glColorMask(red, green, blue, alpha);
}

void GL_APIENTRY glCompileShader(GLuint shader) {
  Context::GetCurrentContext().glCompileShader(shader);
}

void GL_APIENTRY glCompressedTexImage2D(GLenum target, GLint level,
                                        GLenum internalformat, GLsizei width,
                                        GLsizei height, GLint border,
                                        GLsizei imageSize, const GLvoid *data) {
  Context::GetCurrentContext().glCompressedTexImage2D(
      target, level, internalformat, width, height, border, imageSize, data);
}

void GL_APIENTRY glCompressedTexSubImage2D(GLenum target, GLint level,
                                           GLint xoffset, GLint yoffset,
                                           GLsizei width, GLsizei height,
                                           GLenum format, GLsizei imageSize,
                                           const GLvoid *data) {
  Context::GetCurrentContext().glCompressedTexSubImage2D(
      target, level, xoffset, yoffset, width, height, format, imageSize, data);
}

void GL_APIENTRY glCopyTexImage2D(GLenum target, GLint level,
                                  GLenum internalformat, GLint x, GLint y,
                                  GLsizei width, GLsizei height, GLint border) {
  Context::GetCurrentContext().glCopyTexImage2D(target, level, internalformat,
                                                x, y, width, height, border);
}

void GL_APIENTRY glCopyTexSubImage2D(GLenum target, GLint level, GLint xoffset,
                                     GLint yoffset, GLint x, GLint y,
                                     GLsizei width, GLsizei height) {
  Context::GetCurrentContext().glCopyTexSubImage2D(
      target, level, xoffset, yoffset, x, y, width, height);
}

GLuint GL_APIENTRY glCreateProgram(void) {
  return Context::GetCurrentContext().glCreateProgram();
}

GLuint GL_APIENTRY glCreateShader(GLenum type) {
  return Context::GetCurrentContext().glCreateShader(type);
}

void GL_APIENTRY glCullFace(GLenum mode) {
  Context::GetCurrentContext().glCullFace(mode);
}

void GL_APIENTRY glDeleteBuffers(GLsizei n, const GLuint *buffers) {
  Context::GetCurrentContext().glDeleteBuffers(n, buffers);
}

void GL_APIENTRY glDeleteFramebuffers(GLsizei n, const GLuint *framebuffers) {
  Context::GetCurrentContext().glDeleteFramebuffers(n, framebuffers);
}

void GL_APIENTRY glDeleteProgram(GLuint program) {
  Context::GetCurrentContext().glDeleteProgram(program);
}

void GL_APIENTRY glDeleteRenderbuffers(GLsizei n, const GLuint *renderbuffers) {
  Context::GetCurrentContext().glDeleteRenderbuffers(n, renderbuffers);
}

void GL_APIENTRY glDeleteShader(GLuint shader) {
  Context::GetCurrentContext().glDeleteShader(shader);
}

void GL_APIENTRY glDeleteTextures(GLsizei n, const GLuint *textures) {
  Context::GetCurrentContext().glDeleteTextures(n, textures);
}

void GL_APIENTRY glDepthFunc(GLenum func) {
  Context::GetCurrentContext().glDepthFunc(func);
}

void GL_APIENTRY glDepthMask(GLboolean flag) {
  Context::GetCurrentContext().glDepthMask(flag);
}

void GL_APIENTRY glDepthRangef(GLclampf zNear, GLclampf zFar) {
  Context::GetCurrentContext().glDepthRangef(zNear, zFar);
}

void GL_APIENTRY glDetachShader(GLuint program, GLuint shader) {
  Context::GetCurrentContext().glDetachShader(program, shader);
}

void GL_APIENTRY glDisable(GLenum cap) {
  Context::GetCurrentContext().glDisable(cap);
}

void GL_APIENTRY glDisableVertexAttribArray(GLuint index) {
  Context::GetCurrentContext().glDisableVertexAttribArray(index);
}

void GL_APIENTRY glDrawArrays(GLenum mode, GLint first, GLsizei count) {
  Context::GetCurrentContext().glDrawArrays(mode, first, count);
}

void GL_APIENTRY glDrawElements(GLenum mode, GLsizei count, GLenum type,
                                const GLvoid *indices) {
  Context::GetCurrentContext().glDrawElements(mode, count, type, indices);
}

void GL_APIENTRY glEnable(GLenum cap) {
  Context::GetCurrentContext().glEnable(cap);
}

void GL_APIENTRY glEnableVertexAttribArray(GLuint index) {
  Context::GetCurrentContext().glEnableVertexAttribArray(index);
}

void GL_APIENTRY glFinish(void) { Context::GetCurrentContext().glFinish(); }

void GL_APIENTRY glFlush(void) { Context::GetCurrentContext().glFlush(); }

void GL_APIENTRY glFramebufferRenderbuffer(GLenum target, GLenum attachment,
                                           GLenum renderbuffertarget,
                                           GLuint renderbuffer) {
  Context::GetCurrentContext().glFramebufferRenderbuffer(
      target, attachment, renderbuffertarget, renderbuffer);
}

void GL_APIENTRY glFramebufferTexture2D(GLenum target, GLenum attachment,
                                        GLenum textarget, GLuint texture,
                                        GLint level) {
  Context::GetCurrentContext().glFramebufferTexture2D(
      target, attachment, textarget, texture, level);
}

void GL_APIENTRY glFrontFace(GLenum mode) {
  Context::GetCurrentContext().glFrontFace(mode);
}

void GL_APIENTRY glGenBuffers(GLsizei n, GLuint *buffers) {
  Context::GetCurrentContext().glGenBuffers(n, buffers);
}

void GL_APIENTRY glGenerateMipmap(GLenum target) {
  Context::GetCurrentContext().glGenerateMipmap(target);
}

void GL_APIENTRY glGenFramebuffers(GLsizei n, GLuint *framebuffers) {
  Context::GetCurrentContext().glGenFramebuffers(n, framebuffers);
}

void GL_APIENTRY glGenRenderbuffers(GLsizei n, GLuint *renderbuffers) {
  Context::GetCurrentContext().glGenRenderbuffers(n, renderbuffers);
}

void GL_APIENTRY glGenTextures(GLsizei n, GLuint *textures) {
  Context::GetCurrentContext().glGenTextures(n, textures);
}

void GL_APIENTRY glGetActiveAttrib(GLuint program, GLuint index,
                                   GLsizei bufsize, GLsizei *length,
                                   GLint *size, GLenum *type, GLchar *name) {
  Context::GetCurrentContext().glGetActiveAttrib(program, index, bufsize,
                                                 length, size, type, name);
}

void GL_APIENTRY glGetActiveUniform(GLuint program, GLuint index,
                                    GLsizei bufsize, GLsizei *length,
                                    GLint *size, GLenum *type, GLchar *name) {
  Context::GetCurrentContext().glGetActiveUniform(program, index, bufsize,
                                                  length, size, type, name);
}

void GL_APIENTRY glGetAttachedShaders(GLuint program, GLsizei maxcount,
                                      GLsizei *count, GLuint *shaders) {
  Context::GetCurrentContext().glGetAttachedShaders(program, maxcount, count,
                                                    shaders);
}

int GL_APIENTRY glGetAttribLocation(GLuint program, const GLchar *name) {
  return Context::GetCurrentContext().glGetAttribLocation(program, name);
}

void GL_APIENTRY glGetBooleanv(GLenum pname, GLboolean *params) {
  Context::GetCurrentContext().glGetBooleanv(pname, params);
}

void GL_APIENTRY
    glGetBufferParameteriv(GLenum target, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetBufferParameteriv(target, pname, params);
}

GLenum GL_APIENTRY glGetError(void) {
  return Context::GetCurrentContext().glGetError();
}

void GL_APIENTRY glGetFloatv(GLenum pname, GLfloat *params) {
  Context::GetCurrentContext().glGetFloatv(pname, params);
}

void GL_APIENTRY
    glGetFramebufferAttachmentParameteriv(GLenum target, GLenum attachment,
                                          GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetFramebufferAttachmentParameteriv(
      target, attachment, pname, params);
}

void GL_APIENTRY glGetIntegerv(GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetIntegerv(pname, params);
}

void GL_APIENTRY glGetProgramiv(GLuint program, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetProgramiv(program, pname, params);
}

void GL_APIENTRY glGetProgramInfoLog(GLuint program, GLsizei bufsize,
                                     GLsizei *length, GLchar *infolog) {
  Context::GetCurrentContext().glGetProgramInfoLog(program, bufsize, length,
                                                   infolog);
}

void GL_APIENTRY
    glGetRenderbufferParameteriv(GLenum target, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetRenderbufferParameteriv(target, pname,
                                                            params);
}

void GL_APIENTRY glGetShaderiv(GLuint shader, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetShaderiv(shader, pname, params);
}

void GL_APIENTRY glGetShaderInfoLog(GLuint shader, GLsizei bufsize,
                                    GLsizei *length, GLchar *infolog) {
  Context::GetCurrentContext().glGetShaderInfoLog(shader, bufsize, length,
                                                  infolog);
}

void GL_APIENTRY glGetShaderPrecisionFormat(GLenum shadertype,
                                            GLenum precisiontype, GLint *range,
                                            GLint *precision) {
  Context::GetCurrentContext().glGetShaderPrecisionFormat(
      shadertype, precisiontype, range, precision);
}

void GL_APIENTRY glGetShaderSource(GLuint shader, GLsizei bufsize,
                                   GLsizei *length, GLchar *source) {
  Context::GetCurrentContext().glGetShaderSource(shader, bufsize, length,
                                                 source);
}

const GLubyte *GL_APIENTRY glGetString(GLenum name) {
  return Context::GetCurrentContext().glGetString(name);
}

void GL_APIENTRY
    glGetTexParameterfv(GLenum target, GLenum pname, GLfloat *params) {
  Context::GetCurrentContext().glGetTexParameterfv(target, pname, params);
}

void GL_APIENTRY
    glGetTexParameteriv(GLenum target, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetTexParameteriv(target, pname, params);
}

void GL_APIENTRY
    glGetUniformfv(GLuint program, GLint location, GLfloat *params) {
  Context::GetCurrentContext().glGetUniformfv(program, location, params);
}

void GL_APIENTRY glGetUniformiv(GLuint program, GLint location, GLint *params) {
  Context::GetCurrentContext().glGetUniformiv(program, location, params);
}

int GL_APIENTRY glGetUniformLocation(GLuint program, const GLchar *name) {
  return Context::GetCurrentContext().glGetUniformLocation(program, name);
}

void GL_APIENTRY
    glGetVertexAttribfv(GLuint index, GLenum pname, GLfloat *params) {
  Context::GetCurrentContext().glGetVertexAttribfv(index, pname, params);
}

void GL_APIENTRY
    glGetVertexAttribiv(GLuint index, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetVertexAttribiv(index, pname, params);
}

void GL_APIENTRY
    glGetVertexAttribPointerv(GLuint index, GLenum pname, GLvoid **pointer) {
  Context::GetCurrentContext().glGetVertexAttribPointerv(index, pname, pointer);
}

void GL_APIENTRY glHint(GLenum target, GLenum mode) {
  Context::GetCurrentContext().glHint(target, mode);
}

GLboolean GL_APIENTRY glIsBuffer(GLuint buffer) {
  return Context::GetCurrentContext().glIsBuffer(buffer);
}

GLboolean GL_APIENTRY glIsEnabled(GLenum cap) {
  return Context::GetCurrentContext().glIsEnabled(cap);
}

GLboolean GL_APIENTRY glIsFramebuffer(GLuint framebuffer) {
  return Context::GetCurrentContext().glIsFramebuffer(framebuffer);
}

GLboolean GL_APIENTRY glIsProgram(GLuint program) {
  return Context::GetCurrentContext().glIsProgram(program);
}

GLboolean GL_APIENTRY glIsRenderbuffer(GLuint renderbuffer) {
  return Context::GetCurrentContext().glIsRenderbuffer(renderbuffer);
}

GLboolean GL_APIENTRY glIsShader(GLuint shader) {
  return Context::GetCurrentContext().glIsShader(shader);
}

GLboolean GL_APIENTRY glIsTexture(GLuint texture) {
  return Context::GetCurrentContext().glIsTexture(texture);
}

void GL_APIENTRY glLineWidth(GLfloat width) {
  Context::GetCurrentContext().glLineWidth(width);
}

void GL_APIENTRY glLinkProgram(GLuint program) {
  Context::GetCurrentContext().glLinkProgram(program);
}

void GL_APIENTRY glPixelStorei(GLenum pname, GLint param) {
  Context::GetCurrentContext().glPixelStorei(pname, param);
}

void GL_APIENTRY glPolygonOffset(GLfloat factor, GLfloat units) {
  Context::GetCurrentContext().glPolygonOffset(factor, units);
}

void GL_APIENTRY glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height,
                              GLenum format, GLenum type, GLvoid *pixels) {
  Context::GetCurrentContext().glReadPixels(x, y, width, height, format, type,
                                            pixels);
}

void GL_APIENTRY glReleaseShaderCompiler(void) {
  Context::GetCurrentContext().glReleaseShaderCompiler();
}

void GL_APIENTRY glRenderbufferStorage(GLenum target, GLenum internalformat,
                                       GLsizei width, GLsizei height) {
  Context::GetCurrentContext().glRenderbufferStorage(target, internalformat,
                                                     width, height);
}

void GL_APIENTRY glSampleCoverage(GLclampf value, GLboolean invert) {
  Context::GetCurrentContext().glSampleCoverage(value, invert);
}

void GL_APIENTRY glScissor(GLint x, GLint y, GLsizei width, GLsizei height) {
  Context::GetCurrentContext().glScissor(x, y, width, height);
}

void GL_APIENTRY glShaderBinary(GLsizei n, const GLuint *shaders,
                                GLenum binaryformat, const GLvoid *binary,
                                GLsizei length) {
  Context::GetCurrentContext().glShaderBinary(n, shaders, binaryformat, binary,
                                              length);
}

void GL_APIENTRY glShaderSource(GLuint shader, GLsizei count,
                                const GLchar **string, const GLint *length) {
  Context::GetCurrentContext().glShaderSource(shader, count, string, length);
}

void GL_APIENTRY glStencilFunc(GLenum func, GLint ref, GLuint mask) {
  Context::GetCurrentContext().glStencilFunc(func, ref, mask);
}

void GL_APIENTRY
    glStencilFuncSeparate(GLenum face, GLenum func, GLint ref, GLuint mask) {
  Context::GetCurrentContext().glStencilFuncSeparate(face, func, ref, mask);
}

void GL_APIENTRY glStencilMask(GLuint mask) {
  Context::GetCurrentContext().glStencilMask(mask);
}

void GL_APIENTRY glStencilMaskSeparate(GLenum face, GLuint mask) {
  Context::GetCurrentContext().glStencilMaskSeparate(face, mask);
}

void GL_APIENTRY glStencilOp(GLenum fail, GLenum zfail, GLenum zpass) {
  Context::GetCurrentContext().glStencilOp(fail, zfail, zpass);
}

void GL_APIENTRY
    glStencilOpSeparate(GLenum face, GLenum fail, GLenum zfail, GLenum zpass) {
  Context::GetCurrentContext().glStencilOpSeparate(face, fail, zfail, zpass);
}

void GL_APIENTRY glTexImage2D(GLenum target, GLint level, GLint internalformat,
                              GLsizei width, GLsizei height, GLint border,
                              GLenum format, GLenum type,
                              const GLvoid *pixels) {
  Context::GetCurrentContext().glTexImage2D(target, level, internalformat,
                                            width, height, border, format, type,
                                            pixels);
}

void GL_APIENTRY glTexImage3D(GLenum target, GLint level, GLint internalformat,
                              GLsizei width, GLsizei height, GLsizei depth,
                              GLint border, GLenum format, GLenum type,
                              const GLvoid *pixels) {
  Context::GetCurrentContext().glTexImage3D(target, level, internalformat,
                                            width, height, depth, border,
                                            format, type, pixels);
}

void GL_APIENTRY glTexParameterf(GLenum target, GLenum pname, GLfloat param) {
  Context::GetCurrentContext().glTexParameterf(target, pname, param);
}

void GL_APIENTRY
    glTexParameterfv(GLenum target, GLenum pname, const GLfloat *params) {
  Context::GetCurrentContext().glTexParameterfv(target, pname, params);
}

void GL_APIENTRY glTexParameteri(GLenum target, GLenum pname, GLint param) {
  Context::GetCurrentContext().glTexParameteri(target, pname, param);
}

void GL_APIENTRY
    glTexParameteriv(GLenum target, GLenum pname, const GLint *params) {
  Context::GetCurrentContext().glTexParameteriv(target, pname, params);
}

void GL_APIENTRY glTexSubImage2D(GLenum target, GLint level, GLint xoffset,
                                 GLint yoffset, GLsizei width, GLsizei height,
                                 GLenum format, GLenum type,
                                 const GLvoid *pixels) {
  Context::GetCurrentContext().glTexSubImage2D(
      target, level, xoffset, yoffset, width, height, format, type, pixels);
}

void GL_APIENTRY glUniform1f(GLint location, GLfloat x) {
  Context::GetCurrentContext().glUniform1f(location, x);
}

void GL_APIENTRY glUniform1fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform1fv(location, count, v);
}

void GL_APIENTRY glUniform1i(GLint location, GLint x) {
  Context::GetCurrentContext().glUniform1i(location, x);
}

void GL_APIENTRY glUniform1iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform1iv(location, count, v);
}

void GL_APIENTRY glUniform2f(GLint location, GLfloat x, GLfloat y) {
  Context::GetCurrentContext().glUniform2f(location, x, y);
}

void GL_APIENTRY glUniform2fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform2fv(location, count, v);
}

void GL_APIENTRY glUniform2i(GLint location, GLint x, GLint y) {
  Context::GetCurrentContext().glUniform2i(location, x, y);
}

void GL_APIENTRY glUniform2iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform2iv(location, count, v);
}

void GL_APIENTRY glUniform3f(GLint location, GLfloat x, GLfloat y, GLfloat z) {
  Context::GetCurrentContext().glUniform3f(location, x, y, z);
}

void GL_APIENTRY glUniform3fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform3fv(location, count, v);
}

void GL_APIENTRY glUniform3i(GLint location, GLint x, GLint y, GLint z) {
  Context::GetCurrentContext().glUniform3i(location, x, y, z);
}

void GL_APIENTRY glUniform3iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform3iv(location, count, v);
}

void GL_APIENTRY
    glUniform4f(GLint location, GLfloat x, GLfloat y, GLfloat z, GLfloat w) {
  Context::GetCurrentContext().glUniform4f(location, x, y, z, w);
}

void GL_APIENTRY glUniform4fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform4fv(location, count, v);
}

void GL_APIENTRY
    glUniform4i(GLint location, GLint x, GLint y, GLint z, GLint w) {
  Context::GetCurrentContext().glUniform4i(location, x, y, z, w);
}

void GL_APIENTRY glUniform4iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform4iv(location, count, v);
}

void GL_APIENTRY glUniformMatrix2fv(GLint location, GLsizei count,
                                    GLboolean transpose, const GLfloat *value) {
  Context::GetCurrentContext().glUniformMatrix2fv(location, count, transpose,
                                                  value);
}

void GL_APIENTRY glUniformMatrix3fv(GLint location, GLsizei count,
                                    GLboolean transpose, const GLfloat *value) {
  Context::GetCurrentContext().glUniformMatrix3fv(location, count, transpose,
                                                  value);
}

void GL_APIENTRY glUniformMatrix4fv(GLint location, GLsizei count,
                                    GLboolean transpose, const GLfloat *value) {
  Context::GetCurrentContext().glUniformMatrix4fv(location, count, transpose,
                                                  value);
}

void GL_APIENTRY glUseProgram(GLuint program) {
  Context::GetCurrentContext().glUseProgram(program);
}

void GL_APIENTRY glValidateProgram(GLuint program) {
  Context::GetCurrentContext().glValidateProgram(program);
}

void GL_APIENTRY glVertexAttrib1f(GLuint indx, GLfloat x) {
  Context::GetCurrentContext().glVertexAttrib1f(indx, x);
}

void GL_APIENTRY glVertexAttrib1fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib1fv(indx, values);
}

void GL_APIENTRY glVertexAttrib2f(GLuint indx, GLfloat x, GLfloat y) {
  Context::GetCurrentContext().glVertexAttrib2f(indx, x, y);
}

void GL_APIENTRY glVertexAttrib2fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib2fv(indx, values);
}

void GL_APIENTRY
    glVertexAttrib3f(GLuint indx, GLfloat x, GLfloat y, GLfloat z) {
  Context::GetCurrentContext().glVertexAttrib3f(indx, x, y, z);
}

void GL_APIENTRY glVertexAttrib3fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib3fv(indx, values);
}

void GL_APIENTRY
    glVertexAttrib4f(GLuint indx, GLfloat x, GLfloat y, GLfloat z, GLfloat w) {
  Context::GetCurrentContext().glVertexAttrib4f(indx, x, y, z, w);
}

void GL_APIENTRY glVertexAttrib4fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib4fv(indx, values);
}

void GL_APIENTRY glVertexAttribPointer(GLuint indx, GLint size, GLenum type,
                                       GLboolean normalized, GLsizei stride,
                                       const GLvoid *ptr) {
  Context::GetCurrentContext().glVertexAttribPointer(indx, size, type,
                                                     normalized, stride, ptr);
}

void GL_APIENTRY glViewport(GLint x, GLint y, GLsizei width, GLsizei height) {
  Context::GetCurrentContext().glViewport(x, y, width, height);
}

//
// LSGL extenstions.
//
GL_APICALL void GL_APIENTRY lsglSetPointSize(GLfloat size) {
  Context::GetCurrentContext().lsglSetPointSize(size);
}

GL_APICALL void GL_APIENTRY
    lsglSetPointSizev(GLsizei num, const GLfloat *size) {
  Context::GetCurrentContext().lsglSetPointSizev(num, size);
}

GL_APICALL void GL_APIENTRY
    lsglSetCamera(GLfloat *eye, GLfloat *target, GLfloat *up, GLfloat fov) {
  Context::GetCurrentContext().lsglSetCamera(eye, target, up, fov);
}

void GL_APIENTRY lsglInvalidateBuffer(GLenum target) {
  Context::GetCurrentContext().lsglInvalidateBuffer(target);
}

//
// LSGL test functions. Tentative, to be removed in the future
//
void GL_APIENTRY lsglEvalFragmentShader() {
  Context::GetCurrentContext().lsglEvalFragmentShader();
}

void GL_APIENTRY lsglEvalSingleFragmentShader() {
  Context::GetCurrentContext().lsglEvalSingleFragmentShader();
}

// } // namespace lsgl

//}
