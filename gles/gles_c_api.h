/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_GLES_C_API_H__
#define __LSGL_GLES_C_API_H__

#include "GLES2/gl2.h"
#include "gles_context.h"

using namespace lsgl;
extern "C" {

void glActiveTexture(GLenum texture) {
  Context::GetCurrentContext().glActiveTexture(texture);
}

void glAttachShader(GLuint program, GLuint shader) {
  Context::GetCurrentContext().glAttachShader(program, shader);
}

void glBindAttribLocation(GLuint program, GLuint index, const GLchar *name) {
  Context::GetCurrentContext().glBindAttribLocation(program, index, name);
}

void glBindBuffer(GLenum target, GLuint buffer) {
  Context::GetCurrentContext().glBindBuffer(target, buffer);
}

void glBindFramebuffer(GLenum target, GLuint framebuffer) {
  Context::GetCurrentContext().glBindFramebuffer(target, framebuffer);
}

void glBindRenderbuffer(GLenum target, GLuint renderbuffer) {
  Context::GetCurrentContext().glBindRenderbuffer(target, renderbuffer);
}

void glBindTexture(GLenum target, GLuint texture) {
  Context::GetCurrentContext().glBindTexture(target, texture);
}

void glBlendColor(GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha) {
  Context::GetCurrentContext().glBlendColor(red, green, blue, alpha);
}

void glBlendEquation(GLenum mode) {
  Context::GetCurrentContext().glBlendEquation(mode);
}

void glBlendEquationSeparate(GLenum modeRGB, GLenum modeAlpha) {
  Context::GetCurrentContext().glBlendEquationSeparate(modeRGB, modeAlpha);
}

void glBlendFunc(GLenum sfactor, GLenum dfactor) {
  Context::GetCurrentContext().glBlendFunc(sfactor, dfactor);
}

void glBlendFuncSeparate(GLenum srcRGB, GLenum dstRGB, GLenum srcAlpha,
                         GLenum dstAlpha) {
  Context::GetCurrentContext().glBlendFuncSeparate(srcRGB, dstRGB, srcAlpha,
                                                   dstAlpha);
}

void glBufferData(GLenum target, GLsizeiptr size, const GLvoid *data,
                  GLenum usage) {
  Context::GetCurrentContext().glBufferData(target, size, data, usage);
}

void glBufferSubData(GLenum target, GLintptr offset, GLsizeiptr size,
                     const GLvoid *data) {
  Context::GetCurrentContext().glBufferSubData(target, offset, size, data);
}

GLenum glCheckFramebufferStatus(GLenum target) {
  return Context::GetCurrentContext().glCheckFramebufferStatus(target);
}

void glClear(GLbitfield mask) { Context::GetCurrentContext().glClear(mask); }

void glClearColor(GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha) {
  Context::GetCurrentContext().glClearColor(red, green, blue, alpha);
}

void glClearDepthf(GLclampf depth) {
  Context::GetCurrentContext().glClearDepthf(depth);
}

void glClearStencil(GLint s) { Context::GetCurrentContext().glClearStencil(s); }

void glColorMask(GLboolean red, GLboolean green, GLboolean blue,
                 GLboolean alpha) {
  Context::GetCurrentContext().glColorMask(red, green, blue, alpha);
}

void glCompileShader(GLuint shader) {
  Context::GetCurrentContext().glCompileShader(shader);
}

void glCompressedTexImage2D(GLenum target, GLint level, GLenum internalformat,
                            GLsizei width, GLsizei height, GLint border,
                            GLsizei imageSize, const GLvoid *data) {
  Context::GetCurrentContext().glCompressedTexImage2D(
      target, level, internalformat, width, height, border, imageSize, data);
}

void glCompressedTexSubImage2D(GLenum target, GLint level, GLint xoffset,
                               GLint yoffset, GLsizei width, GLsizei height,
                               GLenum format, GLsizei imageSize,
                               const GLvoid *data) {
  Context::GetCurrentContext().glCompressedTexSubImage2D(
      target, level, xoffset, yoffset, width, height, format, imageSize, data);
}

void glCopyTexImage2D(GLenum target, GLint level, GLenum internalformat,
                      GLint x, GLint y, GLsizei width, GLsizei height,
                      GLint border) {
  Context::GetCurrentContext().glCopyTexImage2D(target, level, internalformat,
                                                x, y, width, height, border);
}

void glCopyTexSubImage2D(GLenum target, GLint level, GLint xoffset,
                         GLint yoffset, GLint x, GLint y, GLsizei width,
                         GLsizei height) {
  Context::GetCurrentContext().glCopyTexSubImage2D(
      target, level, xoffset, yoffset, x, y, width, height);
}

GLuint glCreateProgram(void) {
  return Context::GetCurrentContext().glCreateProgram();
}

GLuint glCreateShader(GLenum type) {
  return Context::GetCurrentContext().glCreateShader(type);
}

void glCullFace(GLenum mode) { Context::GetCurrentContext().glCullFace(mode); }

void glDeleteBuffers(GLsizei n, const GLuint *buffers) {
  Context::GetCurrentContext().glDeleteBuffers(n, buffers);
}

void glDeleteFramebuffers(GLsizei n, const GLuint *framebuffers) {
  Context::GetCurrentContext().glDeleteFramebuffers(n, framebuffers);
}

void glDeleteProgram(GLuint program) {
  Context::GetCurrentContext().glDeleteProgram(program);
}

void glDeleteRenderbuffers(GLsizei n, const GLuint *renderbuffers) {
  Context::GetCurrentContext().glDeleteRenderbuffers(n, renderbuffers);
}

void glDeleteShader(GLuint shader) {
  Context::GetCurrentContext().glDeleteShader(shader);
}

void glDeleteTextures(GLsizei n, const GLuint *textures) {
  Context::GetCurrentContext().glDeleteTextures(n, textures);
}

void glDepthFunc(GLenum func) {
  Context::GetCurrentContext().glDepthFunc(func);
}

void glDepthMask(GLboolean flag) {
  Context::GetCurrentContext().glDepthMask(flag);
}

void glDepthRangef(GLclampf zNear, GLclampf zFar) {
  Context::GetCurrentContext().glDepthRangef(zNear, zFar);
}

void glDetachShader(GLuint program, GLuint shader) {
  Context::GetCurrentContext().glDetachShader(program, shader);
}

void glDisable(GLenum cap) { Context::GetCurrentContext().glDisable(cap); }

void glDisableVertexAttribArray(GLuint index) {
  Context::GetCurrentContext().glDisableVertexAttribArray(index);
}

void glDrawArrays(GLenum mode, GLint first, GLsizei count) {
  Context::GetCurrentContext().glDrawArrays(mode, first, count);
}

void glDrawElements(GLenum mode, GLsizei count, GLenum type,
                    const GLvoid *indices) {
  Context::GetCurrentContext().glDrawElements(mode, count, type, indices);
}

void glEnable(GLenum cap) { Context::GetCurrentContext().glEnable(cap); }

void glEnableVertexAttribArray(GLuint index) {
  Context::GetCurrentContext().glEnableVertexAttribArray(index);
}

void glFinish(void) { Context::GetCurrentContext().glFinish(); }

void glFlush(void) { Context::GetCurrentContext().glFlush(); }

void glFramebufferRenderbuffer(GLenum target, GLenum attachment,
                               GLenum renderbuffertarget, GLuint renderbuffer) {
  Context::GetCurrentContext().glFramebufferRenderbuffer(
      target, attachment, renderbuffertarget, renderbuffer);
}

void glFramebufferTexture2D(GLenum target, GLenum attachment, GLenum textarget,
                            GLuint texture, GLint level) {
  Context::GetCurrentContext().glFramebufferTexture2D(
      target, attachment, textarget, texture, level);
}

void glFrontFace(GLenum mode) {
  Context::GetCurrentContext().glFrontFace(mode);
}

void glGenBuffers(GLsizei n, GLuint *buffers) {
  Context::GetCurrentContext().glGenBuffers(n, buffers);
}

void glGenerateMipmap(GLenum target) {
  Context::GetCurrentContext().glGenerateMipmap(target);
}

void glGenFramebuffers(GLsizei n, GLuint *framebuffers) {
  Context::GetCurrentContext().glGenFramebuffers(n, framebuffers);
}

void glGenRenderbuffers(GLsizei n, GLuint *renderbuffers) {
  Context::GetCurrentContext().glGenRenderbuffers(n, renderbuffers);
}

void glGenTextures(GLsizei n, GLuint *textures) {
  Context::GetCurrentContext().glGenTextures(n, textures);
}

void glGetActiveAttrib(GLuint program, GLuint index, GLsizei bufsize,
                       GLsizei *length, GLint *size, GLenum *type,
                       GLchar *name) {
  Context::GetCurrentContext().glGetActiveAttrib(program, index, bufsize,
                                                 length, size, type, name);
}

void glGetActiveUniform(GLuint program, GLuint index, GLsizei bufsize,
                        GLsizei *length, GLint *size, GLenum *type,
                        GLchar *name) {
  Context::GetCurrentContext().glGetActiveUniform(program, index, bufsize,
                                                  length, size, type, name);
}

void glGetAttachedShaders(GLuint program, GLsizei maxcount, GLsizei *count,
                          GLuint *shaders) {
  Context::GetCurrentContext().glGetAttachedShaders(program, maxcount, count,
                                                    shaders);
}

int glGetAttribLocation(GLuint program, const GLchar *name) {
  return Context::GetCurrentContext().glGetAttribLocation(program, name);
}

void glGetBooleanv(GLenum pname, GLboolean *params) {
  Context::GetCurrentContext().glGetBooleanv(pname, params);
}

void glGetBufferParameteriv(GLenum target, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetBufferParameteriv(target, pname, params);
}

GLenum glGetError(void) { return Context::GetCurrentContext().glGetError(); }

void glGetFloatv(GLenum pname, GLfloat *params) {
  Context::GetCurrentContext().glGetFloatv(pname, params);
}

void glGetFramebufferAttachmentParameteriv(GLenum target, GLenum attachment,
                                           GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetFramebufferAttachmentParameteriv(
      target, attachment, pname, params);
}

void glGetIntegerv(GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetIntegerv(pname, params);
}

void glGetProgramiv(GLuint program, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetProgramiv(program, pname, params);
}

void glGetProgramInfoLog(GLuint program, GLsizei bufsize, GLsizei *length,
                         GLchar *infolog) {
  Context::GetCurrentContext().glGetProgramInfoLog(program, bufsize, length,
                                                   infolog);
}

void glGetRenderbufferParameteriv(GLenum target, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetRenderbufferParameteriv(target, pname,
                                                            params);
}

void glGetShaderiv(GLuint shader, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetShaderiv(shader, pname, params);
}

void glGetShaderInfoLog(GLuint shader, GLsizei bufsize, GLsizei *length,
                        GLchar *infolog) {
  Context::GetCurrentContext().glGetShaderInfoLog(shader, bufsize, length,
                                                  infolog);
}

void glGetShaderPrecisionFormat(GLenum shadertype, GLenum precisiontype,
                                GLint *range, GLint *precision) {
  Context::GetCurrentContext().glGetShaderPrecisionFormat(
      shadertype, precisiontype, range, precision);
}

void glGetShaderSource(GLuint shader, GLsizei bufsize, GLsizei *length,
                       GLchar *source) {
  Context::GetCurrentContext().glGetShaderSource(shader, bufsize, length,
                                                 source);
}

const GLubyte *glGetString(GLenum name) {
  return Context::GetCurrentContext().glGetString(name);
}

void glGetTexParameterfv(GLenum target, GLenum pname, GLfloat *params) {
  Context::GetCurrentContext().glGetTexParameterfv(target, pname, params);
}

void glGetTexParameteriv(GLenum target, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetTexParameteriv(target, pname, params);
}

void glGetUniformfv(GLuint program, GLint location, GLfloat *params) {
  Context::GetCurrentContext().glGetUniformfv(program, location, params);
}

void glGetUniformiv(GLuint program, GLint location, GLint *params) {
  Context::GetCurrentContext().glGetUniformiv(program, location, params);
}

int glGetUniformLocation(GLuint program, const GLchar *name) {
  return Context::GetCurrentContext().glGetUniformLocation(program, name);
}

void glGetVertexAttribfv(GLuint index, GLenum pname, GLfloat *params) {
  Context::GetCurrentContext().glGetVertexAttribfv(index, pname, params);
}

void glGetVertexAttribiv(GLuint index, GLenum pname, GLint *params) {
  Context::GetCurrentContext().glGetVertexAttribiv(index, pname, params);
}

void glGetVertexAttribPointerv(GLuint index, GLenum pname, GLvoid **pointer) {
  Context::GetCurrentContext().glGetVertexAttribPointerv(index, pname, pointer);
}

void glHint(GLenum target, GLenum mode) {
  Context::GetCurrentContext().glHint(target, mode);
}

GLboolean glIsBuffer(GLuint buffer) {
  return Context::GetCurrentContext().glIsBuffer(buffer);
}

GLboolean glIsEnabled(GLenum cap) {
  return Context::GetCurrentContext().glIsEnabled(cap);
}

GLboolean glIsFramebuffer(GLuint framebuffer) {
  return Context::GetCurrentContext().glIsFramebuffer(framebuffer);
}

GLboolean glIsProgram(GLuint program) {
  return Context::GetCurrentContext().glIsProgram(program);
}

GLboolean glIsRenderbuffer(GLuint renderbuffer) {
  return Context::GetCurrentContext().glIsRenderbuffer(renderbuffer);
}

GLboolean glIsShader(GLuint shader) {
  return Context::GetCurrentContext().glIsShader(shader);
}

GLboolean glIsTexture(GLuint texture) {
  return Context::GetCurrentContext().glIsTexture(texture);
}

void glLineWidth(GLfloat width) {
  Context::GetCurrentContext().glLineWidth(width);
}

void glLinkProgram(GLuint program) {
  Context::GetCurrentContext().glLinkProgram(program);
}

void glPixelStorei(GLenum pname, GLint param) {
  Context::GetCurrentContext().glPixelStorei(pname, param);
}

void glPolygonOffset(GLfloat factor, GLfloat units) {
  Context::GetCurrentContext().glPolygonOffset(factor, units);
}

void glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height,
                  GLenum format, GLenum type, GLvoid *pixels) {
  Context::GetCurrentContext().glReadPixels(x, y, width, height, format, type,
                                            pixels);
}

void glReleaseShaderCompiler(void) {
  Context::GetCurrentContext().glReleaseShaderCompiler();
}

void glRenderbufferStorage(GLenum target, GLenum internalformat, GLsizei width,
                           GLsizei height) {
  Context::GetCurrentContext().glRenderbufferStorage(target, internalformat,
                                                     width, height);
}

void glSampleCoverage(GLclampf value, GLboolean invert) {
  Context::GetCurrentContext().glSampleCoverage(value, invert);
}

void glScissor(GLint x, GLint y, GLsizei width, GLsizei height) {
  Context::GetCurrentContext().glScissor(x, y, width, height);
}

void glShaderBinary(GLsizei n, const GLuint *shaders, GLenum binaryformat,
                    const GLvoid *binary, GLsizei length) {
  Context::GetCurrentContext().glShaderBinary(n, shaders, binaryformat, binary,
                                              length);
}

void glShaderSource(GLuint shader, GLsizei count, const GLchar **string,
                    const GLint *length) {
  Context::GetCurrentContext().glShaderSource(shader, count, string, length);
}

void glStencilFunc(GLenum func, GLint ref, GLuint mask) {
  Context::GetCurrentContext().glStencilFunc(func, ref, mask);
}

void glStencilFuncSeparate(GLenum face, GLenum func, GLint ref, GLuint mask) {
  Context::GetCurrentContext().glStencilFuncSeparate(face, func, ref, mask);
}

void glStencilMask(GLuint mask) {
  Context::GetCurrentContext().glStencilMask(mask);
}

void glStencilMaskSeparate(GLenum face, GLuint mask) {
  Context::GetCurrentContext().glStencilMaskSeparate(face, mask);
}

void glStencilOp(GLenum fail, GLenum zfail, GLenum zpass) {
  Context::GetCurrentContext().glStencilOp(fail, zfail, zpass);
}

void glStencilOpSeparate(GLenum face, GLenum fail, GLenum zfail, GLenum zpass) {
  Context::GetCurrentContext().glStencilOpSeparate(face, fail, zfail, zpass);
}

void glTexImage2D(GLenum target, GLint level, GLint internalformat,
                  GLsizei width, GLsizei height, GLint border, GLenum format,
                  GLenum type, const GLvoid *pixels) {
  Context::GetCurrentContext().glTexImage2D(target, level, internalformat,
                                            width, height, border, format, type,
                                            pixels);
}

void glTexImage3D(GLenum target, GLint level, GLint internalformat,
                  GLsizei width, GLsizei height, GLsizei depth, GLint border,
                  GLenum format, GLenum type, const GLvoid *pixels) {
  Context::GetCurrentContext().glTexImage3D(target, level, internalformat,
                                            width, height, depth, border,
                                            format, type, pixels);
}

void glTexParameterf(GLenum target, GLenum pname, GLfloat param) {
  Context::GetCurrentContext().glTexParameterf(target, pname, param);
}

void glTexParameterfv(GLenum target, GLenum pname, const GLfloat *params) {
  Context::GetCurrentContext().glTexParameterfv(target, pname, params);
}

void glTexParameteri(GLenum target, GLenum pname, GLint param) {
  Context::GetCurrentContext().glTexParameteri(target, pname, param);
}

void glTexParameteriv(GLenum target, GLenum pname, const GLint *params) {
  Context::GetCurrentContext().glTexParameteriv(target, pname, params);
}

void glTexSubImage2D(GLenum target, GLint level, GLint xoffset, GLint yoffset,
                     GLsizei width, GLsizei height, GLenum format, GLenum type,
                     const GLvoid *pixels) {
  Context::GetCurrentContext().glTexSubImage2D(
      target, level, xoffset, yoffset, width, height, format, type, pixels);
}

void glTexSubImage3D(GLenum target, GLint level, GLint xoffset, GLint yoffset,
                     GLint zoffset, GLsizei width, GLsizei height,
                     GLsizei depth, GLenum format, GLenum type,
                     const GLvoid *pixels) {
  Context::GetCurrentContext().glTexSubImage3D(target, level, xoffset, yoffset,
                                               zoffset, width, height, depth,
                                               format, type, pixels);
}

void glUniform1f(GLint location, GLfloat x) {
  Context::GetCurrentContext().glUniform1f(location, x);
}

void glUniform1fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform1fv(location, count, v);
}

void glUniform1i(GLint location, GLint x) {
  Context::GetCurrentContext().glUniform1i(location, x);
}

void glUniform1iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform1iv(location, count, v);
}

void glUniform2f(GLint location, GLfloat x, GLfloat y) {
  Context::GetCurrentContext().glUniform2f(location, x, y);
}

void glUniform2fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform2fv(location, count, v);
}

void glUniform2i(GLint location, GLint x, GLint y) {
  Context::GetCurrentContext().glUniform2i(location, x, y);
}

void glUniform2iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform2iv(location, count, v);
}

void glUniform3f(GLint location, GLfloat x, GLfloat y, GLfloat z) {
  Context::GetCurrentContext().glUniform3f(location, x, y, z);
}

void glUniform3fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform3fv(location, count, v);
}

void glUniform3i(GLint location, GLint x, GLint y, GLint z) {
  Context::GetCurrentContext().glUniform3i(location, x, y, z);
}

void glUniform3iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform3iv(location, count, v);
}

void glUniform4f(GLint location, GLfloat x, GLfloat y, GLfloat z, GLfloat w) {
  Context::GetCurrentContext().glUniform4f(location, x, y, z, w);
}

void glUniform4fv(GLint location, GLsizei count, const GLfloat *v) {
  Context::GetCurrentContext().glUniform4fv(location, count, v);
}

void glUniform4i(GLint location, GLint x, GLint y, GLint z, GLint w) {
  Context::GetCurrentContext().glUniform4i(location, x, y, z, w);
}

void glUniform4iv(GLint location, GLsizei count, const GLint *v) {
  Context::GetCurrentContext().glUniform4iv(location, count, v);
}

void glUniformMatrix2fv(GLint location, GLsizei count, GLboolean transpose,
                        const GLfloat *value) {
  Context::GetCurrentContext().glUniformMatrix2fv(location, count, transpose,
                                                  value);
}

void glUniformMatrix3fv(GLint location, GLsizei count, GLboolean transpose,
                        const GLfloat *value) {
  Context::GetCurrentContext().glUniformMatrix3fv(location, count, transpose,
                                                  value);
}

void glUniformMatrix4fv(GLint location, GLsizei count, GLboolean transpose,
                        const GLfloat *value) {
  Context::GetCurrentContext().glUniformMatrix4fv(location, count, transpose,
                                                  value);
}

void glUseProgram(GLuint program) {
  Context::GetCurrentContext().glUseProgram(program);
}

void glValidateProgram(GLuint program) {
  Context::GetCurrentContext().glValidateProgram(program);
}

void glVertexAttrib1f(GLuint indx, GLfloat x) {
  Context::GetCurrentContext().glVertexAttrib1f(indx, x);
}

void glVertexAttrib1fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib1fv(indx, values);
}

void glVertexAttrib2f(GLuint indx, GLfloat x, GLfloat y) {
  Context::GetCurrentContext().glVertexAttrib2f(indx, x, y);
}

void glVertexAttrib2fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib2fv(indx, values);
}

void glVertexAttrib3f(GLuint indx, GLfloat x, GLfloat y, GLfloat z) {
  Context::GetCurrentContext().glVertexAttrib3f(indx, x, y, z);
}

void glVertexAttrib3fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib3fv(indx, values);
}

void glVertexAttrib4f(GLuint indx, GLfloat x, GLfloat y, GLfloat z, GLfloat w) {
  Context::GetCurrentContext().glVertexAttrib4f(indx, x, y, z, w);
}

void glVertexAttrib4fv(GLuint indx, const GLfloat *values) {
  Context::GetCurrentContext().glVertexAttrib4fv(indx, values);
}

void glVertexAttribPointer(GLuint indx, GLint size, GLenum type,
                           GLboolean normalized, GLsizei stride,
                           const GLvoid *ptr) {
  Context::GetCurrentContext().glVertexAttribPointer(indx, size, type,
                                                     normalized, stride, ptr);
}

void glViewport(GLint x, GLint y, GLsizei width, GLsizei height) {
  Context::GetCurrentContext().glViewport(x, y, width, height);
}

//
// LSGL extenstions.
//

// Set custom progress report function.
// LSGLProgressCallback is defined in gles_common.h
void lsglSetProgressCallback(LSGLProgressCallback func, void *userdata) {
  Context::GetCurrentContext().lsglSetProgressCallback(func, userdata);
}

void lsglSetPixelStep(GLint step) {
  Context::GetCurrentContext().lsglSetPixelStep(step);
}

GLint lsglGetPixelStep() {
  return Context::GetCurrentContext().lsglGetPixelStep();
}

// @fixme { Deprecated }
void lsglSetPointSize(GLfloat size) {
  Context::GetCurrentContext().lsglSetPointSize(size);
}

// @fixme { Deprecated }
void lsglSetPointSizev(GLsizei num, const GLfloat *size) {
  Context::GetCurrentContext().lsglSetPointSizev(num, size);
}

void lsglSetCamera(const GLfloat *eye, const GLfloat *target, const GLfloat *up,
                   GLfloat fov) {
  Context::GetCurrentContext().lsglSetCamera(eye, target, up, fov);
}

void lsglSetStereoEnvCamera(const GLfloat *eye, const GLfloat *target,
                            const GLfloat *up, GLfloat zeroParallax,
                            GLfloat eyeSeparation) {
  Context::GetCurrentContext().lsglSetStereoEnvCamera(
      eye, target, up, zeroParallax, eyeSeparation);
}

void lsglInvalidateBuffer(GLenum target) {
  Context::GetCurrentContext().lsglInvalidateBuffer(target);
}

void lsglEvalFragmentShader() {
  Context::GetCurrentContext().lsglEvalFragmentShader();
}

void lsglEvalSingleFragmentShader() {
  Context::GetCurrentContext().lsglEvalSingleFragmentShader();
}

void lsglSetShaderCompiler(const char *path, const char *options) {
  Context::GetCurrentContext().lsglSetShaderCompiler(path, options);
}

// Zero-copy extension
void lsglBufferDataPointer(GLenum target, GLsizeiptr size, const GLvoid *data,
                           GLenum usage) {
  Context::GetCurrentContext().lsglBufferDataPointer(target, size, data, usage);
}

void lsglTexImage3DPointer(GLenum target, GLint level, GLint internalformat,
                           GLsizei width, GLsizei height, GLsizei depth,
                           GLint border, GLenum format, GLenum type,
                           const GLvoid *pixels) {
  Context::GetCurrentContext().lsglTexImage3DPointer(
      target, level, internalformat, width, height, depth, border, format, type,
      pixels);
}

// Add coordinate remapping feature for non-uniform volume data
void lsglTexCoordRemap(GLenum target, GLenum coord, GLsizei size,
                       GLfloat *coords) {
  Context::GetCurrentContext().lsglTexCoordRemap(target, coord, size, coords);
}

// Alias for glTexPageCommitmentARB
void lsglTexPageCommitment(GLenum target, GLint level, GLint xoffset,
                           GLint yoffset, GLint zoffset, GLsizei width,
                           GLsizei height, GLsizei depth, GLboolean commit) {
  Context::GetCurrentContext().lsglTexPageCommitment(
      target, level, xoffset, yoffset, zoffset, width, height, depth, commit);
}

// Zero-copy version of glTexSubImage3D(for saving memory)
void lsglTexSubImage3DPointer(GLenum target, GLint level, GLint xoffset,
                              GLint yoffset, GLint zoffset, GLsizei width,
                              GLsizei height, GLsizei depth, GLenum format,
                              GLenum type, const GLvoid *pixels) {
  Context::GetCurrentContext().lsglTexSubImage3DPointer(
      target, level, xoffset, yoffset, zoffset, width, height, depth, format,
      type, pixels);
}
}

#endif // __LSGL_GLES_C_API_H__
