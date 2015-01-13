#include "gtest/gtest.h"
#include "gles_context.h"
#include "gles_resource_manager.h"

using namespace lsgl;

TEST(BufferTest, GenNegativeNum) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      
  GLuint buffer;
  ctx.glGenBuffers(-1, &buffer);
  EXPECT_EQ(GL_INVALID_VALUE, glGetError()); 
}

TEST(BufferTest, GenZeroNum) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      
  GLuint buffer;
  ctx.glGenBuffers(0, &buffer);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

TEST(BufferTest, Gen) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      
  GLuint buffer;
  ctx.glGenBuffers(1, &buffer);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  ctx.glDeleteBuffers(1, &buffer);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

TEST(BufferTest, DelNegative) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      
  GLuint buffer;
  ctx.glDeleteBuffers(-1, &buffer);
  EXPECT_EQ(GL_INVALID_VALUE, glGetError()); 
}

TEST(BufferTest, DelNull) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      
  EXPECT_DEATH(ctx.glDeleteBuffers(1, NULL), "buffers are null");
}

TEST(RenderBufferTest, RGBA8) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
  Context& ctx = Context::GetCurrentContext();      
  GLuint colorRenderbuffer;
  ctx.glGenRenderbuffers(1, &colorRenderbuffer);
  ctx.glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  ctx.glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8_OES, 128, 128);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  const Renderbuffer* buf = ctx.resourceManager_.GetRenderbuffer(colorRenderbuffer);
  EXPECT_TRUE(buf);

  EXPECT_EQ(128, buf->GetWidth());
  EXPECT_EQ(128, buf->GetHeight());
  EXPECT_EQ(4, buf->GetBytesPerPixel());
}

TEST(RenderBufferTest, RGBAF32) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
  Context& ctx = Context::GetCurrentContext();      
  GLuint colorRenderbuffer;
  ctx.glGenRenderbuffers(1, &colorRenderbuffer);
  ctx.glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  ctx.glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32F_EXT, 128, 128);
  ASSERT_EQ(GL_NO_ERROR, glGetError());

  const Renderbuffer* buf = ctx.resourceManager_.GetRenderbuffer(colorRenderbuffer);
  EXPECT_TRUE(buf);

  EXPECT_EQ(128, buf->GetWidth());
  EXPECT_EQ(128, buf->GetHeight());
  EXPECT_EQ(16, buf->GetBytesPerPixel());
}

TEST(RenderBufferTest, DEPTH) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
  Context& ctx = Context::GetCurrentContext();      
  GLuint colorRenderbuffer;
  ctx.glGenRenderbuffers(1, &colorRenderbuffer);
  ctx.glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  ctx.glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, 128, 128);
  ASSERT_EQ(GL_NO_ERROR, glGetError());

  const Renderbuffer* buf = ctx.resourceManager_.GetRenderbuffer(colorRenderbuffer);
  EXPECT_TRUE(buf);

  EXPECT_EQ(128, buf->GetWidth());
  EXPECT_EQ(128, buf->GetHeight());
  EXPECT_EQ(4, buf->GetBytesPerPixel()); // float x 1 channel
}

TEST(RenderBufferTest, DEPTH32) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
  Context& ctx = Context::GetCurrentContext();      
  GLuint colorRenderbuffer;
  ctx.glGenRenderbuffers(1, &colorRenderbuffer);
  ctx.glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  ctx.glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32_OES, 128, 128);
  ASSERT_EQ(GL_NO_ERROR, glGetError());

  const Renderbuffer* buf = ctx.resourceManager_.GetRenderbuffer(colorRenderbuffer);
  EXPECT_TRUE(buf);

  EXPECT_EQ(128, buf->GetWidth());
  EXPECT_EQ(128, buf->GetHeight());
  EXPECT_EQ(4, buf->GetBytesPerPixel()); // float x 1 channel
}

TEST(RenderBufferTest, DEPTH16) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
  Context& ctx = Context::GetCurrentContext();      
  GLuint colorRenderbuffer;
  ctx.glGenRenderbuffers(1, &colorRenderbuffer);
  ctx.glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  ctx.glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, 128, 128);
  // 16bit depth is not supported.
  ASSERT_EQ(GL_INVALID_ENUM, glGetError());
}

TEST(RenderBufferTest, Paraameteriv) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
  Context& ctx = Context::GetCurrentContext();      
  GLuint colorRenderbuffer;
  ctx.glGenRenderbuffers(1, &colorRenderbuffer);
  ctx.glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  ctx.glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8_OES, 128, 128);
  ASSERT_EQ(GL_NO_ERROR, glGetError());

  GLint w, h;
  ctx.glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_WIDTH, &w);
  ASSERT_EQ(GL_NO_ERROR, glGetError());

  ctx.glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_HEIGHT, &h);
  ASSERT_EQ(GL_NO_ERROR, glGetError());

  EXPECT_EQ(128, w);
  EXPECT_EQ(128, h);
}
