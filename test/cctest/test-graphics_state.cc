#include "gtest/gtest.h"
#include "gles_context.h"
#include "gles_resource_manager.h"

using namespace lsgl;

static int windowWidth = 512;
static int windowHeight = 512;

TEST(ClearTest, ClearRGBA8) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8_OES, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearColor(1.0f, 0.5f, 0.25f, 0.125f);
  glClear(GL_COLOR_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}

TEST(ClearTest, ClearRGBA32F) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32F_EXT, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearColor(1.0f, 0.5f, 0.25f, 0.125f);
  glClear(GL_COLOR_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(ReadPixelTest, ReadRGBA) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearColor(1.0f, 0.5f, 0.25f, 0.125f);
  glClear(GL_COLOR_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  std::vector<unsigned char> img(windowWidth * windowHeight * 4);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  EXPECT_EQ(255, img[0]);
  EXPECT_EQ(127, img[1]);
  EXPECT_EQ(63, img[2]);
  EXPECT_EQ(31, img[3]);
}

TEST(ReadPixelTest, ReadRGBA8) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8_OES, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearColor(1.0f, 0.5f, 0.25f, 0.125f);
  glClear(GL_COLOR_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  std::vector<unsigned char> img(windowWidth * windowHeight * 4);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  EXPECT_EQ(255, img[0]);
  EXPECT_EQ(127, img[1]);
  EXPECT_EQ(63, img[2]);
  EXPECT_EQ(31, img[3]);

}

TEST(ReadPixelTest, ReadRGBA32F) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA32F_EXT, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearColor(1.0f, 0.5f, 0.25f, 0.125f);
  glClear(GL_COLOR_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  std::vector<float> img(windowWidth * windowHeight * 4);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  EXPECT_FLOAT_EQ(1.0f, img[0]);
  EXPECT_FLOAT_EQ(0.5f, img[1]);
  EXPECT_FLOAT_EQ(0.25f, img[2]);
  EXPECT_FLOAT_EQ(0.125f, img[3]);
}

TEST(ReadPixelTest, ReadDEPTH) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearDepthf(0.25f);
  glClear(GL_DEPTH_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  std::vector<float> img(windowWidth * windowHeight);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_DEPTH_COMPONENT, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  EXPECT_FLOAT_EQ(0.25f, img[0]);
}

TEST(ReadPixelTest, ReadDEPTH32OES) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32_OES, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, colorRenderbuffer);

  glClearDepthf(0.25f);
  glClear(GL_DEPTH_BUFFER_BIT);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  std::vector<float> img(windowWidth * windowHeight);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_DEPTH_COMPONENT, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  EXPECT_FLOAT_EQ(0.25f, img[0]);
}

TEST(ReadPixelTest, ReadUnattachedReadBuf) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);

  std::vector<unsigned char> img(windowWidth * windowHeight);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, &img.at(0));
  EXPECT_EQ(GL_INVALID_FRAMEBUFFER_OPERATION, glGetError());
}

TEST(ReadPixelTest, ReadUnattachedDepthBuf) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  Context &ctx = Context::GetCurrentContext();

  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);

  std::vector<float> img(windowWidth * windowHeight);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_DEPTH_COMPONENT, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_INVALID_FRAMEBUFFER_OPERATION, glGetError());
}
