#include "gtest/gtest.h"
#include "gles_context.h"

using namespace lsgl;

TEST(TextureTest, TextureLuminanceByte) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  int w = 512;
  int h = 512;

  std::vector<char> img(w*h);

  GLuint tex;
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, w, h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(TextureTest, TextureRGBAByte) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  int w = 512;
  int h = 512;

  std::vector<char> img(w*h*4);

  GLuint tex;
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(TextureTest, TextureLumianceFloat) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  int w = 512;
  int h = 512;

  std::vector<float> img(w*h);

  GLuint tex;
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, w, h, 0, GL_LUMINANCE, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(TextureTest, TextureRGBAFloat) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  int w = 512;
  int h = 512;

  std::vector<float> img(w*h*4);

  GLuint tex;
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(TextureTest, Texture3DRGBByte) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  int w = 32;
  int h = 32;
  int d = 32;

  std::vector<unsigned char> img(w*h*d*3);

  GLuint tex;
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_3D, tex);
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, w, h, d, 0, GL_RGB, GL_UNSIGNED_BYTE, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(TextureTest, Texture3DRGBFloat) {

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  int w = 32;
  int h = 32;
  int d = 32;

  std::vector<float> img(w*h*d*3);

  GLuint tex;
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_3D, tex);
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB, w, h, d, 0, GL_RGB, GL_FLOAT, &img.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}
