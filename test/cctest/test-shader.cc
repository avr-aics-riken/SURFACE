#include "gtest/gtest.h"
#include "gles_context.h"

using namespace lsgl;

static const char* fragShaderCode0 = 
"#ifdef GL_ES\n"
"precision mediump float;\n"
"#endif\n"
"\n"
"uniform sampler2D tex0;\n"
"uniform vec2      resolution;\n"
"varying vec3      normal;\n"
"\n"
"void main( void ) {\n"
"     vec3 color = vec3(0);\n"
"     const int MAXITER = 30;\n"
"     for (int i = 0; i < MAXITER; i++) {\n"
"         color += 0.1*(float(i)/float(MAXITER));\n"
"     }\n"
"     gl_FragColor = vec4(gl_FragCoord.xyz / vec3(resolution, 1.0), 1.0);\n"
"}\n";

static const char* fragShaderCode1 = 
"#ifdef GL_ES\n"
"precision mediump float;\n"
"#endif\n"
"\n"
"uniform sampler2D tex0;\n"
"uniform vec2      resolution;\n"
"varying vec3      normal;\n"
"\n"
"void main( void ) {\n"
"     vec3 color = vec3(0);\n"
"     const int MAXITER = 30;\n"
"     for (int i = 0; i < MAXITER; i++) {\n"
"         color += 0.1*(float(i)/float(MAXITER));\n"
"     }\n"
"     gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);\n"
"}\n";

TEST(ShaderTest, UseProgram0) {

  Context &ctx = Context::GetCurrentContext();

  ctx.glUseProgram(0);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

#if defined(__sparc__)
// @todo
#else
TEST(ShaderTest, CompileShader) {

  Context &ctx = Context::GetCurrentContext();

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  const GLchar *src = fragShaderCode0;
    
  GLuint fragShader = ctx.glCreateShader(GL_FRAGMENT_SHADER);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  ctx.glShaderSource(fragShader, 1, &src, NULL);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
  ctx.glCompileShader(fragShader);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
  int val;
  ctx.glGetShaderiv(fragShader, GL_COMPILE_STATUS, &val);
  EXPECT_EQ(GL_TRUE, val);
}

TEST(ShaderTest, CompileShader2) {

  Context &ctx = Context::GetCurrentContext();

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  const GLchar *src = fragShaderCode0;
    
  GLuint fragShader0 = ctx.glCreateShader(GL_FRAGMENT_SHADER);
  EXPECT_EQ(GL_NO_ERROR, glGetError());


  ctx.glShaderSource(fragShader0, 1, &src, NULL);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
  ctx.glCompileShader(fragShader0);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
  int val;
  ctx.glGetShaderiv(fragShader0, GL_COMPILE_STATUS, &val);
  EXPECT_EQ(GL_TRUE, val);

  GLuint fragShader1 = ctx.glCreateShader(GL_FRAGMENT_SHADER);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  ctx.glShaderSource(fragShader1, 1, &src, NULL);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
  ctx.glCompileShader(fragShader1);
  EXPECT_EQ(GL_NO_ERROR, glGetError());
  ctx.glGetShaderiv(fragShader1, GL_COMPILE_STATUS, &val);
  EXPECT_EQ(GL_TRUE, val);

}
#endif

TEST(ShaderTest, ShaderSource) {

  Context &ctx = Context::GetCurrentContext();

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  const GLchar *src = fragShaderCode0;
    
  GLuint sid = -1;
  ctx.glShaderSource(sid, 1, &src, NULL);
  EXPECT_EQ(GL_INVALID_VALUE, glGetError());

  GLuint fragShader0 = ctx.glCreateShader(GL_FRAGMENT_SHADER);
  ctx.glShaderSource(fragShader0, -1, &src, NULL);
  EXPECT_EQ(GL_INVALID_VALUE, glGetError());

  ctx.glFinish();

}
