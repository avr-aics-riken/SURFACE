#include "gtest/gtest.h"
#include "gles_context.h"

using namespace lsgl;

static const char* fragShaderCode = 
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
										
										

static bool
CompileShader(
  GLuint& prog,
  GLuint& fragShader,
  const char* fragShaderSource)
{
  GLint val = 0;

  // free old shader/program
  if (prog != 0)   glDeleteProgram(prog);
  if (fragShader != 0) glDeleteShader(fragShader);

  static const GLchar *src = fragShaderSource;
    
  fragShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragShader, 1, &src, NULL);
  glCompileShader(fragShader);
  glGetShaderiv(fragShader, GL_COMPILE_STATUS, &val);
  assert(val == GL_TRUE && "failed to compile shader");

  prog = glCreateProgram();
  glAttachShader(prog, fragShader);
  glLinkProgram(prog);
  glGetProgramiv(prog, GL_LINK_STATUS, &val);
  assert(val == GL_TRUE && "failed to link shader");

  return true;
}

#if defined(__sparc__)
// @todo
#else

TEST(LSGLExtTest, UniformLSGL) {

  // lsgl-ext uniforms
  //
  // lsgl_World : 4x4 float
  // lsgl_View  : 4x4 float
  // lsgl_Proj  : 4x4 float

  EXPECT_EQ(GL_NO_ERROR, glGetError());

  GLuint prog = 0, fragShader = 0;
  bool ret = CompileShader(prog, fragShader, fragShaderCode);
  assert(ret);

  glUseProgram(prog);


  int loc = glGetUniformLocation(prog, "lsgl_World");
  EXPECT_NE(-1, loc);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  float m[4][4];
  float n[4][4];

  glGetUniformfv(prog, loc, (GLfloat*)n);
  EXPECT_FLOAT_EQ(1.0f, n[0][0]);
  EXPECT_FLOAT_EQ(0.0f, n[1][0]);
  EXPECT_FLOAT_EQ(0.0f, n[2][0]);
  EXPECT_FLOAT_EQ(0.0f, n[3][0]);
  EXPECT_FLOAT_EQ(1.0f, n[1][1]);
  EXPECT_FLOAT_EQ(1.0f, n[2][2]);
  EXPECT_FLOAT_EQ(1.0f, n[3][3]);

  m[0][0] = 1.0f;
  m[1][0] = 2.0f;
  m[2][0] = 0.0f;
  glUniformMatrix4fv(loc, 1, 0, (const GLfloat*)m);

  glGetUniformfv(prog, loc, (GLfloat*)n);
  EXPECT_FLOAT_EQ(1.0f, n[0][0]);
  EXPECT_FLOAT_EQ(2.0f, n[1][0]);
  EXPECT_FLOAT_EQ(0.0f, n[2][0]);

  loc = glGetUniformLocation(prog, "lsgl_View");
  EXPECT_NE(-1, loc);
  m[0][0] = 1.1f;
  m[1][0] = 2.1f;
  m[2][0] = 0.1f;
  glUniformMatrix4fv(loc, 1, 0, (const GLfloat*)m);
  glGetUniformfv(prog, loc, (GLfloat*)n);
  EXPECT_FLOAT_EQ(1.1f, n[0][0]);
  EXPECT_FLOAT_EQ(2.1f, n[1][0]);
  EXPECT_FLOAT_EQ(0.1f, n[2][0]);

  loc = glGetUniformLocation(prog, "lsgl_Proj");
  EXPECT_NE(-1, loc);

  m[0][0] = 1.2f;
  m[1][0] = 2.2f;
  m[2][0] = 0.2f;
  glUniformMatrix4fv(loc, 1, 0, (const GLfloat*)m);
  glGetUniformfv(prog, loc, (GLfloat*)n);
  EXPECT_FLOAT_EQ(1.2f, n[0][0]);
  EXPECT_FLOAT_EQ(2.2f, n[1][0]);
  EXPECT_FLOAT_EQ(0.2f, n[2][0]);

}
#endif

TEST(LSGLExtTest, CompilerPath) {

  Context &ctx = Context::GetCurrentContext();

  char* org_variable = strdup("");
  const char* org_compiler_path = getenv("GLSL_COMPILER");
  if (org_compiler_path) {
    org_variable = strdup(org_compiler_path);
  }

  setenv("GLSL_COMPILER", "", 1);
  const char* tmpval = getenv("GLSL_COMPILER");
  EXPECT_STREQ("", tmpval);

  const char* path = "/my/bin/glslc";
  const char* opts = "-O2";
  ctx.lsglSetShaderCompiler(path, opts);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  const char* compiler_path = getenv("GLSL_COMPILER");
  EXPECT_STREQ(path, compiler_path);

  setenv("GLSL_COMPILER", org_compiler_path, 1);

  free(org_variable);
}

TEST(LSGLExtTest, TexImage3DPointerScalarByte) {

  Context &ctx = Context::GetCurrentContext();

  std::vector<unsigned char> data(4*4*4);
  GLuint texid;
  ctx.glGenTextures(1, &texid);
  ctx.glActiveTexture(GL_TEXTURE0);
  ctx.glBindTexture(GL_TEXTURE_3D, texid);
  ctx.lsglTexImage3DPointer(GL_TEXTURE_3D, 0, GL_LUMINANCE, 4, 4, 4, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, &data.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());
}

TEST(LSGLExtTest, TexImage3DPointerScalar) {

  Context &ctx = Context::GetCurrentContext();

  std::vector<float> data(4*4*4);
  GLuint texid;
  ctx.glGenTextures(1, &texid);
  ctx.glActiveTexture(GL_TEXTURE0);
  ctx.glBindTexture(GL_TEXTURE_3D, texid);
  ctx.lsglTexImage3DPointer(GL_TEXTURE_3D, 0, GL_LUMINANCE, 4, 4, 4, 0, GL_LUMINANCE, GL_FLOAT, &data.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}

TEST(LSGLExtTest, TexImage3DPointerScalarDouble) {

  Context &ctx = Context::GetCurrentContext();

  std::vector<double> data(4*4*4);
  GLuint texid;
  ctx.glGenTextures(1, &texid);
  ctx.glActiveTexture(GL_TEXTURE0);
  ctx.glBindTexture(GL_TEXTURE_3D, texid);
  ctx.lsglTexImage3DPointer(GL_TEXTURE_3D, 0, GL_LUMINANCE, 4, 4, 4, 0, GL_LUMINANCE, GL_DOUBLE, &data.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}

TEST(LSGLExtTest, TexImage3DPointerVectorByte) {

  Context &ctx = Context::GetCurrentContext();

  std::vector<unsigned char> data(4*4*4*3);
  GLuint texid;
  ctx.glGenTextures(1, &texid);
  ctx.glActiveTexture(GL_TEXTURE0);
  ctx.glBindTexture(GL_TEXTURE_3D, texid);
  ctx.lsglTexImage3DPointer(GL_TEXTURE_3D, 0, GL_RGB, 4, 4, 4, 0, GL_RGB, GL_UNSIGNED_BYTE, &data.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}

TEST(LSGLExtTest, TexImage3DPointerVector) {

  Context &ctx = Context::GetCurrentContext();

  std::vector<float> data(4*4*4*3);
  GLuint texid;
  ctx.glGenTextures(1, &texid);
  ctx.glActiveTexture(GL_TEXTURE0);
  ctx.glBindTexture(GL_TEXTURE_3D, texid);
  ctx.lsglTexImage3DPointer(GL_TEXTURE_3D, 0, GL_RGB, 4, 4, 4, 0, GL_RGB, GL_FLOAT, &data.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}

TEST(LSGLExtTest, TexImage3DPointerVectorDouble) {

  Context &ctx = Context::GetCurrentContext();

  std::vector<double> data(4*4*4*3);
  GLuint texid;
  ctx.glGenTextures(1, &texid);
  ctx.glActiveTexture(GL_TEXTURE0);
  ctx.glBindTexture(GL_TEXTURE_3D, texid);
  ctx.lsglTexImage3DPointer(GL_TEXTURE_3D, 0, GL_RGB, 4, 4, 4, 0, GL_RGB, GL_DOUBLE, &data.at(0));
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}

TEST(VertexAttribTest, VertexAttribDouble) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  const double cube[] = {  -1.0, -1.0, -1.0,
                          -1.0, -1.0,  1.0,
                          -1.0,  1.0, -1.0,
                          -1.0,  1.0,  1.0,
                           1.0, -1.0, -1.0,
                           1.0, -1.0,  1.0,
                           1.0,  1.0, -1.0,
                           1.0,  1.0,  1.0 };

  const GLuint cube_indices[] = {
      0, 2, 1, 1, 2, 3,
      4, 5, 6, 5, 7, 6,
      0, 1, 4, 1, 5, 4,
      2, 7, 3, 2, 6, 7,
      1, 3, 5, 3, 7, 5,
      0, 4, 2, 2, 4, 6 };

  // Create Vertex Buffers.
  GLuint cubevtx, cubeidx;
  ctx.glGenBuffers(1, &cubevtx);
  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  ctx.glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  ctx.glGenBuffers(1, &cubeidx);
  ctx.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeidx);
  ctx.glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeidx);
  GLuint pos = 0;
  glVertexAttribPointer(pos, 3, GL_DOUBLE, GL_FALSE, sizeof(double) * 3, (void*)0);
  glEnableVertexAttribArray(pos);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

TEST(BufferTest, BufferDataPointer) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  const double cube[] = {  -1.0, -1.0, -1.0,
                          -1.0, -1.0,  1.0,
                          -1.0,  1.0, -1.0,
                          -1.0,  1.0,  1.0,
                           1.0, -1.0, -1.0,
                           1.0, -1.0,  1.0,
                           1.0,  1.0, -1.0,
                           1.0,  1.0,  1.0 };

  const GLuint cube_indices[] = {
      0, 2, 1, 1, 2, 3,
      4, 5, 6, 5, 7, 6,
      0, 1, 4, 1, 5, 4,
      2, 7, 3, 2, 6, 7,
      1, 3, 5, 3, 7, 5,
      0, 4, 2, 2, 4, 6 };

  // Create Vertex Buffers.
  GLuint cubevtx, cubeidx;
  ctx.glGenBuffers(1, &cubevtx);
  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  ctx.lsglBufferDataPointer(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

TEST(LSGLExtTest, BufferDataPointer2) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  const double cube[] = {  -1.0, -1.0, -1.0,
                          -1.0, -1.0,  1.0,
                          -1.0,  1.0, -1.0,
                          -1.0,  1.0,  1.0,
                           1.0, -1.0, -1.0,
                           1.0, -1.0,  1.0,
                           1.0,  1.0, -1.0,
                           1.0,  1.0,  1.0 };

  const GLuint cube_indices[] = {
      0, 2, 1, 1, 2, 3,
      4, 5, 6, 5, 7, 6,
      0, 1, 4, 1, 5, 4,
      2, 7, 3, 2, 6, 7,
      1, 3, 5, 3, 7, 5,
      0, 4, 2, 2, 4, 6 };

  // Create Vertex Buffers.
  GLuint cubevtx, cubevtx2, cubeidx;
  ctx.glGenBuffers(1, &cubevtx);
  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  ctx.lsglBufferDataPointer(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
  ctx.glGenBuffers(1, &cubevtx2);
  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx2);
  ctx.lsglBufferDataPointer(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 
}

TEST(LSGLExtTest, BufferDataPointerPoint) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  const float cube[] = {  -1.0, -1.0, -1.0,
                          -1.0, -1.0,  1.0,
                          -1.0,  1.0, -1.0,
                          -1.0,  1.0,  1.0,
                           1.0, -1.0, -1.0,
                           1.0, -1.0,  1.0,
                           1.0,  1.0, -1.0,
                           1.0,  1.0,  1.0 };

  GLuint prog = 0, fragShader = 0;
  bool ret = CompileShader(prog, fragShader, fragShaderCode);
  assert(ret);

  ctx.glUseProgram(prog);

  // update shader vertex attribute indices
  GLint attrPos = glGetAttribLocation(prog, "position");
  printf("attr = %d\n", attrPos);

  //ctx.lsglSetPointSize(1.0);

  // Create Vertex Buffers.
  GLuint cubevtx;
  ctx.glGenBuffers(1, &cubevtx);
  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  ctx.lsglBufferDataPointer(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);

  ctx.glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  // 0 = position
  ctx.glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  ctx.glEnableVertexAttribArray(attrPos);
  assert(glGetError() == GL_NO_ERROR);
}

TEST(LSGLExtTest, PixelStepTest) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  ctx.lsglSetPixelStep(4);
  assert(glGetError() == GL_NO_ERROR);

  int step = ctx.lsglGetPixelStep();
  EXPECT_EQ(4, step);

}

TEST(LSGLExtTest, LSGLPointSizeUniform) {
    
  EXPECT_EQ(GL_NO_ERROR, glGetError()); 

  Context& ctx = Context::GetCurrentContext();      

  const float cube[] = {  -1.0, -1.0, -1.0,
                          -1.0, -1.0,  1.0,
                          -1.0,  1.0, -1.0,
                          -1.0,  1.0,  1.0,
                           1.0, -1.0, -1.0,
                           1.0, -1.0,  1.0,
                           1.0,  1.0, -1.0,
                           1.0,  1.0,  1.0 };

  GLuint prog = 0, fragShader = 0;
  bool ret = CompileShader(prog, fragShader, fragShaderCode);
  assert(ret);

  ctx.glUseProgram(prog);

  // "lsgl_PointSize" is predefined lsgl ext variable.
  GLint psizePos = glGetUniformLocation(prog, "lsgl_PointSize");
  EXPECT_NE((GLuint)(-1), psizePos);
  assert(glGetError() == GL_NO_ERROR);

  float psize;
  glGetUniformfv(prog, psizePos, &psize);
  EXPECT_FLOAT_EQ(1.0f, psize);

  assert(glGetError() == GL_NO_ERROR);

}
