#include "gtest/gtest.h"
#include "gles_c_api.h" 

using namespace lsgl;

#include <cstdio>
#include <cstdlib>
#include <cassert>

static int windowWidth = 256;
static int windowHeight = 256;

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
										
										
#if defined(__sparc__)
// @todo
#else

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

TEST(DrawTest, Draw0) {

  glViewport(0, 0, windowWidth, windowHeight);

  GLuint prog = 0, fragShader = 0;
  bool ret = CompileShader(prog, fragShader, fragShaderCode);
  assert(ret);

  glUseProgram(prog);

  // update shader vertex attribute indices
  GLint attrNormal = glGetAttribLocation(prog, "normal");
  printf("attrNormal = %d\n", attrNormal);
  GLint attrPos = glGetAttribLocation(prog, "position");
  printf("attrPos = %d\n", attrPos);
  EXPECT_EQ(0, attrPos);

  printf("GenFramebuffer\n");
  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8_OES, windowWidth, windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, colorRenderbuffer);

  // create depth renderbuffer and attach
  GLuint depthRenderbuffer;
  glGenRenderbuffers(1, &depthRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, depthRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32_OES, windowWidth, windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthRenderbuffer);


  glClearColor(1, 1, 1,1 );
  glClear(GL_COLOR_BUFFER_BIT);

  const float cube[] = {  -1.0, -1.0, -1.0,
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
  glGenBuffers(1, &cubevtx);
  glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);

  glGenBuffers(1, &cubeidx);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeidx);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);

  float eye[3] = {2.0f, -2.0f, 5.0f}; 
  float lookat[3] = {0.0f, 0.0f, 0.0f}; 
  float up[3] = {0.0f, 1.0f, 0.0f}; 
  lsglSetCamera(eye, lookat, up, 45.0f);

  // enable subsampling
  glEnable(GL_SAMPLE_COVERAGE);
  glSampleCoverage(1.0f, GL_FALSE);

  //GLfloat scalefactor = 1.0f;
  //glUniform1f(glGetUniformLocation(prog, "scalefactor"), scalefactor);

  GLfloat resolution[2];
  resolution[0] = (float)windowWidth;
  resolution[1] = (float)windowHeight;
  glUniform2fv(glGetUniformLocation(prog, "resolution"), 1, resolution);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  // Use vertex buffers.
  glBindBuffer(GL_ARRAY_BUFFER, cubevtx);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeidx);
  glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  glVertexAttribPointer(attrNormal, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  glEnableVertexAttribArray(attrPos);
  glEnableVertexAttribArray(attrNormal);
  glDrawElements(GL_TRIANGLES, sizeof(cube_indices) / sizeof(GLuint), GL_UNSIGNED_INT, (void*)0); // 0 = use vertex buffer.
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  //glFinish();
  //EXPECT_EQ(GL_NO_ERROR, glGetError());

  //glClearColor(1, 1, 1,1 );
  //glClear(GL_COLOR_BUFFER_BIT);

  Context& ctx = Context::GetCurrentContext();
  EXPECT_EQ(1, ctx.state_.currentDrawStackIndex);

  glDrawElements(GL_TRIANGLES, sizeof(cube_indices) / sizeof(GLuint), GL_UNSIGNED_INT, (void*)0); // 0 = use vertex buffer.
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  glFinish();
  EXPECT_EQ(GL_NO_ERROR, glGetError());

  glDeleteRenderbuffers(1, &colorRenderbuffer);
  glDeleteRenderbuffers(1, &depthRenderbuffer);
  glDeleteFramebuffers(1, &framebuffer);
  EXPECT_EQ(GL_NO_ERROR, glGetError());

}
#endif
