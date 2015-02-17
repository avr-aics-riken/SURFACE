#include "../../gles/gles_c_api.h"
#include <GLES2/gl2ext.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <string>

#include "../common/SimpleTGA.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

int windowWidth = 512;
int windowHeight = 512;

// xorshift
static float randomfloat(void) {
  static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
  unsigned int t=x^x<<11; x=y; y=z; z=w; w^=w>>19^t^t>>8;
  return (float)(w*(1.0/4294967296.0));
}

void GenerateRandomLines(
  float* vertices, float* width, size_t nLines, const float bmin[3], const float bmax[3])
{

  assert(vertices);

  for (size_t i = 0; i < 2*nLines; i++) {

    // [0, 1.0)
    float u0 = randomfloat();
    float u1 = randomfloat();
    float u2 = randomfloat();
    float w  = 10.0f * randomfloat();

    // [bmin, bmax)
    float px = (bmax[0] - bmin[0]) * u0 + bmin[0];
    float py = (bmax[1] - bmin[1]) * u1 + bmin[1];
    float pz = (bmax[2] - bmin[2]) * u2 + bmin[2];

    vertices[3*i+0] = px;
    vertices[3*i+1] = py;
    vertices[3*i+2] = pz;

    width[i] = w;

  }
  
}


bool SaveColorBufferRGBA(const char *savename) {
  void *tgabuffer;
  unsigned char *imgBuf = new unsigned char[windowWidth * windowHeight * 4];
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE,
               imgBuf);
  int tgasize =
      SimpleTGASaverRGBA(&tgabuffer, windowWidth, windowHeight, imgBuf);
  delete[] imgBuf;
  if (!tgasize) {
    printf("Failed save.\n");
    return false;
  }

  FILE *fp = fopen(savename, "wb");
  fwrite(tgabuffer, 1, tgasize, fp);
  fclose(fp);
  free(tgabuffer);

  return true;
}

static bool LoadShader(GLuint &prog, GLuint &fragShader,
                       const char *fragShaderFilename) {
  GLint val = 0;

  // free old shader/program
  if (prog != 0)
    glDeleteProgram(prog);
  if (fragShader != 0)
    glDeleteShader(fragShader);

  static GLchar srcbuf[16384];
  FILE *fp = fopen(fragShaderFilename, "rb");
  if (!fp) {
    return false;
  }

  fseek(fp, 0, SEEK_END);
  size_t len = ftell(fp);
  rewind(fp);
  len = fread(srcbuf, 1, len, fp);
  srcbuf[len] = 0;
  fclose(fp);

  const GLchar *src = srcbuf;

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

static bool LoadBinaryShader(GLuint &prog, GLuint &fragShader,
                             const char *fragShaderBinaryFilename) {
  GLint val = 0;

  // free old shader/program
  if (prog != 0)
    glDeleteProgram(prog);
  if (fragShader != 0)
    glDeleteShader(fragShader);

  std::vector<unsigned char> data;
  FILE *fp = fopen(fragShaderBinaryFilename, "rb");
  if (!fp) {
    fprintf(stderr, "Failed to open file: %s\n", fragShaderBinaryFilename);
    return false;
  }
  fseek(fp, 0, SEEK_END);
  size_t len = ftell(fp);
  rewind(fp);

  data.resize(len);
  len = fread(&data.at(0), 1, len, fp);
  fclose(fp);

  const void *binary = &data.at(0);
  GLenum format = 0; // format is arbitrary at this time.

  fragShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderBinary(1, &fragShader, format, binary, len);

  prog = glCreateProgram();
  glAttachShader(prog, fragShader);
  glLinkProgram(prog);

  glGetShaderiv(fragShader, GL_COMPILE_STATUS, &val);
  assert(val == GL_TRUE && "failed to compile shader");

  return true;
}

int main(int argc, char **argv) {
#ifdef ENABLE_MPI
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  const char *fragShaderFile = "input.frag";
  const char *fragShaderBinFile = "shader.so";
  int dim[3] = {32, 32, 32};

  int numLines = 100;

  if (argc > 1) {
    numLines = atoi(argv[1]);
  }

  std::vector<float> positions(2 * 3 * numLines);
  std::vector<float> widths(2 * numLines);

  float bmin[3] = {-100, -100, -100};
  float bmax[3] = { 100,  100,  100};
  GenerateRandomLines(&positions.at(0), &widths.at(0), numLines, bmin, bmax);

  GLuint prog = 0, fragShader = 0;
  bool ret = LoadShader(prog, fragShader, fragShaderFile);
  if (!ret) {
    fprintf(stderr, "failed to load shader: %s\n", fragShaderFile);
    exit(-1);
  }

  glUseProgram(prog);

  // update shader vertex attribute indices
  GLint attrPos = glGetAttribLocation(prog, "position");
  printf("attr = %d\n", attrPos);

  GLint matPos = glGetAttribLocation(prog, "matID");
  printf("matID = %d\n", matPos);

  printf("GenFramebuffer\n");
  GLuint framebuffer;
  glGenFramebuffers(1, &framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

  // create color renderbuffer and attach
  GLuint colorRenderbuffer;
  glGenRenderbuffers(1, &colorRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, windowWidth, windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                            GL_RENDERBUFFER, colorRenderbuffer);

  // create depth renderbuffer and attach
  GLuint depthRenderbuffer;
  glGenRenderbuffers(1, &depthRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, depthRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, windowWidth,
                        windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, depthRenderbuffer);

  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  glViewport(0, 0, windowWidth, windowHeight);

  GLboolean ok = glIsProgram(prog);
  assert(ok == GL_TRUE);

  // false = dont' cap line primitive at the extent.
  //glUniform1i(glGetUniformLocation(prog, "lsgl_LineCap"), 0);

  // 1. Create Vertex Buffers.
  GLuint lnvtx, lradius;

  glGenBuffers(1, &lnvtx);
  glBindBuffer(GL_ARRAY_BUFFER, lnvtx);
  lsglBufferDataPointer(GL_ARRAY_BUFFER, numLines * sizeof(float) * 3 * 2, &positions.at(0), GL_STATIC_DRAW);

  glGenBuffers(1, &lradius);
  glBindBuffer(GL_ARRAY_BUFFER, lradius);
  lsglBufferDataPointer(GL_ARRAY_BUFFER, numLines * sizeof(float) * 2, &widths.at(0), GL_STATIC_DRAW);

  // Set varying line with using "lsgl_LineWidth" extenstion.
  GLint lsizePos = glGetAttribLocation(prog, "lsgl_LineWidth");
  printf("lsizePos = %d\n", lsizePos);
  glBindBuffer(GL_ARRAY_BUFFER, lradius);
  glVertexAttribPointer(lsizePos, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
  glEnableVertexAttribArray(lsizePos);
  assert(glGetError() == GL_NO_ERROR);

  GLfloat eye[3] = {300.0, 300.0, 300.0};
  GLfloat lookat[3] = {0.0, 0.0, 0.0};
  GLfloat up[3] = {0.0, 1.0, 0.0};

  lsglSetCamera(eye, lookat, up, 45.0f);

  GLfloat resolution[2];
  resolution[0] = (float)windowWidth;
  resolution[1] = (float)windowHeight;
  glUniform2fv(glGetUniformLocation(prog, "resolution"), 1, resolution);

  glBindBuffer(GL_ARRAY_BUFFER, lnvtx);
  glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  glEnableVertexAttribArray(attrPos);

  glDrawArrays(GL_LINES, 0, 2*numLines);
  assert(glGetError() == GL_NO_ERROR);

  glFinish();

  char buf[1024];
#ifdef ENABLE_MPI
  sprintf(buf, "colorbuf_%04d.tga", rank);
#else
  sprintf(buf, "colorbuf.tga");
#endif
  ret = SaveColorBufferRGBA(buf);
  assert(ret);

  glDeleteRenderbuffers(1, &colorRenderbuffer);
  glDeleteRenderbuffers(1, &depthRenderbuffer);
  glDeleteFramebuffers(1, &framebuffer);

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif

  return 0;
}
