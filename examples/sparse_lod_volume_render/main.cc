#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#include "../../gles/gles_c_api.h"
#include <GLES2/gl2ext.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <string>

#include "../common/SimpleTGA.h"

#define USE_DOUBLE_PRECISION  (1)

int windowWidth = 512;
int windowHeight = 512;

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

static bool GenLodSparseVolumeFloatTexture(GLuint &tex,
                                        float **bufAddr,
                                        const unsigned char *voldata,
                                        int dim[3]) {

  float *buf = new float[dim[0] * dim[1] * dim[2]];

  // uchar -> float
  for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
    buf[i] = voldata[i] / 255.0f;
  }

  int nlods = 4; // Up to 8 for demo;
  int loc[8][3];

  // Compute block origin.
  loc[0][0] = 0;
  loc[0][1] = 0;
  loc[0][2] = 0;
  for (int i = 1; i < nlods; i++) {

    loc[i][0] = 0; //loc[i-1][0] + dim[0] * pow(2.0, i-1);
    loc[i][1] = 0; //loc[i-1][1] + dim[1] * pow(2.0, i-1);
    loc[i][2] = loc[i-1][2] + dim[2] * pow(2.0, i-1);
  }
  

  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_3D, tex);

  for (int i = 0; i < nlods; i++) {
  
    int extent[3] = {0,0,0};

    extent[0] = dim[0] * pow(2.0, i);
    extent[1] = dim[1] * pow(2.0, i);
    extent[2] = dim[2] * pow(2.0, i);

    printf("loc    = %d, %d, %d\n", loc[i][0], loc[i][1], loc[i][2]);
    printf("extent = %d, %d, %d\n", extent[0], extent[1], extent[2]);
    printf("dim    = %d, %d, %d\n", dim[0], dim[1], dim[2]);

    lsglTexPageCommitment(GL_TEXTURE_3D, /* lv= */ 0, loc[i][0],
                          loc[i][1], loc[i][2], extent[0], extent[1], extent[2],
                          GL_TRUE);
    assert(glGetError() == GL_NO_ERROR);
    lsglTexSubImage3DPointer(GL_TEXTURE_3D, 0, loc[i][0], loc[i][1], loc[i][2],
                      extent[0], extent[1], extent[2], dim[0], dim[1], dim[2], GL_LUMINANCE, GL_FLOAT, buf);
    assert(glGetError() == GL_NO_ERROR);

  }

  // Should not delete buf here. Return the pointer, and delete it after the rendering
  (*bufAddr) = buf;

  return true;
}

static bool GenLodSparseVolumeDoubleTexture(GLuint &tex,
                                        double **bufAddr,
                                        const unsigned char *voldata,
                                        int dim[3]) {

  double *buf = new double[dim[0] * dim[1] * dim[2]];

  // uchar -> float
  for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
    buf[i] = voldata[i] / 255.0;
  }

  int nlods = 4; // Up to 8 for demo;
  int loc[8][3];

  // Compute block origin.
  loc[0][0] = 0;
  loc[0][1] = 0;
  loc[0][2] = 0;
  for (int i = 1; i < nlods; i++) {

    loc[i][0] = 0; //loc[i-1][0] + dim[0] * pow(2.0, i-1);
    loc[i][1] = 0; //loc[i-1][1] + dim[1] * pow(2.0, i-1);
    loc[i][2] = loc[i-1][2] + dim[2] * pow(2.0, i-1);
  }
  

  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_3D, tex);

  for (int i = 0; i < nlods; i++) {
  
    int extent[3] = {0,0,0};

    extent[0] = dim[0] * pow(2.0, i);
    extent[1] = dim[1] * pow(2.0, i);
    extent[2] = dim[2] * pow(2.0, i);

    printf("loc    = %d, %d, %d\n", loc[i][0], loc[i][1], loc[i][2]);
    printf("extent = %d, %d, %d\n", extent[0], extent[1], extent[2]);
    printf("dim    = %d, %d, %d\n", dim[0], dim[1], dim[2]);

    lsglTexPageCommitment(GL_TEXTURE_3D, /* lv= */ 0, loc[i][0],
                          loc[i][1], loc[i][2], extent[0], extent[1], extent[2],
                          GL_TRUE);
    assert(glGetError() == GL_NO_ERROR);
    lsglTexSubImage3DPointer(GL_TEXTURE_3D, 0, loc[i][0], loc[i][1], loc[i][2],
                      extent[0], extent[1], extent[2], dim[0], dim[1], dim[2], GL_LUMINANCE, GL_DOUBLE, buf);
    assert(glGetError() == GL_NO_ERROR);

  }

  // Should not delete buf here. Return the pointer, and delete it after the rendering
  (*bufAddr) = buf;

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

static unsigned char *LoadRawVolumeTexture(const char *filename, size_t size) {
  FILE *fp = fopen(filename, "rb");
  assert(fp);

  unsigned char *data = new unsigned char[size];
  size_t sz = fread(data, 1, size, fp);
  assert(sz == size);
  assert(data);

  fclose(fp);

  return data;
}

int main(int argc, char **argv) {
#ifdef LSGL_ENABLE_MPI
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  const char *fragShaderFile = "input.frag";
  const char *fragShaderBinFile = "shader.so";
  const char *textureFile = "Bucky.raw";
  int dim[3] = {32, 32, 32};

  size_t sz = dim[0] * dim[1] * dim[2] * sizeof(unsigned char);
  unsigned char *voldata = LoadRawVolumeTexture(textureFile, sz);

  GLuint prog = 0, fragShader = 0;
  bool ret = LoadShader(prog, fragShader, fragShaderFile);
  // bool ret = LoadBinaryShader(prog, fragShader, fragShaderBinFile);
  assert(ret);

  printf("LoadTex\n");
  GLuint tex0;
#if USE_DOUBLE_PRECISION
  double *bufAddr = NULL;
  ret = GenLodSparseVolumeDoubleTexture(tex0, &bufAddr, voldata, dim);
#else
  float *bufAddr = NULL;
  ret = GenLodSparseVolumeFloatTexture(tex0, &bufAddr, voldata, dim);
#endif
  assert(ret);

  glActiveTexture(GL_TEXTURE0);

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

  glViewport(0, 0, windowWidth, windowHeight);

  GLint numTexUnits;
  glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &numTexUnits);
  printf("Max texture units = %d\n", numTexUnits);

  GLboolean ok = glIsProgram(prog);
  assert(ok == GL_TRUE);

  GLint matLoc = glGetUniformLocation(prog, "tmat");
  printf("matLoc: %d\n", matLoc);
  float m[16];
  m[0] = 1.0;
  m[1] = 0.5;
  m[2] = 1.0;
  m[3] = 4.0;
  m[4] = 0.0;
  m[5] = 6.0;

  GLint texture0Loc = glGetUniformLocation(prog, "tex0");
  printf("texture0Loc: %d\n", texture0Loc);

  glUseProgram(prog);

  glUniformMatrix4fv(matLoc, 1, GL_FALSE, m);
  glUniform1i(texture0Loc, 0); // Texture unit 0 is for base images.

  GLfloat eye[3] = {5.0, 0.0, 0.0};
  GLfloat lookat[3] = {0.0, 0.0, 0.0};
  GLfloat up[3] = {0.0, 1.0, 0.0};

  GLfloat resolution[2];
  resolution[0] = (float)windowWidth;
  resolution[1] = (float)windowHeight;
  glUniform2fv(glGetUniformLocation(prog, "resolution"), 1, resolution);

  glUniform3fv(glGetUniformLocation(prog, "eye"), 1, eye);
  glUniform3fv(glGetUniformLocation(prog, "lookat"), 1, lookat);
  glUniform3fv(glGetUniformLocation(prog, "up"), 1, up);

  printf("eval\n");
  lsglEvalFragmentShader();
  printf("eval DONE\n");

  char buf[1024];
#ifdef LSGL_ENABLE_MPI
  sprintf(buf, "colorbuf_%04d.tga", rank);
#else
  sprintf(buf, "colorbuf.tga");
#endif
  ret = SaveColorBufferRGBA(buf);
  assert(ret);

  glDeleteRenderbuffers(1, &colorRenderbuffer);
  glDeleteRenderbuffers(1, &depthRenderbuffer);
  glDeleteFramebuffers(1, &framebuffer);

  delete [] bufAddr;

#ifdef LSGL_ENABLE_MPI
  MPI_Finalize();
#endif

  return 0;
}
