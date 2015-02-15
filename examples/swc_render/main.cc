#include "../../gles/gles_c_api.h"
#include <GLES2/gl2ext.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <string>

#include "../common/SimpleTGA.h"
#include "tiny_swc.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

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

  if (argc < 2) {
    fprintf(stderr, "need input.pdb\n");
    exit(-1);
  }

  tinyswc::TinySWC swc(argv[1]);
  if (!swc.Parse()) {
    fprintf(stderr, "SWC parsing failed: %s \n", argv[1]);
    return -1;
  }

  //const unsigned int npoints = pdb.GetAtoms().size();

  //std::vector<float> positions(3 * npoints);

  //for (unsigned int i = 0; i < npoints; i++) {
  //  positions[3*i+0] = pdb.GetAtoms()[i].GetX();
  //  positions[3*i+1] = pdb.GetAtoms()[i].GetY();
  //  positions[3*i+2] = pdb.GetAtoms()[i].GetZ();
  //}

  // Render Bond as line primitive.
  // @todo { remove duplicated bonds. }
  std::vector<float> lines;
  std::vector<tinyswc::SamplePoint>& segments = swc.GetNeuronSegments();

  for (unsigned int i = 0; i < segments.size() / 2; i++) {

    tinyswc::SamplePoint& s = segments[2*i+0];
    tinyswc::SamplePoint& e = segments[2*i+1];

    lines.push_back(s.GetX()); 
    lines.push_back(s.GetY()); 
    lines.push_back(s.GetZ()); 

    lines.push_back(e.GetX()); 
    lines.push_back(e.GetY()); 
    lines.push_back(e.GetZ()); 
      
  }

  const unsigned int nlines = lines.size() / (3 * 2);

  // gen indices
  //std::vector<unsigned int> pointIndices;
  //for (int i = 0; i < npoints; i++) {
  //  pointIndices.push_back(i);
  //}

  GLuint prog = 0, fragShader = 0;
  bool ret = LoadShader(prog, fragShader, fragShaderFile);
  if (!ret) {
    fprintf(stderr, "failed to load shader: %s\n", fragShaderFile);
    exit(-1);
  }

  glUseProgram(prog);

  // Used AA samples per pixel in LSGL.
  glSampleCoverage(2.0f, GL_FALSE);

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

  // 1. Create Vertex Buffers.
  GLuint ptvtx, ptmat, ptidx, ptxradius;
  GLuint lnvtx, lnidx;

  //glGenBuffers(1, &ptvtx);
  //glBindBuffer(GL_ARRAY_BUFFER, ptvtx);
  //lsglBufferDataPointer(GL_ARRAY_BUFFER, npoints * sizeof(float) * 3, &positions.at(0), GL_STATIC_DRAW);

  glGenBuffers(1, &lnvtx);
  glBindBuffer(GL_ARRAY_BUFFER, lnvtx);
  lsglBufferDataPointer(GL_ARRAY_BUFFER, nlines * sizeof(float) * 3 * 2, &lines.at(0), GL_STATIC_DRAW);

  GLfloat eye[3] = {95.0, 95.0, 95.0};
  GLfloat lookat[3] = {0.0, 0.0, 0.0};
  GLfloat up[3] = {0.0, 1.0, 0.0};

  lsglSetCamera(eye, lookat, up, 45.0f);

  // Set point size using LSGL extension.
  glUniform1f(glGetUniformLocation(prog, "lsgl_PointSize"), 0.8);

  // Set line size using GLES interface.
  glLineWidth(1.0f);

  GLfloat resolution[2];
  resolution[0] = (float)windowWidth;
  resolution[1] = (float)windowHeight;
  glUniform2fv(glGetUniformLocation(prog, "resolution"), 1, resolution);

  // 2.1 draw points
  //glBindBuffer(GL_ARRAY_BUFFER, ptvtx);
  //glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  //glEnableVertexAttribArray(attrPos);
  //assert(glGetError() == GL_NO_ERROR);

  //glDrawArrays(GL_POINTS, 0, npoints);
  //assert(glGetError() == GL_NO_ERROR);

  // 2.2 draw lines
  glBindBuffer(GL_ARRAY_BUFFER, lnvtx);
  glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  glEnableVertexAttribArray(attrPos);
  assert(glGetError() == GL_NO_ERROR);

  glDrawArrays(GL_LINES, 0, 2*nlines);
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
