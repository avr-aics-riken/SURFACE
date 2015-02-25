#ifdef _WIN32
#include <GLES2/gl2.h>
extern "C" {
__declspec(dllimport) void lsglSetCamera(float eye[3], float lookat[3], float up[3], float fov);
//__declspec(dllimport) void lsglSetPointSize(float partsize);
//__declspec(dllimport) void lsglSetPointSizev(int num, const float* partsize);
};
#else
#include "../../gles/gles_c_api.h"  // GLES + LSGL EXT.
//extern void lsglSetCamera(float eye[3], float lookat[3], float up[3], float fov);
//extern void lsglSetPointSize(float partsize);
//extern void lsglSetPointSizev(int num, const float* partsize);
#endif

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>

#include "../common/SimpleTGA.h"
//#include "particle_loader.h"
#include "timerutil.h"
#include "tinymt64.h"

#define USE_BINARY_SHADER (0)

//std::string particlename = "./nomat2.part";
GLfloat scenescale = 1.0;

int windowWidth = 512;
int windowHeight = 512;

template<typename T>
void GenerateRandomTetras(
  T* vertices, size_t n, double bmin[3], double bmax[3], int seed)
{

  assert(vertices);

  tinymt64_t rng;
  rng.mat1 = 0xfa051f40;
  rng.mat2 = 0xffd0fff4;
  rng.tmat = 0x58d02ffeffbfffbcULL;
  tinymt64_init(&rng, seed);

  float f = tinymt64_generate_double(&rng);

  //#ifdef _OPENMP
  //#pragma omp parallel for
  //#endif
  
  double offset[4][3];
  offset[0][0] = 1.0;
  offset[0][1] = 0.0;
  offset[0][2] = 0.0;
  offset[1][0] = 1.0;
  offset[1][1] = 0.0;
  offset[1][2] = 1.0;
  offset[2][0] = 0.0;
  offset[2][1] = 0.0;
  offset[2][2] = 1.0;
  offset[3][0] = 0.0;
  offset[3][1] = 1.0;
  offset[3][2] = 0.0;
  for (size_t i = 0; i < n; i++) {

    // [0, 1.0)
    double u0 = tinymt64_generate_double(&rng);
    double u1 = tinymt64_generate_double(&rng);
    double u2 = tinymt64_generate_double(&rng);

    // [bmin, bmax)
    double px = (bmax[0] - bmin[0]) * u0 + bmin[0];
    double py = (bmax[1] - bmin[1]) * u1 + bmin[1];
    double pz = (bmax[2] - bmin[2]) * u2 + bmin[2];

    for (size_t j = 0; j < 4; j++) {

#if 1
      double scale = 0.05;
      // random offset. [-0.025, 0.025)
      double v0 = scale * (bmax[0] - bmin[0]) * tinymt64_generate_double(&rng);
      double v1 = scale * (bmax[1] - bmin[1]) * tinymt64_generate_double(&rng);
      double v2 = scale * (bmax[2] - bmin[2]) * tinymt64_generate_double(&rng);

      v0 += 0.002 * (bmax[0] - bmin[0]) * offset[j][0];
      v1 += 0.002 * (bmax[1] - bmin[1]) * offset[j][1];
      v2 += 0.002 * (bmax[2] - bmin[2]) * offset[j][2];
#else 
      // fixed offset
      double v0 = 0.002 * (bmax[0] - bmin[0]) * offset[j][0];
      double v1 = 0.002 * (bmax[1] - bmin[1]) * offset[j][1];
      double v2 = 0.002 * (bmax[2] - bmin[2]) * offset[j][2];
#endif

      vertices[3*(4*i+j)+0] = px + v0;
      vertices[3*(4*i+j)+1] = py + v1;
      vertices[3*(4*i+j)+2] = pz + v2;

    }
  }
  
}

bool SaveColorBufferRGBA(const char* savename)
{
  void* tgabuffer;
  unsigned char* imgBuf = new unsigned char[windowWidth * windowHeight * 4];
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, imgBuf);

  int tgasize = SimpleTGASaverRGBA(&tgabuffer, windowWidth, windowHeight, imgBuf);
  delete [] imgBuf;
  if (!tgasize){
    printf("Failed save.\n");
    return false;
  }
 
  FILE* fp = fopen(savename, "wb");
  fwrite(tgabuffer, 1, tgasize, fp);
  fclose(fp);
  free(tgabuffer);
  return true;
}

//static bool
//LoadBinaryShader(
//  GLuint& prog,
//  GLuint& fragShader,
//  const char* fragShaderBinaryFilename)
//{
//  GLint val = 0;
//
//  // free old shader/program
//  if (prog != 0)   glDeleteProgram(prog);
//  if (fragShader != 0) glDeleteShader(fragShader);
//
//  std::vector<unsigned char> data;
//  FILE *fp = fopen(fragShaderBinaryFilename, "rb");
//  if (!fp) {
//    fprintf(stderr, "Failed to open file: %s\n", fragShaderBinaryFilename);
//    return false;
//  }
//  fseek(fp, 0, SEEK_END);
//  size_t len = ftell(fp);
//  rewind(fp);
//
//  data.resize(len);
//  len = fread(&data.at(0), 1, len, fp);
//  fclose(fp);
//
//  const void* binary = &data.at(0);
//  GLenum format = 0; // format is arbitrary at this time.
//
//  fragShader = glCreateShader(GL_FRAGMENT_SHADER);
//  glShaderBinary(1, &fragShader, format, binary, len);
//
//  prog = glCreateProgram();
//  glAttachShader(prog, fragShader);
//  glLinkProgram(prog);
//
//  glGetShaderiv(fragShader, GL_COMPILE_STATUS, &val);
//  assert(val == GL_TRUE && "failed to compile shader");
//
//  return true;
//}

static bool
LoadShader(
  GLuint& prog,
  GLuint& fragShader,
  const char* fragShaderSourceFilename)
{
  GLint val = 0;

  // free old shader/program
  if (prog != 0)   glDeleteProgram(prog);
  if (fragShader != 0) glDeleteShader(fragShader);

  static GLchar srcbuf[16384];
  FILE *fp = fopen(fragShaderSourceFilename, "rb");
  assert(fp);
  fseek(fp, 0, SEEK_END);
  size_t len = ftell(fp);
  rewind(fp);
  len = fread(srcbuf, 1, len, fp);
  srcbuf[len] = 0;
  fclose(fp);
    
  static const GLchar *src = srcbuf;
    
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

int
main(
  int argc,
  char **argv)
{
  int numTetras = 1000;

#ifdef ENABLE_MPI
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("[MPI] rank = %d\n", rank);
#endif

  if (argc > 1) {
    numTetras = atoi(argv[1]); 
  }
  printf("# of tetras = %d\n", numTetras);

  if (argc > 2) {
    scenescale = atof(argv[2]); 
  }
  printf("scenescale = %f\n", scenescale);

  glViewport(0, 0, windowWidth, windowHeight);

  GLuint prog = 0, fragShader = 0;

#if USE_BINARY_SHADER
  const char* fragShaderFile = "shader.so";
  bool ret = LoadBinaryShader(prog, fragShader, fragShaderFile);
#else
  const char* fragShaderFile = "input.frag";
  bool ret = LoadShader(prog, fragShader, fragShaderFile);
#endif

  assert(ret);

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
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, colorRenderbuffer);

  // create depth renderbuffer and attach
  GLuint depthRenderbuffer;
  glGenRenderbuffers(1, &depthRenderbuffer);
  glBindRenderbuffer(GL_RENDERBUFFER, depthRenderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, windowWidth, windowHeight);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthRenderbuffer);


  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);


  std::vector<float> positions;
  int numPoints = 0;

  numPoints = numTetras * 4; // 4 verts
  double bmax[3];
  double bmin[3];
  bmin[0] = bmin[1] = bmin[2] = -scenescale;
  bmax[0] = bmax[1] = bmax[2] =  scenescale;
  positions.resize(numPoints*3);
  GenerateRandomTetras<float>(&positions.at(0), numTetras, bmin, bmax, 1234);

  // gen indices
  std::vector<unsigned int> indices;
  for (int i = 0; i < 4*numTetras; i++) {
    indices.push_back(i);
  }

  // 1. Create Vertex Buffers.
  GLuint ptvtx;

  glGenBuffers(1, &ptvtx);
  glBindBuffer(GL_ARRAY_BUFFER, ptvtx);
  //glBufferData(GL_ARRAY_BUFFER, npoints * sizeof(float) * 3, part->positions, GL_STATIC_DRAW);
  lsglBufferDataPointer(GL_ARRAY_BUFFER, numPoints * sizeof(float) * 3, &positions.at(0), GL_STATIC_DRAW);

  // Camera
  //float eye[3] = {10.0f, 15.0f, 20.0f};   // bunny x 100
  //float eye[3] = {10.0f, 200.0f, 100.0f};   // bunny x 100
  float lookat[3] = {0.0f, 0.0f, 0.0f}; 
  //float eye[3] = {100.0f, 200.0f, 400.0f};  // PDB
  //float up[3] = {0.0f, 1.0f, 0.0f}; 

  // PDB Z up
  //float eye[3] = {50.0f, 50.0f, 150.0f};  // PDB
  float eye[3] = {2.0f, 2.0f, 1.0f}; // tetra
  float up[3] = {0.0f, 1.0f, 0.0f}; 
  lsglSetCamera(eye, lookat, up, 45.0f);

  // enable subsampling
  //glEnable(GL_SAMPLE_COVERAGE);
  //glSampleCoverage(2.0f, GL_FALSE);

  GLfloat resolution[2];
  resolution[0] = (float)windowWidth;
  resolution[1] = (float)windowHeight;
  glUniform2fv(glGetUniformLocation(prog, "resolution"), 1, resolution);

  // 2. Use vertex buffers.

  glBindBuffer(GL_ARRAY_BUFFER, ptvtx);
  glVertexAttribPointer(attrPos, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (void*)0);
  glEnableVertexAttribArray(attrPos);
  assert(glGetError() == GL_NO_ERROR);

  glDrawArrays(GL_TETRAHEDRONS_EXT, 0, 4*numTetras);
  assert(glGetError() == GL_NO_ERROR);

  timerutil t;
  t.start();
  glFinish();
  assert(glGetError() == GL_NO_ERROR);
  t.end();
  printf("Render time: %d msec\n", (int)t.msec());

  ret = SaveColorBufferRGBA("colorbuf.tga");
  assert(glGetError() == GL_NO_ERROR);
  assert(ret);

  glDeleteRenderbuffers(1, &colorRenderbuffer);
  glDeleteRenderbuffers(1, &depthRenderbuffer);
  glDeleteFramebuffers(1, &framebuffer);

  return 0;
}
