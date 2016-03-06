#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#include "../../gles/gles_c_api.h"  // GLES + LSGL EXT.

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>

#include "../common/SimpleTGA.h"
#include "timerutil.h"
#include "tinymt64.h"

#define USE_BINARY_SHADER (0)

GLfloat scenescale = 1.0;

int windowWidth = 512;
int windowHeight = 512;

struct myvec4 {
  float x;
  float y;
  float z;
  float w;
};

struct float3 {
  float3() {}
  float3(float xx, float yy, float zz) {
    x = xx;
    y = yy;
    z = zz;
  }
  float3(const float *p) {
    x = p[0];
    y = p[1];
    z = p[2];
  }

  float3 operator*(float f) const { return float3(x * f, y * f, z * f); }
  float3 operator-(const float3 &f2) const {
    return float3(x - f2.x, y - f2.y, z - f2.z);
  }
  float3 operator*(const float3 &f2) const {
    return float3(x * f2.x, y * f2.y, z * f2.z);
  }
  float3 operator+(const float3 &f2) const {
    return float3(x + f2.x, y + f2.y, z + f2.z);
  }
  float3 &operator+=(const float3 &f2) {
    x += f2.x;
    y += f2.y;
    z += f2.z;
    return (*this);
  }
  float3 operator/(const float3 &f2) const {
    return float3(x / f2.x, y / f2.y, z / f2.z);
  }
  float operator[](int i) const { return (&x)[i]; }
  float &operator[](int i) { return (&x)[i]; }

  float3 neg() { return float3(-x, -y, -z); }

  float length() { return sqrtf(x * x + y * y + z * z); }

  void normalize() {
    float len = length();
    if (fabs(len) > 1.0e-6f) {
      float inv_len = 1.0 / len;
      x *= inv_len;
      y *= inv_len;
      z *= inv_len;
    }
  }

  float x, y, z;
  // float pad;  // for alignment
};

inline float3 operator*(float f, const float3 &v) {
  return float3(v.x * f, v.y * f, v.z * f);
}

inline float3 vcross(float3 a, float3 b) {
  float3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

inline float vdot(float3 a, float3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline float rand_(float size){
    return (float)rand() / RAND_MAX * size;
}

inline float rand_(){
    return (float)rand() / RAND_MAX;
}

inline float3 rand_vertices(){
    return float3(rand_(), rand_(), rand_());
}

static void SpaceDivision(float *bmin, float *bmax,int number, std::vector< std::vector<float > > &out_d_min, std::vector< std::vector<float > > &out_d_max){
    
    int i = 1,j = 1, k = 1;
    while (1){
        i++;
        if(i*j*k >= number) break;
        j++;
        if(i*j*k >= number) break;
        k++;
        if(i*j*k >= number) break;
    }
    
    int n = 0, ir, jr, kr;
    bool *room = new bool[i*j*k];
    
    for(int m = 0; m < i*j*k; m ++) room[m] = false;
    
    out_d_min.resize(number);
    out_d_max.resize(number);
    
    while(n != number){
        ir = rand_(i);
        jr = rand_(j);
        kr = rand_(k);
        
        if(room[(ir * j + jr) * k + kr]) continue;
        
        room[(ir * j + jr) * k + kr] = true;
        
        out_d_min[n].resize(3);
        out_d_min[n][0] = bmin[0] + (bmax[0] - bmin[0]) / i * ir;
        out_d_min[n][1] = bmin[1] + (bmax[1] - bmin[1]) / j * jr;
        out_d_min[n][2] = bmin[2] + (bmax[2] - bmin[2]) / k * kr;
        
        out_d_max[n].resize(3);
        out_d_max[n][0] = bmin[0] + (bmax[0] - bmin[0]) / i * (ir + 1);
        out_d_max[n][1] = bmin[1] + (bmax[1] - bmin[1]) / j * (jr + 1);
        out_d_max[n][2] = bmin[2] + (bmax[2] - bmin[2]) / k * (kr + 1);
        
        n++;
        
    }
    
    delete [] room;
    
//    cout << " i , j , k = " << i << " , " << j << " , " << k << endl;
//    
//    for(int m = 0; m < out_d_max.size(); m ++){
//        cout << "out_d_min[" << m << "] = " << out_d_min[m][0] << " , " << out_d_min[m][1] << " , " << out_d_min[m][2] << endl;
//        cout << "out_d_max[" << m << "] = " << out_d_max[m][0] << " , " << out_d_max[m][1] << " , " << out_d_max[m][2] << endl;
//    }
}

static bool RandomTetraMake(float *vertices, float *bmin, float *bmax){
    
    float scale = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.25f;
    scale = pow(scale, 0.333);
    
    float3 p[4];
    
    p[0] = scale*2*rand_vertices();
    p[1] = scale*2*rand_vertices();
    p[2] = scale*2*rand_vertices();
    p[3] = scale*2*rand_vertices();
    
    while( (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.3f > -vdot( p[3] - p[0], vcross(p[1] - p[0], p[2] - p[1]))){
        p[0] = scale*2*rand_vertices();
        p[1] = scale*2*rand_vertices();
        p[2] = scale*2*rand_vertices();
        p[3] = scale*2*rand_vertices();
//        cout << (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.3f << endl;
    }
    bool f = false;
    
    for(int n = 0; n < 4; n ++ ){
        p[n] = p[n] + float3(bmin);
        if( p[n].x > bmax[0]) f = true;
        if( p[n].y > bmax[1]) f = true;
        if( p[n].z > bmax[2]) f = true;
    }
    
    if(f) return false;
    
    for(int n = 0; n < 4; n ++ ){
        vertices[n*3 + 0] = p[n].x;
        vertices[n*3 + 1] = p[n].y;
        vertices[n*3 + 2] = p[n].z;
    }
    
    return true;
}

static void RandomTetraMake(float *vertices, float *bmin, float *bmax, int PrismN){
    
    std::vector<std::vector<float> > bmin_, bmax_;
    
    SpaceDivision(bmin, bmax, PrismN, bmin_, bmax_);
    
    for(int i = 0; i < PrismN;){
        float vs[12];
        if(RandomTetraMake(vs, &bmin_[i][0], &bmax_[i][0])){
            for(int n = 0; n < 12; n ++)
                vertices[i * 12 + n] = vs[n];
            i++;
        }
    }
}

static bool RandomPyramidMake(float *vertices, float *bmin, float *bmax){
    
    float scale = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.25f;
    scale = pow(scale, 0.333);
    
    float3 p[5];
    
    while(1){
        p[0] = scale*2*rand_vertices();
        p[1] = scale*2*rand_vertices();
        p[2] = scale*2*rand_vertices();
        p[3] = p[2] + (p[2] - p[1]) * rand_(scale * 2) + (p[0] - p[1]) * rand_(scale * 2);
        p[4] = scale*2*rand_vertices();
        if ((bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.3f >
            -vdot( p[3] - p[0], vcross(p[1] - p[0], p[2] - p[1])))
            break;
    }
    bool f = false;
    
    for(int n = 0; n < 5; n ++ ){
        p[n] = p[n] + float3(bmin);
        if( p[n].x > bmax[0]) f = true;
        if( p[n].y > bmax[1]) f = true;
        if( p[n].z > bmax[2]) f = true;
    }
    
    if(f) return false;
    
    for(int n = 0; n < 5; n ++ ){
        vertices[n*3 + 0] = p[n].x;
        vertices[n*3 + 1] = p[n].y;
        vertices[n*3 + 2] = p[n].z;
    }
    
    return true;
}

static void RandomPyramidMake(std::vector<float> &vertices, float *bmin, float *bmax, int PyramidN){
    std::vector<std::vector<float> > bmin_, bmax_;
    
    SpaceDivision(bmin, bmax, PyramidN, bmin_, bmax_);
    
    for(int i = 0; i < PyramidN;){
        float vs[15];
        if(RandomTetraMake(vs, &bmin_[i][0], &bmax_[i][0])){
            for(int n = 0; n < 15; n ++)
                vertices[i * 15 + n] = vs[n];
            i++;
        }
    }
    
}

static bool RandomPrismMake(float *vertices, float *bmin, float *bmax){
//    float3 center = float3(bmin[0]+bmax[0], bmin[1]+bmax[1], bmin[2]+bmax[2]) * 0.5f;
    
//    glm::vec4 center4(center.x,center.y,center.z,0);
    
    float scale = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.25f;
    scale = pow(scale, 0.333);
    
    float3 p[6];
    
    p[0] = float3(rand_(scale),rand_(scale),rand_(scale));
    p[1] = float3(rand_(scale),rand_(scale),rand_(scale));
    p[2] = float3(rand_(scale),rand_(scale),rand_(scale));
    float3 p_ = float3(rand_(scale),rand_(scale),rand_(scale));
    
    while( (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.1f > -vdot( p_ - p[0], vcross(p[1] - p[0], p[2] - p[1]))){
        p[0] = float3(rand_(scale),rand_(scale),rand_(scale));
        p[1] = float3(rand_(scale),rand_(scale),rand_(scale));
        p[2] = float3(rand_(scale),rand_(scale),rand_(scale));
        p_ = float3(rand_(scale),rand_(scale),rand_(scale));
    }
    
    p[3] = p[0] + (p_ - p[0]) * rand_();
    p[4] = p[2] + (p_ - p[2]) * rand_();
    p[5] = p[1] + (p_ - p[1]) * rand_();
    
    for(int n = 0; n < 6; n ++ ){
        p[n] = p[n] + float3(bmin);
        if( p[n].x > bmax[0] ) return false;
        if( p[n].y > bmax[1] ) return false;
        if( p[n].z > bmax[2] ) return false;
    }
    
    for(int n = 0; n < 6; n ++ ){
        vertices[n*3 + 0] = p[n].x;
        vertices[n*3 + 1] = p[n].y;
        vertices[n*3 + 2] = p[n].z;
    }
    
    return true;
    
    
}

static void RandomPrismMake(std::vector<float> &vertices, float *bmin, float *bmax, int PrismN){
    
    std::vector<std::vector<float> > bmin_, bmax_;
    
    SpaceDivision(bmin, bmax, PrismN, bmin_, bmax_);
    
    for(int i = 0; i < PrismN;){
        float vs[18];
        if(RandomPrismMake(vs, &bmin_[i][0], &bmax_[i][0])){
            for(int n = 0; n < 18; n ++)
                vertices[i * 18 + n] = vs[n];
            i++;
        }
    }
    
        
}

bool RandomHexaMake(float *vertices, float *bmin, float *bmax){
    
    float3 min(bmin);
    
    float scale = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1])*0.1f;
    scale = pow(scale, 0.333);
    
    float3 verts[8];
    
    while (1){
        verts[0].x = 0;
        verts[0].y = 0;
        verts[0].z = -0.5;
        
        verts[1].x = 0;
        verts[1].y = 1;
        verts[1].z = -0.5;
        
        verts[2].x = 1;
        verts[2].y = 1;
        verts[2].z = -0.5;
        
        verts[3].x = 1;
        verts[3].y = 0;
        verts[3].z = -0.5;
        
        verts[4].x = 0;
        verts[4].y = 0;
        verts[4].z = 0.5;
        
        verts[5].x = 1;
        verts[5].y = 0;
        verts[5].z = 0.5;
        
        verts[6].x = 1;
        verts[6].y = 1;
        verts[6].z = 0.5;
        
        verts[7].x = 0;
        verts[7].y = 1;
        verts[7].z = 0.5;
        
        float mat[4][3] = {
            {rand_(scale),rand_(scale),rand_(scale)},
            {rand_(scale),rand_(scale),rand_(scale)},
            {rand_(scale),rand_(scale),rand_(scale)},
            {rand_(scale),rand_(scale),rand_(scale)},
        };
        
        
        
        bool f = false;
        
        for(int n = 0; n < 8; n++){
            
            verts[n] = float3( mat[0][0] * verts[n].x + mat[1][0] * verts[n].y + mat[2][0] * verts[n].z + mat[3][0],
                              mat[0][1] * verts[n].x + mat[1][1] * verts[n].y + mat[2][1] * verts[n].z + mat[3][1],
                              mat[0][2] * verts[n].x + mat[1][2] * verts[n].y + mat[2][2] * verts[n].z + mat[3][2])
            + min;
            
            if( verts[n].x > bmax[0] || verts[n].x < bmin[0]) f = true;
            if( verts[n].y > bmax[1] || verts[n].y < bmin[1]) f = true;
            if( verts[n].z > bmax[2] || verts[n].z < bmin[2]) f = true;
        }
        
        if(f) return false;
        
        if((bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) * 0.06f
           < - vdot( float3(verts[4].x - verts[0].x, verts[4].y - verts[0].y, verts[4].z - verts[0].z), vcross(float3(verts[1].x - verts[0].x, verts[1].y - verts[0].y, verts[1].z - verts[0].z ), float3(verts[2].x - verts[1].x, verts[2].y - verts[1].y, verts[2].z - verts[2].z))))break;
    }
    
    for(int n = 0; n < 8; n++){
        vertices[n*3] = verts[n].x;
        vertices[n*3+1] = verts[n].y;
        vertices[n*3+2] = verts[n].z;
    }
    
    return true;
}

void RandomHexaMake(std::vector<float> &vertices, float *bmin, float *bmax, int HexaN){
    
    std::vector<std::vector<float> > bmin_, bmax_;
    
    SpaceDivision(bmin, bmax, HexaN, bmin_, bmax_);
    
    for(int i = 0; i < HexaN;){
        float vs[24];
        if(RandomHexaMake(vs, &bmin_[i][0], &bmax_[i][0])){
            for(int n = 0; n < 24; n ++)
                vertices[i * 24 + n] = vs[n];
            i++;
        }
    }
}

#if 0
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
#endif

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
  if (!fp) {
    printf("File not found: %s\n", fragShaderSourceFilename);
    exit(-1);
  }
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
  int numSolids = 1000;
  int solidType = 5;

#ifdef LSGL_ENABLE_MPI
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("[MPI] rank = %d\n", rank);
#endif

  printf("Usage: %s <solidType> <numSolids> <sceneScale>\n", argv[0]);
  printf("    solidType = 5(Pyramid), 6(Prism), 8(Hexa)\n");

  if (argc > 1) {
    solidType = atoi(argv[1]); 
  }
  if ((solidType == 5) || (solidType == 6) || (solidType == 8)) {
    // OK
  } else {
    printf("Invalid solidType: %d\n", solidType);
    exit(-1);
  }

  if (argc > 2) {
    numSolids = atoi(argv[2]); 
  }
  printf("# of solids = %d\n", numSolids);

  if (argc > 3) {
    scenescale = atof(argv[3]); 
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

  numPoints = numSolids * solidType; // solidType == numVert
  float bmax[3];
  float bmin[3];
  bmin[0] = bmin[1] = bmin[2] = -scenescale;
  bmax[0] = bmax[1] = bmax[2] =  scenescale;

  positions.resize(numSolids * 3 * solidType); // solidType == numVert
  if (solidType == 5) {
    printf("GenPyramid...\n");
    RandomPyramidMake(positions, bmin, bmax, numSolids);
  } else if (solidType == 6) {
    printf("GenPrism...\n");
    RandomPrismMake(positions, bmin, bmax, numSolids);
  } else if (solidType == 8) {
    printf("GenHexa...\n");
    RandomHexaMake(positions, bmin, bmax, numSolids);
  }
  printf("Done generating test data...\n");

  // gen indices
  std::vector<unsigned int> indices;
  for (int i = 0; i < solidType*numSolids; i++) {
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

  if (solidType == 5) {
    glDrawArrays(GL_PYRAMIDS_EXT, 0, 5*numSolids);
  } else if (solidType == 6) {
    glDrawArrays(GL_PRISMS_EXT, 0, 6*numSolids);
  } else if (solidType == 8) {
    glDrawArrays(GL_HEXAHEDRONS_EXT, 0, 8*numSolids);
  }

  assert(glGetError() == GL_NO_ERROR);

  timerutil t;
  t.start();
  glFinish();
  assert(glGetError() == GL_NO_ERROR);
  t.end();
  printf("Render time: %d msec\n", (int)t.msec());

  char buf[1024];
#ifdef LSGL_ENABLE_MPI
  sprintf(buf, "colorbuf_%06d.tga", rank);
#else
  sprintf(buf, "colorbuf.tga");
#endif
  ret = SaveColorBufferRGBA("colorbuf.tga");
  assert(glGetError() == GL_NO_ERROR);
  assert(ret);

  glDeleteRenderbuffers(1, &colorRenderbuffer);
  glDeleteRenderbuffers(1, &depthRenderbuffer);
  glDeleteFramebuffers(1, &framebuffer);

#ifdef LSGL_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
