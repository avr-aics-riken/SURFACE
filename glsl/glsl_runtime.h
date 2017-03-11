#ifndef __GLSL_CC_SHADER_RUNTIME_H__
#define __GLSL_CC_SHADER_RUNTIME_H__

#include <stdio.h>
#include <stdlib.h>

#define MAX_TEXTURE_UNITS               (32) // FIXME(syoyo): Currently this must be same with MAX_FRAGMENT_UNIFORM_VECTORS
#define MAX_FRAGMENT_UNIFORM_VECTORS    (32)
#define MAX_FRAGMENT_VARYING_VARIABLES  (8)
#define MAX_DRAW_BUFFERS                (4)

#define GLSL_QUALIFIER_IN               (1 << 1)
#define GLSL_QUALIFIER_OUT              (1 << 2)
#define GLSL_QUALIFIER_UNIFORM          (1 << 3)

#if defined(_WIN32) && defined(LSGL_EXPORT)
#define GLSL_RT_DLLEXPORT __declspec(dllexport)
#else
#define GLSL_RT_DLLEXPORT
#endif

#define GLSL_ENABLE_ASSERTION (1)

#if GLSL_ENABLE_ASSERTION
#define GLSL_ASSERT(cond) { \
  if (!(cond)) { \
    fprintf(stderr, "[GLSL] " #cond " : %s(%s:%d)\n", __FUNCTION__, __FILE__, __LINE__); \
    abort(); \
  } \
}
#else
#define GLSL_ASSERT(cond) (void)(cond)
#endif

#define OPTIMIZED_TEXTURING (0)

#if OPTIMIZED_TEXTURING
#include "../render/texture.h"
#endif

extern "C" {

typedef struct
{
    const char* name;       ///< Name of type.
    int         n;          ///< Vector length of type.
} GLSLType;

typedef struct
{
    GLSLType    type;       ///< Type of variable.
    int         qualifier;  ///< Qualifier of variable.
    const char* name;       ///< Name of variable.
} GLSLUniformInfo;

typedef struct
{
    GLSLType    type;       ///< Type of variable.
    int         qualifier;  ///< Qualifier of variable.
    const char* name;       ///< Name of variable.
} GLSLVaryingInfo;

typedef struct
{
    GLSLUniformInfo uniformInfos[MAX_FRAGMENT_UNIFORM_VECTORS];
    int             numUniforms;

    GLSLVaryingInfo varyingInfos[MAX_FRAGMENT_VARYING_VARIABLES];
    int             numVaryings;
} FragmentConfig;


typedef struct
{
    float v[2];
} vec2;

typedef struct
{
    float v[3];
} vec3;

typedef struct
{
    float v[4];
} vec4;

typedef struct
{
    double v[2];
} dvec2;

typedef struct
{
    double v[3];
} dvec3;

typedef struct
{
    double v[4];
} dvec4;

typedef struct
{
    int v[2];
} ivec2;

typedef struct
{
    int v[3];
} ivec3;

typedef struct
{
    int v[4];
} ivec4;

typedef struct
{
    char v[2];
} bvec2;

typedef struct
{
    char v[3];
} bvec3;

typedef struct
{
    char v[4];
} bvec4;

typedef struct
{
    float v[2][2];
} mat2;

typedef struct
{
    float v[3][3];
} mat3;

typedef struct
{
    float v[4][4];
} mat4;

// Function pointer for builtin functions.
typedef vec4 (*texture2DFuncPtr)(unsigned long long samplerID, vec2* coords);
typedef vec4 (*texture3DFuncPtr)(unsigned long long samplerID, vec3* coords);
typedef float (*shadowFuncPtr)(void* frag, vec3* pos, vec3* dir);
typedef float (*traceFuncPtr)(void* frag, vec3* pos, vec3* dir, vec4* col, float* userattrib);
typedef float (*randomFuncPtr)(int threadID);

// Texture sampler
typedef unsigned long long sampler2D;
typedef unsigned long long sampler3D;

// Builtin fragment variables.
typedef struct
{
    // input
    float texcoord[4];
    float depth;                        // gl_Depth
    float fragCoord[4];                 // gl_FragCoord
    int   frontFacing;                  // gl_FrontFacing
    float pointCoord[2];                // gl_PointCoord

    //
    // output
    //
    float fragColor[4];                 // gl_FragColor
    float fragData[MAX_DRAW_BUFFERS];   // gl_FragData
    int   fragDiscarded;                // 1 if discarded

    // For trace() function
    float position[3];
    float normal[3];
    float geometricNormal[3];
    float tangent[3];
    float binormal[3];
    float indir[3];
    float barycentric[2];
    float hitdist;
    int   raydepth;
    float px, py;
    int   doubleSided;
    unsigned int prev_prim_id;
    const unsigned char* prev_node; 
    float prev_hit_t; 
    float prev_hit_normal[3];
  
    float rayattrib;                    // User ray attribute

    float cameraFrame[3][3];            // (eye, lookat, up)
    float cameraFov;                    // field of view

    int   threadID;                     // Thread ID

    texture2DFuncPtr    texture2D;
    texture3DFuncPtr    texture3D;
    shadowFuncPtr       shadow;
    traceFuncPtr        trace;
    randomFuncPtr       random;
} Fragment;

// Uniform variable
typedef struct
{
    unsigned char* data;                // opeque data
} FragmentUniform;

// Varying variable
typedef struct
{
    unsigned char* data;                // opeque data
} FragmentVarying;

typedef struct
{
    // FIXME(syoyo): Currently texture slot index is shared with uniforms, thus textures must have same array size of uniforms, not MAX_TEXTURE_UNITS.
    unsigned long long textures[MAX_FRAGMENT_UNIFORM_VECTORS];
    FragmentUniform uniforms[MAX_FRAGMENT_UNIFORM_VECTORS];
    FragmentVarying varyings[MAX_FRAGMENT_VARYING_VARIABLES];
} FragmentState;

//
// Builtin GLSL function decls
// 

extern GLSL_RT_DLLEXPORT vec4 texture2D(unsigned long long samplerID, vec2* coords);
extern GLSL_RT_DLLEXPORT vec4 texture3D(unsigned long long samplerID, vec3* coords);

//
// Builtin LSGL function decls.
//
extern GLSL_RT_DLLEXPORT float shadow(void* frag, vec3* pos, vec3* dir);
extern GLSL_RT_DLLEXPORT float trace(void* frag, vec3* pos, vec3* dir, vec4* col, float* userattrib);
extern GLSL_RT_DLLEXPORT int   rayoption(void* frag, int* prev_double_sided, int double_sided);
//extern int   raydepth();
//extern float isectinfo(vec3* p, vec3* n, vec3* d);

}   // extern "C" 


class ShaderRuntime
{
  public:
    ShaderRuntime();
    ~ShaderRuntime();
};

// Include inline math functions
#include "glsl_math_func.h"

//
// Bridge functions
//
static inline vec4 __glsl_texture2D(Fragment* frag, unsigned long long samplerID, vec2 coords)
{
    vec2 _coord = coords;
    return frag->texture2D(samplerID, &_coord);
}

#if OPTIMIZED_TEXTURING
#define FORCEINLINE __attribute__((always_inline))
#include <algorithm>

FORCEINLINE __vec4_d vclamp(const __vec4_d& v, const __vec4_d& low, const __vec4_d& high) {
  return vmin(vmax(v, low), high);
}

inline double dlerp(double t, double a, double b) { return (1.0 - t) * a + t * b; }

FORCEINLINE __vec4_d vlerp4(const __vec4_d& t, const __vec4_d& a, const __vec4_d& b) {
  const __vec4_d vone(1.0);
  return (vone - t) * a + t * b;
}

FORCEINLINE __vec2_d vlerp2(const __vec2_d& t, const __vec2_d& a, const __vec2_d& b) {
  const __vec2_d vone(1.0);
  return (vone - t) * a + t * b;
}

static inline int fasterfloorf( const float x ) {
  if (x >= 0) {
    return (int)x;
  }

  int y = (int)x;
  if (fabs(x - y) <= 1.0e-6f) {
    // Do nothing.
  } else {
    y = y - 1;
  }

  return y;
}

static inline float Lerp(float t, float a, float b) { return (1.f - t) * a + t * b; }

static inline int Clamp(int v, int low, int high) {
  return std::min(std::max(v, low), high);
}

static inline float ADDR3D(int x, int y, int z, int comp, int idx, int dim[3],
               const float density[]) {
  return density[comp * (z * dim[0] * dim[1] + y * dim[0] + x) + idx];
}

static inline void FilterTex3D(vec4* rgba, const Texture3D* tex, float u, float v, float r) {

  rgba->v[0] = rgba->v[1] = rgba->v[2] = rgba->v[3] = 0.0f;

  if (!tex) return;

  // @fixme { REPEAT only }
  float u01 = u - fasterfloorf(u);
  float v01 = v - fasterfloorf(v);
  float r01 = r - fasterfloorf(r);

  int dim[3];
  dim[0] = tex->width;
  dim[1] = tex->height;
  dim[2] = tex->depth;

  const float *density = tex->image;

  int comps = tex->components;

#if 0

  float vox[3];
  //vox[0] = u01 * dim[0] - .5f;
  //vox[1] = v01 * dim[1] - .5f;
  //vox[2] = r01 * dim[2] - .5f;
  vox[0] = u01 * (dim[0] - 1.0f);
  vox[1] = v01 * (dim[1] - 1.0f);
  vox[2] = r01 * (dim[2] - 1.0f);


  int vx = (int)(vox[0]), vy = (int)(vox[1]), vz = (int)(vox[2]);
  vx = Clamp(vx, 0, dim[0] -1);
  vy = Clamp(vy, 0, dim[1] -1);
  vz = Clamp(vz, 0, dim[2] -1);
  int vx1 = std::min(vx+1, dim[0] - 1);
  int vy1 = std::min(vy+1, dim[1] - 1);
  int vz1 = std::min(vz+1, dim[2] - 1);

  float dx = vox[0] - vx, dy = vox[1] - vy, dz = vox[2] - vz;
  //assert(dx <= 1.0);
  //assert(dy <= 1.0);
  //assert(dz <= 1.0);
  //assert(dx >= 0.0f);
  //assert(dy >= 0.0f);
  //assert(dz >= 0.0f);

  float vals[4] = { 0.0f, 0.0f, 0.0f, 0.0f }; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const float D000 = ADDR3D(vx,  vy , vz , comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy , vz , comps, i, dim, density);
    const float D010 = ADDR3D(vx,  vy1, vz , comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz , comps, i, dim, density);
    const float D001 = ADDR3D(vx,  vy , vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy , vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx,  vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D(vx1, vy1, vz1, comps, i, dim, density);

    // Trilinearly interpolate density values to compute local density
    //float d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
    //                 D(vx + 1, vy, vz, comps, i, dim, density));
    //float d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
    //                 D(vx + 1, vy + 1, vz, comps, i, dim, density));
    //float d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
    //                 D(vx + 1, vy, vz + 1, comps, i, dim, density));
    //float d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
    //                 D(vx + 1, vy + 1, vz + 1, comps, i, dim, density));
    float d00 = Lerp(dx, D000, D100);
    float d10 = Lerp(dx, D010, D110);
    float d01 = Lerp(dx, D001, D101);
    float d11 = Lerp(dx, D011, D111);
    float d0 = Lerp(dy, d00, d10);
    float d1 = Lerp(dy, d01, d11);

    float dlerp = Lerp(dz, d0, d1);
    vals[i] = dlerp;
  }

  if (comps == 1) { // scalar
    rgba->v[0] = vals[0];
    rgba->v[1] = vals[0];
    rgba->v[2] = vals[0];
    rgba->v[3] = vals[0];
  } else {
    rgba->v[0] = vals[0];
    rgba->v[1] = vals[1];
    rgba->v[2] = vals[2];
    rgba->v[3] = vals[3];
  }

#else

  __vec4_d vzero(0.0);
  __vec4_d vone(1.0);
  __vec4_d vdim(dim[0], dim[1], dim[2], 0.0);
  __vec4_d vcoord(u01, v01, r01, 0.0);

  __vec4_d vox = vcoord * (vdim - vone);

  const int vx0 = (int)(vox.u.v[0]);
  const int vy0 = (int)(vox.u.v[1]);
  const int vz0 = (int)(vox.u.v[2]);
  const int vx = Clamp(vx0, 0, dim[0] -1);
  const int vy = Clamp(vy0, 0, dim[1] -1);
  const int vz = Clamp(vz0, 0, dim[2] -1);
  const int vx1 = std::min(vx0+1, dim[0] -1);
  const int vy1 = std::min(vy0+1, dim[1] -1);
  const int vz1 = std::min(vz0+1, dim[2] -1);
  const __vec4_d vxyz0(vx0, vy0, vz0, 0.0);
  const __vec4_d vxyz = vclamp(vxyz0, vzero, vdim - vone);

  double dx = vox.u.v[0] - vxyz.u.v[0];
  double dy = vox.u.v[1] - vxyz.u.v[1];
  double dz = vox.u.v[2] - vxyz.u.v[2];

  double vals[4] = { 0.0, 0.0, 0.0, 0.0 }; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const float D000 = ADDR3D(vx,  vy , vz , comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy , vz , comps, i, dim, density);
    const float D010 = ADDR3D(vx,  vy1, vz , comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz , comps, i, dim, density);
    const float D001 = ADDR3D(vx,  vy , vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy , vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx,  vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D(vx1, vy1, vz1, comps, i, dim, density);

    // (d00, d01, d10, d11)
    const __vec4_d vsrc0(D000, D001, D010, D011);
    const __vec4_d vsrc1(D100, D101, D110, D111);

    __vec4_d vd4 = vlerp4(__vec4_d(dx), vsrc0, vsrc1);
    const __vec2_d vd01 = vlerp2(__vec2_d(dy), vd4.u.v0, vd4.u.v1);
    //const __vec2_d vd01 = __vec2_d(dy) + vd4.u.v0 + vd4.u.v1;

    double dret= dlerp(dz, vd01.u.v[0], vd01.u.v[1]);
    //printf("dret = %f\n", dret);
    vals[i] = dret;
  }

  if (comps == 1) { // scalar
    rgba->v[0] = vals[0];
    rgba->v[1] = vals[0];
    rgba->v[2] = vals[0];
    rgba->v[3] = vals[0];
  } else {
    rgba->v[0] = vals[0];
    rgba->v[1] = vals[1];
    rgba->v[2] = vals[2];
    rgba->v[3] = vals[3];
  }
#endif
}
#endif

static inline vec4 __glsl_texture3D(Fragment* frag, unsigned long long samplerID, vec3 coords)
{
    vec3 _coord = coords;
#if OPTIMIZED_TEXTURING // opt
    const Texture3D* tex = reinterpret_cast<const Texture3D*>(samplerID); 
    vec4 rgba;
    FilterTex3D(&rgba, tex, coords.v[0], coords.v[1], coords.v[2]);
    return rgba;
#else
    return frag->texture3D(samplerID, &_coord);
#endif
}

static inline float __glsl_trace(Fragment* frag, vec3 org, vec3 dir)
{
    vec3 o = org;
    vec3 d = dir;
    return frag->shadow((void*)frag, &o, &d);
}

static inline float __glsl_trace(Fragment* frag, vec3 org, vec3 dir, vec4& col)
{
    vec3 o = org;
    vec3 d = dir;
    vec4 c;             // out
    c.v[0] = c.v[1] = c.v[2] = c.v[3] = 0.0f; // for safety
#if 0
    float tt = trace((void*)frag, &o, &d, &c, NULL);
#else
    float tt = frag->trace((void*)frag, &o, &d, &c, NULL);
#endif
    col = c;
    return tt;
}

static inline float __glsl_trace(Fragment* frag, vec3 org, vec3 dir, vec4& col, float user)
{
    vec3 o = org;
    vec3 d = dir;
    vec4 c;             // out
    c.v[0] = c.v[1] = c.v[2] = c.v[3] = 0.0f; // for safety
    float val = user; // user attribute
#if 0
    float tt = trace((void*)frag, &o, &d, &c, &val);
#else
    float tt = frag->trace((void*)frag, &o, &d, &c, &val);
#endif
    col = c;
    return tt;
}

static inline int __glsl_raydepth(Fragment* frag, int& depth)
{
    depth = frag->raydepth;
    return frag->raydepth;
}

static inline float __glsl_rayattrib(Fragment* frag, float& val)
{
    val = frag->rayattrib;
    return frag->rayattrib;
}

static inline int __glsl_rayoption(Fragment* frag, int& prev_double_sided, int double_sided)
{
  prev_double_sided = frag->doubleSided;
  frag->doubleSided = double_sided;
  return prev_double_sided;
}

static inline float __glsl_isectinfo(Fragment* frag, vec3 &p, vec3& n, vec3& dir)
{
    p.v[0] = frag->position[0];
    p.v[1] = frag->position[1];
    p.v[2] = frag->position[2];

    n.v[0] = frag->normal[0];
    n.v[1] = frag->normal[1];
    n.v[2] = frag->normal[2];

    dir.v[0] = frag->indir[0];
    dir.v[1] = frag->indir[1];
    dir.v[2] = frag->indir[2];

    return frag->hitdist;
}

static inline int __glsl_numIntersects(Fragment* frag, int& n)
{
    // @todo
    n = 1;
    return n;
}

static inline float __glsl_queryIntersect(Fragment* frag, int i, vec3& p, vec3& normal, vec3& geom_normal, vec3& tangent, vec3& binormal, vec3& indir, vec2& barycentric_coord)
{
    // @fixe { look up i'th intersection point. }

    p.v[0] = frag->position[0];
    p.v[1] = frag->position[1];
    p.v[2] = frag->position[2];

    normal.v[0] = frag->normal[0];
    normal.v[1] = frag->normal[1];
    normal.v[2] = frag->normal[2];

    geom_normal.v[0] = frag->geometricNormal[0];
    geom_normal.v[1] = frag->geometricNormal[1];
    geom_normal.v[2] = frag->geometricNormal[2];

    tangent.v[0] = frag->tangent[0];
    tangent.v[1] = frag->tangent[1];
    tangent.v[2] = frag->tangent[2];

    binormal.v[0] = frag->binormal[0];
    binormal.v[1] = frag->binormal[1];
    binormal.v[2] = frag->binormal[2];

    indir.v[0] = frag->indir[0];
    indir.v[1] = frag->indir[1];
    indir.v[2] = frag->indir[2];

    barycentric_coord.v[0] = frag->barycentric[0];
    barycentric_coord.v[1] = frag->barycentric[1];

    return frag->hitdist;
}

static inline float __glsl_camerainfo(Fragment* frag, vec3 &eye, vec3& lookat, vec3& up)
{
    eye.v[0] = frag->cameraFrame[0][0];
    eye.v[1] = frag->cameraFrame[0][1];
    eye.v[2] = frag->cameraFrame[0][2];

    lookat.v[0] = frag->cameraFrame[1][0];
    lookat.v[1] = frag->cameraFrame[1][1];
    lookat.v[2] = frag->cameraFrame[1][2];

    up.v[0] = frag->cameraFrame[2][0];
    up.v[1] = frag->cameraFrame[2][1];
    up.v[2] = frag->cameraFrame[2][2];

    return frag->cameraFov;
}

// To prevent glsl compiler optimize away random() function, return random value as an 'out' variable.
static inline float __glsl_random(Fragment* frag, float& value)
{
  float ret = frag->random(frag->threadID);
  value = ret;
  return ret;
}

static inline void __glsl_discard(Fragment* frag)
{
  frag->fragDiscarded = 1;
  return;
}

//
// Util functions
//
static inline vec2 __make_vec2(float x, float y) {
    vec2 v;
    v.v[0] = x; 
    v.v[1] = y; 
    return v;
}

static inline vec3 __make_vec3(float x, float y, float z) {
    vec3 v;
    v.v[0] = x; 
    v.v[1] = y; 
    v.v[2] = z; 
    return v;
}

static inline vec4 __make_vec4(float x, float y, float z, float w) {
    vec4 v;
    v.v[0] = x; 
    v.v[1] = y; 
    v.v[2] = z; 
    v.v[3] = w; 
    return v;
}

static inline mat2 __make_mat2(
  float x0, float y0,
  float x1, float y1)
{
    mat2 m;
    m.v[0][0] = x0; 
    m.v[0][1] = y0; 
    m.v[1][0] = x1; 
    m.v[1][1] = y1; 
    return m;
}

static inline mat3 __make_mat3(
  float x0, float y0, float z0,
  float x1, float y1, float z1,
  float x2, float y2, float z2)
{
    mat3 m;
    m.v[0][0] = x0; 
    m.v[0][1] = y0; 
    m.v[0][2] = z0; 
    m.v[1][0] = x1; 
    m.v[1][1] = y1; 
    m.v[1][2] = z1; 
    m.v[2][0] = x2; 
    m.v[2][1] = y2; 
    m.v[2][2] = z2; 
    return m;
}

static inline mat4 __make_mat4(
  float x0, float y0, float z0, float w0,
  float x1, float y1, float z1, float w1,
  float x2, float y2, float z2, float w2,
  float x3, float y3, float z3, float w3)
{
    mat4 m;
    m.v[0][0] = x0; 
    m.v[0][1] = y0; 
    m.v[0][2] = z0; 
    m.v[0][3] = w0; 
    m.v[1][0] = x1; 
    m.v[1][1] = y1; 
    m.v[1][2] = z1; 
    m.v[1][3] = w1; 
    m.v[2][0] = x2; 
    m.v[2][1] = y2; 
    m.v[2][2] = z2; 
    m.v[2][3] = w2; 
    m.v[3][0] = x3; 
    m.v[3][1] = y3; 
    m.v[3][2] = z3; 
    m.v[3][3] = w3; 
    return m;
}

static inline float __swizzle(int x, vec2 v)
{
    return v.v[x];
}

static inline float __swizzle(int x, vec3 v)
{
    return v.v[x];
}

static inline float __swizzle(int x, vec4 v)
{
    return v.v[x];
}

static inline vec2 __swizzle(int x, int y, float f)
{
    vec2 v2;
    v2.v[0] = f;
    v2.v[1] = f;
    return v2;
}

static inline vec2 __swizzle(int x, int y, vec2 v)
{
    vec2 v2;
    v2.v[0] = v.v[x];
    v2.v[1] = v.v[y];
    return v2;
}

static inline vec2 __swizzle(int x, int y, vec3 v)
{
    vec2 v2;
    v2.v[0] = v.v[x];
    v2.v[1] = v.v[y];
    return v2;
}

static inline vec2 __swizzle(int x, int y, vec4 v)
{
    vec2 v2;
    v2.v[0] = v.v[x];
    v2.v[1] = v.v[y];
    return v2;
}

static inline vec3 __swizzle(int x, int y, int z, float f)
{
    vec3 v3;
    v3.v[0] = f;
    v3.v[1] = f;
    v3.v[2] = f;
    return v3;
}

static inline vec3 __swizzle(int x, int y, int z, vec3 v)
{
    vec3 v3;
    v3.v[0] = v.v[x];
    v3.v[1] = v.v[y];
    v3.v[2] = v.v[z];
    return v3;
}

static inline vec3 __swizzle(int x, int y, int z, vec4 v)
{
    vec3 v3;
    v3.v[0] = v.v[x];
    v3.v[1] = v.v[y];
    v3.v[2] = v.v[z];
    return v3;
}

static inline float __swizzle(int x, float v[4])
{
    return v[x];
}

static inline vec2 __swizzle(int x, int y, float v[4])
{
    vec2 v2;
    v2.v[0] = v[x];
    v2.v[1] = v[y];
    return v2;
}

static inline vec3 __swizzle(int x, int y, int z, float v[4])
{
    vec3 v3;
    v3.v[0] = v[x];
    v3.v[1] = v[y];
    v3.v[2] = v[z];
    return v3;
}

static inline vec4 __swizzle(int x, int y, int z, int w, float f)
{
    vec4 v4;
    v4.v[0] = f;
    v4.v[1] = f;
    v4.v[2] = f;
    v4.v[3] = f;
    return v4;
}

static inline vec4 __swizzle(int x, int y, int z, int w, float v[4])
{
    vec4 v4;
    v4.v[0] = v[x];
    v4.v[1] = v[y];
    v4.v[2] = v[z];
    v4.v[3] = v[w];
    return v4;
}

static inline void __assign(float a[4], vec4 b)
{
    a[0] = b.v[0];
    a[1] = b.v[1];
    a[2] = b.v[2];
    a[3] = b.v[3];
}

static inline void __assign(vec4& a, vec4 b)
{
    a = b;
}

static inline float __add(float a, float b)
{
    return a + b;
}

static inline vec2 __add(vec2 a, vec2 b)
{
    vec2 c;
    c.v[0] = a.v[0] + b.v[0];
    c.v[1] = a.v[1] + b.v[1];
    return c;
}

static inline vec2 __add(vec2 a, float b)
{
    vec2 c;
    c.v[0] = a.v[0] + b;
    c.v[1] = a.v[1] + b;
    return c;
}

static inline vec2 __add(float b, vec2 a)
{
    vec2 c;
    c.v[0] = a.v[0] + b;
    c.v[1] = a.v[1] + b;
    return c;
}

static inline vec3 __add(vec3 a, vec3 b)
{
    vec3 c;
    c.v[0] = a.v[0] + b.v[0];
    c.v[1] = a.v[1] + b.v[1];
    c.v[2] = a.v[2] + b.v[2];
    return c;
}

static inline vec3 __add(vec3 a, float b)
{
    vec3 c;
    c.v[0] = a.v[0] + b;
    c.v[1] = a.v[1] + b;
    c.v[2] = a.v[2] + b;
    return c;
}

static inline vec3 __add(float b, vec3 a)
{
    vec3 c;
    c.v[0] = a.v[0] + b;
    c.v[1] = a.v[1] + b;
    c.v[2] = a.v[2] + b;
    return c;
}

static inline vec4 __add(vec4 a, vec4 b)
{
    vec4 c;
    c.v[0] = a.v[0] + b.v[0];
    c.v[1] = a.v[1] + b.v[1];
    c.v[2] = a.v[2] + b.v[2];
    c.v[3] = a.v[3] + b.v[3];
    return c;
}

static inline vec4 __add(vec4 a, float b)
{
    vec4 c;
    c.v[0] = a.v[0] + b;
    c.v[1] = a.v[1] + b;
    c.v[2] = a.v[2] + b;
    c.v[3] = a.v[3] + b;
    return c;
}

static inline vec4 __add(float b, vec4 a)
{
    vec4 c;
    c.v[0] = a.v[0] + b;
    c.v[1] = a.v[1] + b;
    c.v[2] = a.v[2] + b;
    c.v[3] = a.v[3] + b;
    return c;
}

static inline float __sub(float a, float b)
{
    return a - b;
}

static inline vec2 __sub(vec2 a, vec2 b)
{
    vec2 c;
    c.v[0] = a.v[0] - b.v[0];
    c.v[1] = a.v[1] - b.v[1];
    return c;
}

static inline vec2 __sub(vec2 a, float b)
{
    vec2 c;
    c.v[0] = a.v[0] - b;
    c.v[1] = a.v[1] - b;
    return c;
}

static inline vec2 __sub(float b, vec2 a)
{
    vec2 c;
    c.v[0] = b - a.v[0];
    c.v[1] = b - a.v[1];
    return c;
}

static inline vec3 __sub(vec3 a, vec3 b)
{
    vec3 c;
    c.v[0] = a.v[0] - b.v[0];
    c.v[1] = a.v[1] - b.v[1];
    c.v[2] = a.v[2] - b.v[2];
    return c;
}

static inline vec3 __sub(vec3 a, float b)
{
    vec3 c;
    c.v[0] = a.v[0] - b;
    c.v[1] = a.v[1] - b;
    c.v[2] = a.v[2] - b;
    return c;
}

static inline vec3 __sub(float b, vec3 a)
{
    vec3 c;
    c.v[0] = b - a.v[0];
    c.v[1] = b - a.v[1];
    c.v[2] = b - a.v[2];
    return c;
}

static inline vec4 __sub(vec4 a, vec4 b)
{
    vec4 c;
    c.v[0] = a.v[0] - b.v[0];
    c.v[1] = a.v[1] - b.v[1];
    c.v[2] = a.v[2] - b.v[2];
    c.v[3] = a.v[3] - b.v[3];
    return c;
}

static inline vec4 __sub(vec4 a, float b)
{
    vec4 c;
    c.v[0] = a.v[0] - b;
    c.v[1] = a.v[1] - b;
    c.v[2] = a.v[2] - b;
    c.v[3] = a.v[3] - b;
    return c;
}

static inline vec4 __sub(float b, vec4 a)
{
    vec4 c;
    c.v[0] = b - a.v[0];
    c.v[1] = b - a.v[1];
    c.v[2] = b - a.v[2];
    c.v[3] = b - a.v[3];
    return c;
}

static inline float __mul(float a, float b)
{
    return a * b;
}

static inline vec2 __mul(vec2 a, vec2 b)
{
    vec2 c;
    c.v[0] = a.v[0] * b.v[0];
    c.v[1] = a.v[1] * b.v[1];
    return c;
}

static inline vec2 __mul(vec2 a, float b)
{
    vec2 c;
    c.v[0] = a.v[0] * b;
    c.v[1] = a.v[1] * b;
    return c;
}

static inline vec2 __mul(float b, vec2 a)
{
    vec2 c;
    c.v[0] = a.v[0] * b;
    c.v[1] = a.v[1] * b;
    return c;
}

static inline vec3 __mul(vec3 a, vec3 b)
{
    vec3 c;
    c.v[0] = a.v[0] * b.v[0];
    c.v[1] = a.v[1] * b.v[1];
    c.v[2] = a.v[2] * b.v[2];
    return c;
}

static inline mat3 __mul(mat3 a, mat3 b)
{
    mat3 c;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            c.v[i][j] = 0.;
            for (int k = 0; k < 3; ++k) {
                c.v[i][j] += a.v[i][k] * b.v[k][j];
            }
        }
    }
    return c;
}

static inline vec3 __mul(vec3 a, float b)
{
    vec3 c;
    c.v[0] = a.v[0] * b;
    c.v[1] = a.v[1] * b;
    c.v[2] = a.v[2] * b;
    return c;
}

static inline vec3 __mul(float b, vec3 a)
{
    vec3 c;
    c.v[0] = a.v[0] * b;
    c.v[1] = a.v[1] * b;
    c.v[2] = a.v[2] * b;
    return c;
}

static inline vec4 __mul(vec4 a, vec4 b)
{
    vec4 c;
    c.v[0] = a.v[0] * b.v[0];
    c.v[1] = a.v[1] * b.v[1];
    c.v[2] = a.v[2] * b.v[2];
    c.v[3] = a.v[3] * b.v[3];
    return c;
}

static inline vec4 __mul(vec4 a, float b)
{
    vec4 c;
    c.v[0] = a.v[0] * b;
    c.v[1] = a.v[1] * b;
    c.v[2] = a.v[2] * b;
    c.v[3] = a.v[3] * b;
    return c;
}

static inline vec4 __mul(float b, vec4 a)
{
    vec4 c;
    c.v[0] = a.v[0] * b;
    c.v[1] = a.v[1] * b;
    c.v[2] = a.v[2] * b;
    c.v[3] = a.v[3] * b;
    return c;
}

static inline vec3 __mul(mat3 a, vec3 b)
{
    vec3 c;
    c.v[0] = a.v[0][0] * b.v[0] + a.v[1][0] * b.v[1] + a.v[2][0] * b.v[2];
    c.v[1] = a.v[0][1] * b.v[0] + a.v[1][1] * b.v[1] + a.v[2][1] * b.v[2];
    c.v[2] = a.v[0][2] * b.v[0] + a.v[1][2] * b.v[1] + a.v[2][2] * b.v[2];
    return c;
}

static inline vec4 __mul(mat4 a, vec4 b)
{
    vec4 c;
    c.v[0] = a.v[0][0] * b.v[0] + a.v[1][0] * b.v[1] + a.v[2][0] * b.v[2] + a.v[3][0] * b.v[3];
    c.v[1] = a.v[0][1] * b.v[0] + a.v[1][1] * b.v[1] + a.v[2][1] * b.v[2] + a.v[3][1] * b.v[3];
    c.v[2] = a.v[0][2] * b.v[0] + a.v[1][2] * b.v[1] + a.v[2][2] * b.v[2] + a.v[3][2] * b.v[3];
    c.v[3] = a.v[0][3] * b.v[0] + a.v[1][3] * b.v[1] + a.v[2][3] * b.v[2] + a.v[3][3] * b.v[3];
    return c;
}

static inline mat4 __mul(mat4 a, float b)
{
  mat4 c;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      c.v[i][j] = a.v[i][j] * b;
    }
  }
  return c;
}

static inline float __div(float a, float b)
{
    return a / b;
}

static inline vec2 __div(vec2 a, vec2 b)
{
    vec2 c;
    c.v[0] = a.v[0] / b.v[0];
    c.v[1] = a.v[1] / b.v[1];
    return c;
}

static inline vec2 __div(vec2 a, float b)
{
    vec2 c;
    c.v[0] = a.v[0] / b;
    c.v[1] = a.v[1] / b;
    return c;
}

static inline vec2 __div(float b, vec2 a)
{
    vec2 c;
    c.v[0] = b / a.v[0];
    c.v[1] = b / a.v[1];
    return c;
}

static inline vec3 __div(vec3 a, vec3 b)
{
    vec3 c;
    c.v[0] = a.v[0] / b.v[0];
    c.v[1] = a.v[1] / b.v[1];
    c.v[2] = a.v[2] / b.v[2];
    return c;
}

static inline vec3 __div(vec3 a, float b)
{
    vec3 c;
    c.v[0] = a.v[0] / b;
    c.v[1] = a.v[1] / b;
    c.v[2] = a.v[2] / b;
    return c;
}

static inline vec3 __div(float b, vec3 a)
{
    vec3 c;
    c.v[0] = b / a.v[0];
    c.v[1] = b / a.v[1];
    c.v[2] = b / a.v[2];
    return c;
}

static inline vec4 __div(vec4 a, vec4 b)
{
    vec4 c;
    c.v[0] = a.v[0] / b.v[0];
    c.v[1] = a.v[1] / b.v[1];
    c.v[2] = a.v[2] / b.v[2];
    c.v[3] = a.v[3] / b.v[3];
    return c;
}

static inline vec4 __div(vec4 a, float b)
{
    vec4 c;
    c.v[0] = a.v[0] / b;
    c.v[1] = a.v[1] / b;
    c.v[2] = a.v[2] / b;
    c.v[3] = a.v[3] / b;
    return c;
}

static inline vec4 __div(float b, vec4 a)
{
    vec4 c;
    c.v[0] = b / a.v[0];
    c.v[1] = b / a.v[1];
    c.v[2] = b / a.v[2];
    c.v[3] = b / a.v[3];
    return c;
}

static inline float __rcp(float a)
{
    return 1.0f / a;
}

static inline vec2 __rcp(vec2 a)
{
    vec2 c;
    c.v[0] = 1.0f / a.v[0];
    c.v[1] = 1.0f / a.v[1];
    return c;
}

static inline vec3 __rcp(vec3 a)
{
    vec3 c;
    c.v[0] = 1.0f / a.v[0];
    c.v[1] = 1.0f / a.v[1];
    c.v[2] = 1.0f / a.v[2];
    return c;
}

static inline vec4 __rcp(vec4 a)
{
    vec4 c;
    c.v[0] = 1.0f / a.v[0];
    c.v[1] = 1.0f / a.v[1];
    c.v[2] = 1.0f / a.v[2];
    c.v[3] = 1.0f / a.v[3];
    return c;
}

static inline float __rsq(float a)
{
    return 1.0f / sqrtf(a);
}

static inline vec2 __rsq(vec2 a)
{
    vec2 c;
    c.v[0] = 1.0f / sqrtf(a.v[0]);
    c.v[1] = 1.0f / sqrtf(a.v[1]);
    return c;
}

static inline vec3 __rsq(vec3 a)
{
    vec3 c;
    c.v[0] = 1.0f / sqrtf(a.v[0]);
    c.v[1] = 1.0f / sqrtf(a.v[1]);
    c.v[2] = 1.0f / sqrtf(a.v[2]);
    return c;
}

static inline vec4 __rsq(vec4 a)
{
    vec4 c;
    c.v[0] = 1.0f / sqrtf(a.v[0]);
    c.v[1] = 1.0f / sqrtf(a.v[1]);
    c.v[2] = 1.0f / sqrtf(a.v[2]);
    c.v[3] = 1.0f / sqrtf(a.v[3]);
    return c;
}

static inline float __neg(float a)
{
    return -a;
}

static inline vec2 __neg(vec2 a)
{
    vec2 c;
    c.v[0] = -a.v[0];
    c.v[1] = -a.v[1];
    return c;
}

static inline vec3 __neg(vec3 a)
{
    vec3 c;
    c.v[0] = -a.v[0];
    c.v[1] = -a.v[1];
    c.v[2] = -a.v[2];
    return c;
}

static inline vec4 __neg(vec4 a)
{
    vec4 c;
    c.v[0] = -a.v[0];
    c.v[1] = -a.v[1];
    c.v[2] = -a.v[2];
    c.v[3] = -a.v[3];
    return c;
}

static inline float __i2f(int a)
{
    return (float)a;
}

static inline vec2 __i2f(ivec2 a)
{
    vec2 c;
    c.v[0] = (float)a.v[0];
    c.v[1] = (float)a.v[1];
    return c;
}

static inline vec3 __i2f(ivec3 a)
{
    vec3 c;
    c.v[0] = (float)a.v[0];
    c.v[1] = (float)a.v[1];
    c.v[2] = (float)a.v[2];
    return c;
}

static inline vec4 __i2f(ivec4 a)
{
    vec4 c;
    c.v[0] = (float)a.v[0];
    c.v[1] = (float)a.v[1];
    c.v[2] = (float)a.v[2];
    c.v[3] = (float)a.v[3];
    return c;
}

static inline float __b2f(bool a)
{
    return (float)a;
}

static inline vec2 __b2f(bvec2 a)
{
    vec2 c;
    c.v[0] = (float)a.v[0];
    c.v[1] = (float)a.v[1];
    return c;
}

static inline vec3 __b2f(bvec3 a)
{
    vec3 c;
    c.v[0] = (float)a.v[0];
    c.v[1] = (float)a.v[1];
    c.v[2] = (float)a.v[2];
    return c;
}

static inline vec4 __b2f(bvec4 a)
{
    vec4 c;
    c.v[0] = (float)a.v[0];
    c.v[1] = (float)a.v[1];
    c.v[2] = (float)a.v[2];
    c.v[3] = (float)a.v[3];
    return c;
}

static inline int __f2i(float a)
{
    return (int)a;
}

static inline ivec2 __f2i(vec2 a)
{
    ivec2 c;
    c.v[0] = (int)a.v[0];
    c.v[1] = (int)a.v[1];
    return c;
}

static inline ivec3 __f2i(vec3 a)
{
    ivec3 c;
    c.v[0] = (int)a.v[0];
    c.v[1] = (int)a.v[1];
    c.v[2] = (int)a.v[2];
    return c;
}

static inline ivec4 __f2i(vec4 a)
{
    ivec4 c;
    c.v[0] = (int)a.v[0];
    c.v[1] = (int)a.v[1];
    c.v[2] = (int)a.v[2];
    c.v[3] = (int)a.v[3];
    return c;
}

static inline bool __and(bool a, bool b)
{
    return a && b;
}

static inline bool __or(bool a, bool b)
{
    return a || b;
}

static inline bool __eq(bool a, bool b)
{
    return a == b;
}

static inline bool __neq(bool a, bool b)
{
    return a != b;
}

static inline bool __all_eq(bool a, bool b)
{
    return a == b;
}

static inline bool __all_eq(int a, int b)
{
    return a == b;
}

static inline bool __all_eq(float a, float b)
{
    return a == b; // @todo { use float eq. }
}


static inline bool __any_neq(bool a, bool b)
{
    return a != b;
}

static inline bool __any_neq(int a, int b)
{
    return a != b;
}

static inline bool __any_neq(float a, float b)
{
    return a != b; // @todo { use float neq }
}

static inline bool __gt(int a, int b)
{
    return a > b;
}

static inline bool __gt(float a, float b)
{
    return a > b;
}

static inline bool __ge(int a, int b)
{
    return a >= b;
}


static inline bool __ge(float a, float b)
{
    return a >= b;
}

static inline bool __lt(int a, int b)
{
    return a < b;
}

static inline bool __lt(float a, float b)
{
    return a < b;
}

static inline bool __le(int a, int b)
{
    return a <= b;
}

static inline bool __le(float a, float b)
{
    return a <= b;
}

static inline bool __not(int a)
{
    return !a;
}

static inline bool __not(float a)
{
    return !a; // ???
}

#endif  // __GLSL_CC_SHADER_RUNTIME_H__
