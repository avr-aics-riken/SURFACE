/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2016 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <limits>
#include <functional>
#include <algorithm>

#include "render_accel_solid.h"
#include "render_prefix_tree_util.h"

using namespace lsgl::render;

#ifdef __FUJITSU
#define FORCEINLINE __attribute__((always_inline))
#else
#define FORCEINLINE inline
#endif

#define MAX_LEAF_ELEMENTS (4)
#define MAX_TREE_DEPTH_32BIT                                                   \
  (22) // FYI, log2(1G/16) ~= 25.897, log2(1G/32) ~= 21

#define ENABLE_TRACE_PRINT (0)
#define ENABLE_DEBUG_PRINT (0)

#define ENABLE_TRAVERSAL_STATISTICS (0)

#define trace(f, ...)                                                          \
  {                                                                            \
    if (ENABLE_TRACE_PRINT)                                                    \
      printf(f, __VA_ARGS__);                                                  \
  }

#if ENABLE_DEBUG_PRINT
#define debug(f, ...)                                                          \
  { printf(f, __VA_ARGS__); }
#else
#define debug(f, ...)
#endif

namespace {

// assume real == float
void PrintReal(std::string msg, real v) {
  printf("%s : %f(0x%08x)\n", msg.c_str(), v, *((unsigned int *)&v));
}

void PrintVec3(std::string msg, real3 v) {
  printf("%s : %f(0x%08x), %f(0x%08x), %f(0x%08x)\n", msg.c_str(), v[0],
         *((unsigned int *)&v[0]), v[1], *((unsigned int *)&v[1]), v[2],
         *((unsigned int *)&v[2]));
}

FORCEINLINE double3 vcrossd(double3 a, double3 b) {
  double3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

FORCEINLINE double vdotd(double3 a, double3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

FORCEINLINE double length(double3 a) {
  return sqrt(a.x() * a.x() + a.y() * a.y() + a.z() * a.z());
}

FORCEINLINE double3 normalize(double3 a) { return a / length(a); }

FORCEINLINE real3 toreal3(double3 a) { return real3(a.x(), a.y(), a.z()); }

FORCEINLINE int fsign(const double x) {
  // real eps = std::numeric_limits<real>::epsilon() * 16;
  double eps = DBL_EPSILON;
  return (x > eps ? 1 : (x < -eps ? -1 : 0));
}

// Vertex order definitial table. Counter clock-wise.
const int kTetraFaces[4][3] = {{0,2,1}, {1,2,3}, {0,3,2}, {0,1,3} };
const int kPyramidFaces[5][4] = {{0,3,2,1}, {0,1,4,-1}, {1,2,4,-1}, {2,3,4,-1}, {0,4,3,-1} };
const int kPrismFaces[5][4] = {{0,2,1,-1}, {3,4,5,-1}, {0,3,5,2}, {0,1,4,3}, {1,2,5,4} };
const int kHexaFaces[6][4] = {{0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {0,4,7,3}, };
  
//
// Simple Pluecker coordinate class
//
class Pluecker {
public:
  double3 d; // direction
  double3 c; // cross
  
  Pluecker(const double3 &v0, const double3 &v1)
  : d(v1 - v0), c(vcrossd(v1, v0)) {}
};
  
// Inner product
FORCEINLINE double operator*(const Pluecker &p0, const Pluecker &p1) {
  return vdotd(p0.d, p1.c) + vdotd(p1.d, p0.c);
}
  
// Up to 12 edges(Hexahedron)
void GetEdges(double3 *edges, int solidType, const double3 *vertices)
{
  if (solidType == 5) { // Pyramid
    edges[0] = vertices[1]-vertices[0];
    edges[1] = vertices[2]-vertices[1];
    edges[2] = vertices[3]-vertices[2];
    edges[3] = vertices[0]-vertices[3];
    edges[4] = vertices[4]-vertices[0];
    edges[5] = vertices[4]-vertices[1];
    edges[6] = vertices[4]-vertices[2];
    edges[7] = vertices[4]-vertices[3];
  } else if (solidType == 6) { // Prism
    edges[0] = vertices[1]-vertices[0];
    edges[1] = vertices[2]-vertices[1];
    edges[2] = vertices[0]-vertices[2];
    edges[3] = vertices[4]-vertices[3];
    edges[4] = vertices[5]-vertices[4];
    edges[5] = vertices[3]-vertices[5];
    edges[6] = vertices[3]-vertices[0];
    edges[7] = vertices[4]-vertices[1];
    edges[8] = vertices[5]-vertices[2];
  } else if (solidType == 8) { // Hexa
    edges[0] = vertices[1]-vertices[0];
    edges[1] = vertices[2]-vertices[1];
    edges[2] = vertices[3]-vertices[2];
    edges[3] = vertices[0]-vertices[3];
    edges[4] = vertices[5]-vertices[4];
    edges[5] = vertices[6]-vertices[5];
    edges[6] = vertices[7]-vertices[6];
    edges[7] = vertices[4]-vertices[7];
    edges[8] = vertices[4]-vertices[0];
    edges[9] = vertices[5]-vertices[1];
    edges[10] = vertices[6]-vertices[2];
    edges[11] = vertices[7]-vertices[3];
  }
}
  
void GetEdges(real3 *edges, int solidType, const real3 *vertices)
{
  if (solidType == 5) { // Pyramid
    edges[0] = vertices[1]-vertices[0];
    edges[1] = vertices[2]-vertices[1];
    edges[2] = vertices[3]-vertices[2];
    edges[3] = vertices[0]-vertices[3];
    edges[4] = vertices[4]-vertices[0];
    edges[5] = vertices[4]-vertices[1];
    edges[6] = vertices[4]-vertices[2];
    edges[7] = vertices[4]-vertices[3];
  } else if (solidType == 6) { // Prism
    edges[0] = vertices[1]-vertices[0];
    edges[1] = vertices[2]-vertices[1];
    edges[2] = vertices[0]-vertices[2];
    edges[3] = vertices[4]-vertices[3];
    edges[4] = vertices[5]-vertices[4];
    edges[5] = vertices[3]-vertices[5];
    edges[6] = vertices[3]-vertices[0];
    edges[7] = vertices[4]-vertices[1];
    edges[8] = vertices[5]-vertices[2];
  } else if (solidType == 8) { // Hexa
    edges[0] = vertices[1]-vertices[0];
    edges[1] = vertices[2]-vertices[1];
    edges[2] = vertices[3]-vertices[2];
    edges[3] = vertices[0]-vertices[3];
    edges[4] = vertices[5]-vertices[4];
    edges[5] = vertices[6]-vertices[5];
    edges[6] = vertices[7]-vertices[6];
    edges[7] = vertices[4]-vertices[7];
    edges[8] = vertices[4]-vertices[0];
    edges[9] = vertices[5]-vertices[1];
    edges[10] = vertices[6]-vertices[2];
    edges[11] = vertices[7]-vertices[3];
  }
}
  
  inline bool contain(const int &value, const int list[4]){
    if (list[0] == value) return true;
    if (list[1] == value) return true;
    if (list[2] == value) return true;
    if (list[3] == value) return true;
    return false;
  }
  
  inline int edge_n (const int &type){
    switch (type) {
      case 4:
        return 6;
      case 5:
        return 8;
      case 6:
        return 9;
      case 8:
        return 12;
      default:
        break;
    }
    return 0;
  }
  
  void ComputeEdgePcs (real3 *OutPc, const real3 *vertices,
                       const real3 *edges, const int solid_type){
    
    int xxx = 1;
    switch (solid_type) {
      case 4:
        xxx = 3;
        break;
      case 5:
        xxx = 4;
        break;
      case 6:
        xxx = 6;
        break;
      case 8:
        xxx = 8;
        break;
    }
    
    for (int i = 0; i < edge_n(solid_type); i++)
      OutPc[i] = cross(edges[i], vertices[i%xxx]);
  }
  
  inline real3 getCP_sq(const real3 &v0, const real3 &v1, const real3 &v2,
                        const double &ws0, const double &ws1, const double &ws2, const double &ws3,
                        const real3 &edge2, const real3 &edge3, const real3 &raydir){
    double w = ws2 + ws3 + dot(cross(edge2, edge3),raydir);
    return (v2*ws0 + v0*ws1 + v1*w) / (ws0+ws1+w);
  }
  
void interpolate(float d_out[12], float p[3], const int solid_type, const double3 vertices[12], int interpolate_mode = 0) {
  static const int faces[][6][4] = {
    {{}}, {{}}, {{}}, {{}},        // 0,1,2,3
    {{0,2,1,-1}, {1,2,3,-1}, {0,3,2,-1}, {0,1,3,-1} },  //tetra
    {{0,3,2,1}, {0,1,4,-1}, {1,2,4,-1}, {2,3,4,-1}, {0,4,3,-1} },  //pyramid
    {{0,2,1,-1}, {3,4,5,-1}, {0,3,5,2}, {0,1,4,3}, {1,2,5,4} },   //prism
    {{}},       // 7
    {{0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {0,4,7,3}, },  //hexa
  };
  
  real3 pt (p);
  real3 verts[8];
  for (int i = 0; i < solid_type; i++) {
    verts[i] = toreal3(vertices[i]);
  }
  
  pt = pt - verts[0];
  for (int i = solid_type - 1; i >= 0; i--) {
    verts[i] = verts[i] - verts[0];
  }
  
  real3 edges[16];
  GetEdges(edges,solid_type, verts);
  
  float kEPS;     //@fixme
  for (int i = 0; i < edge_n(solid_type); i++)
    kEPS += edges[i].length();
  
  kEPS = (kEPS * 1e-6) / edge_n(solid_type);
    
  bool cw_f[2][16]; // Alloc enough size. edges.size()
  float ws[16]; // Alloc enough size. edges.size()
  
  real3 edgepc[16]; // Alloc enough size. edges.size()
  ComputeEdgePcs(edgepc, verts, edges, solid_type);
  
  for(int v = 0; v < solid_type; v++){
    
    real3 rayorg (verts[v]);
    real3 raydir (pt - rayorg);
    real3 raypc (cross(raydir, rayorg));
    
    float d = (pt - rayorg).length();
    
    // cross flag of ray vs edges.
    bool cross_rve[edge_n(solid_type)];
    
    for(int i = 0; i < edge_n(solid_type) ; i++){
      ws[i] = dot(raydir, edgepc[i]) + dot(edges[i], raypc);
      if(ws[i] >= -kEPS)  cw_f[0][i] = true;
      else cw_f[0][i] = false;
      if(ws[i] <= kEPS)  cw_f[1][i] = true;
      else cw_f[1][i] = false;
      if (cw_f[0][i] && cw_f[1][i]) cross_rve[i] = 1;
    }
    
    real3 cp(1e16,1e16,1e16);
    
    switch(solid_type){
      case 4:     //tetra
        if( cw_f[0][0] && cw_f[0][1] && cw_f[0][2] && !contain(v, faces[solid_type][0]) )
          cp = (verts[2]*ws[0] + verts[0]*ws[1] + verts[1]*ws[2]) / (ws[0]+ws[1]+ws[2]);
        else if( cw_f[1][1] && cw_f[0][4] && cw_f[1][5] && !contain(v, faces[solid_type][1]) )
          cp = (verts[3]*ws[1] + verts[1]*ws[5] - verts[2]*ws[4]) / (ws[5]+ws[1]-ws[4]);
        else if( cw_f[1][2] && cw_f[1][3] && cw_f[0][5] && !contain(v, faces[solid_type][2]) )
          cp = (verts[3]*ws[2] + verts[2]*ws[3] - verts[0]*ws[5]) / (ws[2]+ws[3]-ws[5]);
        else if( cw_f[1][0] && cw_f[0][3] && cw_f[1][4] && !contain(v, faces[solid_type][3]) )
          cp = (verts[3]*ws[0] + verts[0]*ws[4] - verts[1]*ws[3]) / (ws[0]+ws[4]-ws[3]);
        else {
          d_out[v] = 0;
          continue;
        }
        break;
      case 5:     //Pyramid
        if(cw_f[1][0] && cw_f[1][1] && cw_f[1][2] && cw_f[1][3] && !contain(v, faces[solid_type][0])){
          cp = getCP_sq(verts[0], verts[1], verts[2], ws[0], ws[1], ws[2], ws[3],
                        edges[2], edges[3], raydir);
        }else if (cw_f[0][0] && cw_f[1][4] && cw_f[0][5] && !contain(v, faces[solid_type][1])){
          cp = (verts[4]*-ws[0] + verts[1]*ws[4] + verts[0]*-ws[5]) / (-ws[0]+ws[4]-ws[5]);
        }else if (cw_f[0][1] && cw_f[1][5] && cw_f[0][6] && !contain(v, faces[solid_type][2])){
          cp = (verts[4]*-ws[1] + verts[2]*ws[5] + verts[1]*-ws[6]) / (-ws[1]+ws[5]-ws[6]);
        }else if (cw_f[0][2] && cw_f[1][6] && cw_f[0][7] && !contain(v, faces[solid_type][3])){
          cp = (verts[4]*-ws[2] + verts[3]*ws[6] + verts[2]*-ws[7]) / (-ws[2]+ws[6]-ws[7]);
        }else if (cw_f[0][3] && cw_f[1][7] && cw_f[0][4] && !contain(v, faces[solid_type][4])){
          cp = (verts[4]*-ws[3] + verts[0]*ws[7] + verts[3]*-ws[4]) / (-ws[3]+ws[7]-ws[4]);
        }else {
          d_out[v] = 0;
          continue;
        }
        break;
      case 6:     //Prism
        if(cw_f[0][0] && cw_f[0][1] && cw_f[0][2] && !contain(v, faces[solid_type][0])){
            cp = (verts[2]*ws[0] + verts[0]*ws[1] + verts[1]*ws[2]) / (ws[0]+ws[1]+ws[2]);
        } else if(cw_f[1][3] && cw_f[1][4] && cw_f[1][5] &&
                  !contain(v, faces[solid_type][1])){
            cp = (verts[5]*ws[3] + verts[3]*ws[4] + verts[4]*ws[5]) / (ws[3]+ws[4]+ws[5]);
        }else if(cw_f[1][2] && cw_f[1][6] && cw_f[0][5] && cw_f[0][8] &&
                 !contain(v, faces[solid_type][2])){
            cp = getCP_sq(verts[2], verts[0], verts[3],
                        ws[2], ws[6], -ws[5], -ws[8], edges[5], edges[8], raydir);
        }else if(cw_f[1][0] && cw_f[1][7] && cw_f[0][3] && cw_f[0][6] &&
                 !contain(v, faces[solid_type][3])){
            cp = getCP_sq(verts[0], verts[1], verts[4],
                        ws[0], ws[7], -ws[3], -ws[6], edges[3], edges[6], raydir);
        }else if(cw_f[1][1] && cw_f[1][8] && cw_f[0][4] && cw_f[0][7] &&
                 !contain(v, faces[solid_type][4])){
            cp = getCP_sq(verts[1], verts[2], verts[5],
                        ws[1], ws[8], -ws[4], -ws[7], edges[4], edges[7], raydir);
        }else {
          d_out[v] = 0;
          continue;
        }
        break;
      case 8:     //hexahedron
        if (cw_f[0][0] && cw_f[0][1] && cw_f[0][2] && cw_f[0][3] &&
            !contain(v, faces[solid_type][0])){
          cp = getCP_sq(verts[0], verts[1], verts[2], ws[0], ws[1], ws[2], ws[3],
                        edges[2], edges[3], raydir);
        }else if (cw_f[1][4] && cw_f[1][5] && cw_f[1][6] && cw_f[1][7] &&
                  !contain(v, faces[solid_type][1])){
          cp = getCP_sq(verts[4], verts[5], verts[6], ws[4], ws[5], ws[6], ws[7],
                        edges[6], edges[7],raydir);
        }else if (cw_f[1][0] && cw_f[1][9] && cw_f[0][4] && cw_f[0][8] &&
                 !contain(v, faces[solid_type][2])){
          cp = getCP_sq(verts[0], verts[1], verts[5], ws[0], ws[9], -ws[4], -ws[8],
                        edges[4], edges[8],raydir);
        }else if (cw_f[1][1] && cw_f[1][10] && cw_f[0][5] && cw_f[0][9] &&
                 !contain(v, faces[solid_type][3])){
          cp = getCP_sq(verts[1], verts[2], verts[6], ws[1], ws[10], -ws[5], -ws[9],
                        edges[5], edges[9],raydir);
        }else if (cw_f[1][2] && cw_f[1][11] && cw_f[0][6] && cw_f[0][10] &&
                 !contain(v, faces[solid_type][4])){
          cp = getCP_sq(verts[2], verts[3], verts[7], ws[2], ws[11], -ws[6], -ws[10],
                        edges[6], edges[10],raydir);
        }else if (cw_f[1][3] && cw_f[1][8] && cw_f[0][7] && cw_f[0][11] &&
                 !contain(v, faces[solid_type][5])){
          cp = getCP_sq(verts[3], verts[0], verts[4], ws[3], ws[8], -ws[7], -ws[11],
                        edges[7], edges[11],raydir);
        }else {
          d_out[v] = 0;
          continue;
        }
        break;
    }
    d_out[v] = d / (cp - rayorg).length();
    
    if (d_out[v] < 0) d_out[v] = 0;
    
    else if (!(d_out[v] < 1) && !(d_out[v] > 0)){
      d_out[v] = 1;
    }
    
    switch (interpolate_mode) {
      case 0:
        d_out[v] = 1 - d_out[v];
        break;
      case 1:
        d_out[v] = d_out[v]*d_out[v]*(d_out[v]-1.5) + 0.5;
        break;
      case 2:
        d_out[v] = cos(3.141592*d_out[v])+1;
        break;
      default:
        d_out[v] = -pow(cos(0.5*3.141592*(d_out[v])), interpolate_mode);
        break;
    }
    
  }
  
  float ds = 0, dss = 0, d[solid_type], w = 0;
  for (int i = 0; i < solid_type; i++) {
    ds += (pt - verts[i]).length();
  }
  
  for (int i = 0; i < solid_type; i++) {
    d[i] = ds/(pt - verts[i]).length();
    dss += d[i];
  }
  
  for (int i = 0; i < solid_type; i++) {
    d_out[i] *= d[i];
    w += d_out[i];
  }
  
  for (int i = 0; i < solid_type; i++) {
    d_out[i] = d_out[i] / w;
  }
  
}

// commpute ray-square_face cross point
FORCEINLINE void SetCrossPoint_Sq(real3 &point,const double3 &v0,const double3 &v1,
                                  const double3 &v2, const double &ws0, const double &ws1,
                                  const double &ws2, const double &ws3,
                                  double3 &edge2, double3 &edge3,const double3& raydir){
  double w = ws2 + ws3 + vdotd(vcrossd(edge2, edge3),raydir);
  point = toreal3( (v2*ws0 + v0*ws1 + v1*w) / (ws0+ws1+w) );
}
  
FORCEINLINE void SetCrossPoint_Sq(real3 &point,const real3 &v0,const real3 &v1,
                                  const real3 &v2, const float &ws0, const float &ws1,
                                  const float &ws2, const float &ws3, real3 &edge2,
                                  real3 &edge3,const real3& raydir){
  double w = ws2 + ws3 + dot(cross(edge2, edge3),raydir);
  point = (v2*ws0 + v0*ws1 + v1*w) / (ws0+ws1+w);
}
  
//Pyramid
bool IntersectPyramidD(const double3& rayorg, const double3& raydir,
                       const double3* vertices, Intersection *isects)
{
  bool cw_f[2][8];
  double ws[8];
  
  double3 edges[8];
  GetEdges(edges, 5, vertices);
  
  double3 raypc = vcrossd(raydir, rayorg);
  double3 edgepc[8];
  
  for(int i = 0; i < 8; i ++){
    edgepc[i] = vcrossd(edges[i],vertices[i%4]);
    ws[i] = vdotd(raydir, edgepc[i]) + vdotd(edges[i], raypc);
    if(ws[i] >= 0)  cw_f[0][i] = true;
    else cw_f[0][i] = false;
    if(ws[i] <= 0)  cw_f[1][i] = true;
    else cw_f[1][i] = false;
  }
    
  for (int i = 1,n = 0; n < 2; i--,n++)
    if(cw_f[i][0] && cw_f[i][1] && cw_f[i][2] && cw_f[i][3]){
      SetCrossPoint_Sq(isects[n].position, vertices[0], vertices[1], vertices[2],
                       ws[0], ws[1], ws[2], ws[3], edges[2], edges[3], raydir);
      isects[n].normal = toreal3(normalize(vcrossd(edges[1], edges[0])));
    }else if (cw_f[n][0] && cw_f[i][4] && cw_f[n][5]){
      isects[n].position = toreal3(vertices[4]*-ws[0] + vertices[1]*ws[4] + vertices[0]*-ws[5]) / (-ws[0]+ws[4]-ws[5]);
      isects[n].normal = toreal3(normalize(vcrossd(edges[4], edges[0])));
    }else if (cw_f[n][1] && cw_f[i][5] && cw_f[n][6]){
      isects[n].position = toreal3(vertices[4]*-ws[1] + vertices[2]*ws[5] + vertices[1]*-ws[6]) / (-ws[1]+ws[5]-ws[6]);
      isects[n].normal = toreal3(normalize(vcrossd(edges[5], edges[1])));
    }else if (cw_f[n][2] && cw_f[i][6] && cw_f[n][7]){
      isects[n].position = toreal3(vertices[4]*-ws[2] + vertices[3]*ws[6] + vertices[2]*-ws[7]) / (-ws[2]+ws[6]-ws[7]);
      isects[n].normal = toreal3(normalize(vcrossd(edges[6], edges[2])));
    }else if (cw_f[n][3] && cw_f[i][7] && cw_f[n][4]){
      isects[n].position = toreal3(vertices[4]*-ws[3] + vertices[0]*ws[7] + vertices[3]*-ws[4]) / (-ws[3]+ws[7]-ws[4]);
      isects[n].normal = toreal3(normalize(vcrossd(edges[7], edges[3])));
    }else
      return false;
    
  for (int j = 0; j < 3; j++)
    if (raydir[j] > 0.1 || raydir[j] < -0.1){
      isects[0].t = (isects[0].position[j] - rayorg[j]) / raydir[j];
      isects[1].t = (isects[1].position[j] - rayorg[j]) / raydir[j];
    }
  
  isects[1].normal = -1 * isects[1].normal;
  
  return true;
}
  
bool IntersectPyramidF(const real3& rayorg, const real3& raydir,
                       const real3* vertices, Intersection *isects)
{
  bool cw_f[2][8];
  real ws[8];
    
  real3 edges[8];
  GetEdges(edges, 5, vertices);
    
  real3 raypc = cross(raydir, rayorg);
  real3 edgepc[8];
  
  for(int i = 0; i < 8; i ++){
    edgepc[i] = cross(edges[i],vertices[i%4]);
    ws[i] = dot(raydir, edgepc[i]) + dot(edges[i], raypc);
    if(ws[i] >= 0)  cw_f[0][i] = true;
    else cw_f[0][i] = false;
    if(ws[i] <= 0)  cw_f[1][i] = true;
    else cw_f[1][i] = false;
  }
    
  for (int i = 1,n = 0; n < 2; i--,n++)
    if(cw_f[i][0] && cw_f[i][1] && cw_f[i][2] && cw_f[i][3]){
      SetCrossPoint_Sq(isects[n].position, vertices[0], vertices[1], vertices[2], ws[0], ws[1], ws[2], ws[3], edges[2], edges[3], raydir);
      isects[i].normal = cross(edges[1], edges[0]).normalize();
    }else if (cw_f[n][0] && cw_f[i][4] && cw_f[n][5]){
      isects[n].position = (vertices[4]*-ws[0] + vertices[1]*ws[4] + vertices[0]*-ws[5]) / (-ws[0]+ws[4]-ws[5]);
      isects[n].normal = cross(edges[4], edges[0]).normalize();
    }else if (cw_f[n][1] && cw_f[i][5] && cw_f[n][6]){
      isects[n].position = (vertices[4]*-ws[1] + vertices[2]*ws[5] + vertices[1]*-ws[6]) / (-ws[1]+ws[5]-ws[6]);
      isects[n].normal = cross(edges[5], edges[1]).normalize();
    }else if (cw_f[n][2] && cw_f[i][6] && cw_f[n][7]){
      isects[n].position = (vertices[4]*-ws[2] + vertices[3]*ws[6] + vertices[2]*-ws[7]) / (-ws[2]+ws[6]-ws[7]);
        isects[n].normal = cross(edges[6], edges[2]).normalize();
      }else if (cw_f[n][3] && cw_f[i][7] && cw_f[n][4]){
        isects[n].position = (vertices[4]*-ws[3] + vertices[0]*ws[7] + vertices[3]*-ws[4]) / (-ws[3]+ws[7]-ws[4]);
        isects[n].normal = cross(edges[7], edges[3]).normalize();
      }else
        return false;
  
  for (int j = 0; j < 3; j++)
    if (raydir[j] > 0.1 || raydir[j] < -0.1){
      isects[0].t = (isects[0].position[j] - rayorg[j]) / raydir[j];
      isects[1].t = (isects[1].position[j] - rayorg[j]) / raydir[j];
    }
  
  isects[1].normal = -1 * isects[1].normal;
    
  return true;
}
  
// Prism
bool IntersectPrismD(const double3& rayorg, const double3& raydir,
                     const double3* verts, Intersection *isects)
{
  bool cw_f[2][9];
  double ws[9];
    
  double3 edges[9];
  GetEdges(edges, 6, verts);
  
  double3 raypc = vcrossd(raydir, rayorg);
  double3 edgepc[9];
  
  // cross flag of ray vs edges.
  bool cross_rve[9];
  
  for(int i = 0; i < 9; i ++){
    edgepc[i] = vcrossd(edges[i], verts[i%6]);
    ws[i] = vdotd(raydir, edgepc[i]) + vdotd(edges[i], raypc);
    if(ws[i] >= 0)  cw_f[0][i] = true;
    else cw_f[0][i] = false;
    if(ws[i] <= 0)  cw_f[1][i] = true;
    else cw_f[1][i] = false;
    if (cw_f[0][i] && cw_f[1][i]) cross_rve[i] = 1;
  }
  
  for (int i = 0,n = 1; i < 2; i++,n--)
    if(cw_f[n][0] && cw_f[n][1] && cw_f[n][2]){
      isects[i].position = toreal3(verts[2]*ws[0] + verts[0]*ws[1] + verts[1]*ws[2]) / (ws[0]+ws[1]+ws[2]);
      isects[i].normal = toreal3( normalize(vcrossd(edges[0], edges[1])));
    } else if(cw_f[i][3] && cw_f[i][4] && cw_f[i][5]){
      isects[i].position = toreal3(verts[5]*ws[3] + verts[3]*ws[4] + verts[4]*ws[5]) / (ws[3]+ws[4]+ws[5]);
      isects[i].normal = toreal3( normalize(vcrossd(edges[3], edges[4])));
    }else if(cw_f[i][2] && cw_f[i][6] && cw_f[n][5] && cw_f[n][8]){
      SetCrossPoint_Sq(isects[i].position, verts[2], verts[0], verts[3], ws[2], ws[6], -ws[5], -ws[8], edges[5], edges[8], raydir);
      isects[i].normal = toreal3( normalize(vcrossd(edges[2],edges[6])));
    }else if(cw_f[i][0] && cw_f[i][7] && cw_f[n][3] && cw_f[n][6]){
      SetCrossPoint_Sq(isects[i].position, verts[0], verts[1], verts[4], ws[0], ws[7], -ws[3], -ws[6], edges[3], edges[6], raydir);
      isects[i].normal = toreal3( normalize(vcrossd(edges[0], edges[7])));
    }else if(cw_f[i][1] && cw_f[i][8] && cw_f[n][4] && cw_f[n][7]){
      SetCrossPoint_Sq(isects[i].position, verts[1], verts[2], verts[5], ws[1], ws[8], -ws[4], -ws[7], edges[4], edges[7], raydir);
      isects[i].normal = toreal3( normalize(vcrossd(edges[1], edges[8])));
    }else
      return false;
    
  for (int j = 0; j < 3; j++)
    if (raydir[j] > 0.1 || raydir[j] < -0.1){
      isects[0].t = (isects[0].position[j] - rayorg[j]) / raydir[j];
      isects[1].t = (isects[1].position[j] - rayorg[j]) / raydir[j];
    }
    
  isects[1].normal = -1 * isects[1].normal;
  return true;
}
  
bool IntersectPrismF(const real3& rayorg, const real3& raydir,
                     const real3* vertices, Intersection *isects)
{
  bool cw_f[2][9];
  real ws[9];
    
  real3 edges[9];
  GetEdges(edges, 6, vertices);
    
  real3 raypc = cross(raydir, rayorg);
  real3 edgepc[9];
    
  for(int i = 0; i < 9; i ++){
    edgepc[i] = cross(edges[i], vertices[i%6]);
    ws[i] = dot(raydir, edgepc[i]) + dot(edges[i], raypc);
    if(ws[i] >= 0)  cw_f[0][i] = true;
    else cw_f[0][i] = false;
    if(ws[i] <= 0)  cw_f[1][i] = true;
    else cw_f[1][i] = false;
  }
  
  for (int i = 0,n = 1; i < 2; i++,n--)
    if(cw_f[n][0] && cw_f[n][1] && cw_f[n][2]){
      isects[i].position = (vertices[2]*ws[0] + vertices[0]*ws[1] + vertices[1]*ws[2]) / (ws[0]+ws[1]+ws[2]);
      isects[i].normal =  cross(edges[0], edges[1]).normalize();
    } else if(cw_f[i][3] && cw_f[i][4] && cw_f[i][5]){
      isects[i].position = (vertices[5]*ws[3] + vertices[3]*ws[4] + vertices[4]*ws[5]) / (ws[3]+ws[4]+ws[5]);
      isects[i].normal = cross(edges[3], edges[4]).normalize();
    }else if(cw_f[i][2] && cw_f[i][6] && cw_f[n][5] && cw_f[n][8]){
      SetCrossPoint_Sq(isects[i].position, vertices[2], vertices[0], vertices[3], ws[2], ws[6], -ws[5], -ws[8], edges[5], edges[8], raydir);
      isects[i].normal = cross(edges[2],edges[6]).normalize();
    }else if(cw_f[i][0] && cw_f[i][7] && cw_f[n][3] && cw_f[n][6]){
      SetCrossPoint_Sq(isects[i].position, vertices[0], vertices[1], vertices[4], ws[0], ws[7], -ws[3], -ws[6], edges[3], edges[6], raydir);
      isects[i].normal = cross(edges[0], edges[7]).normalize();
    }else if(cw_f[i][1] && cw_f[i][8] && cw_f[n][4] && cw_f[n][7]){
      SetCrossPoint_Sq(isects[i].position, vertices[1], vertices[2], vertices[5], ws[1], ws[8], -ws[4], -ws[7], edges[4], edges[7], raydir);
      isects[i].normal = cross(edges[1], edges[8]).normalize();
    }else
      return false;
    
  for (int j = 0; j < 3; j++)
    if (raydir[j] > 0.1 || raydir[j] < -0.1){
      isects[0].t = (isects[0].position[j] - rayorg[j]) / raydir[j];
      isects[1].t = (isects[1].position[j] - rayorg[j]) / raydir[j];
    }
  
  isects[1].normal = -1 * isects[1].normal;
  return true;
}
  
//Hexa
bool IntersectHexaD(const double3& rayorg, const double3& raydir,
                    const double3* vertices, Intersection *isects)
{
  bool cw_f[2][12];
  double ws[12];
    
  double3 edges[12];
  GetEdges(edges, 8, vertices);
    
  double3 raypc = vcrossd(raydir, rayorg);
  double3 edgepc[12];
    
  for(int i = 0; i < 12; i ++){
    edgepc[i] = vcrossd(edges[i], vertices[i%8]);
    ws[i] = vdotd(raydir, edgepc[i]) + vdotd(edges[i], raypc);
    if(ws[i] >= 0)  cw_f[0][i] = true;
    else cw_f[0][i] = false;
    if(ws[i] <= 0)  cw_f[1][i] = true;
    else cw_f[1][i] = false;
  }
    
  for (int i = 0,n = 1; i < 2; i++,n-- )
    if     (cw_f[n][0] && cw_f[n][1] && cw_f[n][2] && cw_f[n][3]){
      SetCrossPoint_Sq(isects[i].position, vertices[0], vertices[1], vertices[2], ws[0], ws[1], ws[2], ws[3], edges[2], edges[3], raydir);
      isects[i].normal = toreal3(normalize(vcrossd(edges[0], edges[1])));
    }else if(cw_f[i][4] && cw_f[i][5] && cw_f[i][6] && cw_f[i][7]){
      SetCrossPoint_Sq(isects[i].position, vertices[4], vertices[5], vertices[6], ws[4], ws[5], ws[6], ws[7], edges[6], edges[7],raydir);
      isects[i].normal = toreal3(normalize(vcrossd(edges[4], edges[5])));
    }else if(cw_f[i][0] && cw_f[i][9] && cw_f[n][4] && cw_f[n][8]){
      SetCrossPoint_Sq(isects[i].position, vertices[0], vertices[1], vertices[5], ws[0], ws[9], -ws[4], -ws[8], edges[4], edges[8],raydir);
      isects[i].normal = toreal3(normalize(vcrossd(edges[0], edges[9])));
    }else if(cw_f[i][1] && cw_f[i][10] && cw_f[n][5] && cw_f[n][9]){
      SetCrossPoint_Sq(isects[i].position, vertices[1], vertices[2], vertices[6], ws[1], ws[10], -ws[5], -ws[9], edges[5], edges[9],raydir);
      isects[i].normal = toreal3(normalize(vcrossd(edges[1], edges[10])));
    }else if(cw_f[i][2] && cw_f[i][11] && cw_f[n][6] && cw_f[n][10]){
      SetCrossPoint_Sq(isects[i].position, vertices[2], vertices[3], vertices[7], ws[2], ws[11], -ws[6], -ws[10], edges[6], edges[10],raydir);
      isects[i].normal = toreal3(normalize(vcrossd(edges[2], edges[11])));
    }else if(cw_f[i][3] && cw_f[i][8] && cw_f[n][7] && cw_f[n][11]){
      SetCrossPoint_Sq(isects[i].position, vertices[3], vertices[0], vertices[4], ws[3], ws[8], -ws[7], -ws[11], edges[7], edges[11],raydir);
      isects[i].normal = toreal3(normalize(vcrossd(edges[3], edges[8])));
    }else
      return false;
    
  for (int j = 0; j < 3; j++)
    if (raydir[j] > 0.1 || raydir[j] < -0.1){
      isects[0].t = (isects[0].position[j] - rayorg[j]) / raydir[j];
      isects[1].t = (isects[1].position[j] - rayorg[j]) / raydir[j];
    }
  isects[1].normal = -1 * isects[1].normal;
  
  return true;
}
  
//Hexa
bool IntersectHexaF(const real3& rayorg, const real3& raydir,
                    const real3* vertices, Intersection *isects)
{
  bool cw_f[2][12];
  real ws[12];
  
  real3 edges[12];
  GetEdges(edges, 8, vertices);
  
  real3 raypc = cross(raydir, rayorg);
  real3 edgepc[12];
  
  for(int i = 0; i < 12; i ++){
    edgepc[i] = cross(edges[i], vertices[i%8]);
    ws[i] = dot(raydir, edgepc[i]) + dot(edges[i], raypc);
    if(ws[i] >= 0)  cw_f[0][i] = true;
    else cw_f[0][i] = false;
    if(ws[i] <= 0)  cw_f[1][i] = true;
    else cw_f[1][i] = false;
  }
    
  for (int i = 0,n = 1; i < 2; i++,n-- )
    if     (cw_f[n][0] && cw_f[n][1] && cw_f[n][2] && cw_f[n][3]){
      SetCrossPoint_Sq(isects[i].position, vertices[0], vertices[1], vertices[2], ws[0], ws[1], ws[2], ws[3], edges[2], edges[3], raydir);
      isects[i].normal = cross(edges[0], edges[1]).normalize();
    }else if(cw_f[i][4] && cw_f[i][5] && cw_f[i][6] && cw_f[i][7]){
      SetCrossPoint_Sq(isects[i].position, vertices[4], vertices[5], vertices[6], ws[4], ws[5], ws[6], ws[7], edges[6], edges[7],raydir);
      isects[i].normal = cross(edges[4], edges[5]).normalize();
    }else if(cw_f[i][0] && cw_f[i][9] && cw_f[n][4] && cw_f[n][8]){
      SetCrossPoint_Sq(isects[i].position, vertices[0], vertices[1], vertices[5], ws[0], ws[9], -ws[4], -ws[8], edges[4], edges[8],raydir);
      isects[i].normal = cross(edges[0], edges[9]).normalize();
    }else if(cw_f[i][1] && cw_f[i][10] && cw_f[n][5] && cw_f[n][9]){
      SetCrossPoint_Sq(isects[i].position, vertices[1], vertices[2], vertices[6], ws[1], ws[10], -ws[5], -ws[9], edges[5], edges[9],raydir);
      isects[i].normal = cross(edges[1], edges[10]).normalize();
    }else if(cw_f[i][2] && cw_f[i][11] && cw_f[n][6] && cw_f[n][10]){
      SetCrossPoint_Sq(isects[i].position, vertices[2], vertices[3], vertices[7], ws[2], ws[11], -ws[6], -ws[10], edges[6], edges[10],raydir);
      isects[i].normal = cross(edges[2], edges[11]).normalize();
    }else if(cw_f[i][3] && cw_f[i][8] && cw_f[n][7] && cw_f[n][11]){
      SetCrossPoint_Sq(isects[i].position, vertices[3], vertices[0], vertices[4], ws[3], ws[8], -ws[7], -ws[11], edges[7], edges[11],raydir);
      isects[i].normal = cross(edges[3], edges[8]).normalize();
    }else
      return false;
    
  for (int j = 0; j < 3; j++)
    if (raydir[j] > 0.1 || raydir[j] < -0.1){
      isects[0].t = (isects[0].position[j] - rayorg[j]) / raydir[j];
      isects[1].t = (isects[1].position[j] - rayorg[j]) / raydir[j];
    }
  isects[1].normal = -1 * isects[1].normal;
  
  
  return true;
}


//
// SAH functions
//

struct BinBuffer {

  BinBuffer(int size) {
    binSize = size;
    bin.resize(2 * 3 * size);
    clear();
  }

  void clear() { memset(&bin[0], 0, sizeof(size_t) * 2 * 3 * binSize); }

  std::vector<size_t> bin; // (min, max) * xyz * binsize
  int binSize;
};

real CalculateSurfaceArea(const real3 &min, const real3 &max) {
  // PrintVec3("sah.bmin", min);
  // PrintVec3("sah.bmax", max);
  real3 box = max - min;
  // PrintVec3("sah.box", box);
  real S = 2.0 * ((box[0] * box[1]) + (box[1] * box[2]) + (box[2] * box[0]));
  // PrintReal("S", S);
  return S;
}

static inline void GetBoundingBoxOfSolid(real3 &bmin, real3 &bmax,
                                         const Solid *solids,
                                         unsigned int index) {
  int numIndices = solids->numVertsPerSolid;

  // real3 p[16]; // Assume numIndices < 16

  bmin = real3(REAL_MAX, REAL_MAX, REAL_MAX);
  bmax = real3(-REAL_MAX, -REAL_MAX, -REAL_MAX);

  for (int j = 0; j < numIndices; j++) {
    unsigned int f = solids->indices[numIndices * index + j];

    real3 p;
    if (solids->isDoublePrecisionPos) {
      p = real3(solids->dvertices[3 * f + 0], solids->dvertices[3 * f + 1],
                solids->dvertices[3 * f + 2]);
    } else {
      p = real3(solids->vertices[3 * f + 0], solids->vertices[3 * f + 1],
                solids->vertices[3 * f + 2]);
    }

    bmin[0] = (std::min)(bmin[0], p[0]);
    bmin[1] = (std::min)(bmin[1], p[1]);
    bmin[2] = (std::min)(bmin[2], p[2]);

    bmax[0] = (std::max)(bmax[0], p[0]);
    bmax[1] = (std::max)(bmax[1], p[1]);
    bmax[2] = (std::max)(bmax[2], p[2]);
  }
}

static void ContributeBinBuffer(BinBuffer *bins, // [out]
                                const real3 &sceneMin, const real3 &sceneMax,
                                const Solid *solids, unsigned int *indices,
                                unsigned int leftIdx, unsigned int rightIdx) {
  static const real EPS = REAL_EPSILON * 16;

  real binSize = (real)bins->binSize;

  // Calculate extent
  real3 sceneSize, sceneInvSize;
  sceneSize = sceneMax - sceneMin;
  for (int i = 0; i < 3; ++i) {
    assert(sceneSize[i] >= 0.0);

    if (sceneSize[i] > EPS) {
      sceneInvSize[i] = binSize / sceneSize[i];
    } else {
      sceneInvSize[i] = 0.0;
    }
  }

  // Clear bin data
  bins->clear();
  // memset(&bins->bin[0], 0, sizeof(2 * 3 * bins->binSize));

  size_t idxBMin[3];
  size_t idxBMax[3];

  for (size_t i = leftIdx; i < rightIdx; i++) {

    //
    // Quantize the position into [0, BIN_SIZE)
    //
    // q[i] = (int)(p[i] - scene_bmin) / scene_size
    //
    real3 bmin;
    real3 bmax;

    GetBoundingBoxOfSolid(bmin, bmax, solids, indices[i]);
    debug("bbox[%d] = %f, %f, %f - %f, %f, %f\n", (int)i, bmin[0], bmin[1],
          bmin[2], bmax[0], bmax[1], bmax[2]);

    real3 quantizedBMin = (bmin - sceneMin) * sceneInvSize;
    real3 quantizedBMax = (bmax - sceneMin) * sceneInvSize;

    // idx is now in [0, BIN_SIZE)
    for (size_t j = 0; j < 3; ++j) {
      idxBMin[j] = (unsigned int)floor(quantizedBMin[j]);
      idxBMax[j] = (unsigned int)floor(quantizedBMax[j]);

      if (idxBMin[j] >= binSize)
        idxBMin[j] = binSize - 1;
      if (idxBMax[j] >= binSize)
        idxBMax[j] = binSize - 1;

      assert(idxBMin[j] < binSize);
      assert(idxBMax[j] < binSize);

      // Increment bin counter
      bins->bin[0 * (bins->binSize * 3) + j * bins->binSize + idxBMin[j]] += 1;
      bins->bin[1 * (bins->binSize * 3) + j * bins->binSize + idxBMax[j]] += 1;
    }
  }

  for (int i = 0; i < 2 * 3 * bins->binSize; i++) {
    if (bins->bin[i] > 0) {
      debug("bins[%d] = %d\n", i, (int)bins->bin[i]);
    }
  }
}

real SAH(size_t ns1, real leftArea, size_t ns2, real rightArea, real invS,
         real Taabb, real Ttri) {
  // const real Taabb = 0.2f;
  // const real Ttri = 0.8f;
  real T;

  T = 2.0f * Taabb + (leftArea * invS) * (real)(ns1)*Ttri +
      (rightArea * invS) * (real)(ns2)*Ttri;

  return T;
}

static bool FindCutFromBinBuffer(real *cutPos,     // [out] xyz
                                 int &minCostAxis, // [out]
                                 const BinBuffer *bins, const real3 &bmin,
                                 const real3 &bmax, size_t numPrims,
                                 real costTaabb) // should be in [0.0, 1.0]
{
  const real eps = REAL_EPSILON * 16;

  size_t left, right;
  real3 bsize, bstep;
  real3 bminLeft, bmaxLeft;
  real3 bminRight, bmaxRight;
  real saLeft, saRight, saTotal;
  real pos;
  real minCost[3];

  real costTtri = 1.0 - costTaabb;

  minCostAxis = 0;

  bsize = bmax - bmin;
  bstep = bsize * (1.0 / bins->binSize);
  saTotal = CalculateSurfaceArea(bmin, bmax);

  real invSaTotal = 0.0;
  if (saTotal > eps) {
    invSaTotal = 1.0 / saTotal;
  }

  for (int j = 0; j < 3; ++j) {

    //
    // Compute SAH cost for right side of each cell of the bbox.
    // Exclude both extreme side of the bbox.
    //
    //  i:      0    1    2    3
    //     +----+----+----+----+----+
    //     |    |    |    |    |    |
    //     +----+----+----+----+----+
    //

    real minCostPos = bmin[j] + 0.5f * bstep[j];
    minCost[j] = REAL_MAX;

    left = 0;
    right = numPrims;
    bminLeft = bminRight = bmin;
    bmaxLeft = bmaxRight = bmax;

    for (int i = 0; i < bins->binSize - 1; ++i) {
      left += bins->bin[0 * (3 * bins->binSize) + j * bins->binSize + i];
      right -= bins->bin[1 * (3 * bins->binSize) + j * bins->binSize + i];

      assert(left <= numPrims);
      assert(right <= numPrims);

      //
      // Split pos bmin + (i + 1) * (bsize / BIN_SIZE)
      // +1 for i since we want a position on right side of the cell.
      //

      pos = bmin[j] + (i + 0.5f) * bstep[j];
      bmaxLeft[j] = pos;
      bminRight[j] = pos;

      saLeft = CalculateSurfaceArea(bminLeft, bmaxLeft);
      saRight = CalculateSurfaceArea(bminRight, bmaxRight);

      real cost =
          SAH(left, saLeft, right, saRight, invSaTotal, costTaabb, costTtri);
      debug("[%d] cost = %f(0x%08x), left = %d, right = %d\n", j, cost,
            *((unsigned int *)&cost), left, right);
      if (cost < minCost[j]) {
        debug("saLeft = %f(0x%08x), saRight = %f(0x%08x0\n", saLeft,
              *((unsigned int *)&saLeft), saRight, *((unsigned int *)&saRight));

        debug("[%d] i = %d, MinCost = %f(0x%08x), cutPos = %f(0x%08x)\n", j, i,
              cost, *((unsigned int *)&cost), pos, *((unsigned int *)&cost));
        //
        // Update the min cost
        //
        minCost[j] = cost;
        minCostPos = pos;
        // minCostAxis = j;
      }
    }

    cutPos[j] = minCostPos;
  }

  // cutAxis = minCostAxis;
  // cutPos = minCostPos;

  // Find min cost axis
  real cost = minCost[0];
  minCostAxis = 0;
  if (cost > minCost[1]) {
    minCostAxis = 1;
    cost = minCost[1];
  }
  if (cost > minCost[2]) {
    minCostAxis = 2;
    cost = minCost[2];
  }

  return true;
}

class SAHPred : public std::unary_function<unsigned int, bool> {
public:
  SAHPred(int axis, real pos, const Solid *solids)
      : axis_(axis), pos_(pos), solids_(solids) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;
    real pos = pos_;

    real center = 0.0;

    for (int j = 0; j < solids_->numVertsPerSolid; j++) {
      unsigned int f = solids_->indices[solids_->numVertsPerSolid * i + j];

      center += solids_->vertices[3 * f + axis];
    }

    return (center < pos * (real)(solids_->numVertsPerSolid));
  }

private:
  int axis_;
  real pos_;
  const Solid *solids_;
};

class SAHPredD : public std::unary_function<unsigned int, bool> {
public:
  SAHPredD(int axis, double pos, const Solid *solids)
      : axis_(axis), pos_(pos), solids_(solids) {}

  bool operator()(unsigned int i) const {
    int axis = axis_;
    double pos = pos_;

    real center = 0.0;

    for (int j = 0; j < solids_->numVertsPerSolid; j++) {
      unsigned int f = solids_->indices[solids_->numVertsPerSolid * i + j];

      center += solids_->dvertices[3 * f + axis];
    }

    return (center < pos * (real)(solids_->numVertsPerSolid));
  }

private:
  int axis_;
  double pos_;
  const Solid *solids_;
};

template <typename T>
static void ComputeBoundingBox(real3 &bmin, real3 &bmax, int numVertsPerSolid,
                               const T *vertices, unsigned int *indices,
                               unsigned int *primIds, unsigned int leftIndex,
                               unsigned int rightIndex) {
  const real kEPS = REAL_EPSILON * 16;

  bmin[0] = REAL_MAX;
  bmin[1] = REAL_MAX;
  bmin[2] = REAL_MAX;
  bmax[0] = -REAL_MAX;
  bmax[1] = -REAL_MAX;
  bmax[2] = -REAL_MAX;

  if (rightIndex <= leftIndex) {
    return;
  }

  size_t i = leftIndex;
  size_t idx = primIds[i];
  bmin[0] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 2] + kEPS;

  for (i = leftIndex; i < rightIndex; i++) { // for each primitives
    size_t idx = primIds[i];
    debug("idx = %d\n", (int)idx);
    for (int j = 0; j < numVertsPerSolid; j++) {
      size_t fid = indices[numVertsPerSolid * idx + j];
      for (int k = 0; k < 3; k++) { // xyz
        real minval = vertices[3 * fid + k] - kEPS;
        real maxval = vertices[3 * fid + k] + kEPS;
        if (bmin[k] > minval)
          bmin[k] = minval;
        if (bmax[k] < maxval)
          bmax[k] = maxval;
      }
    }
  }
}

template <typename T>
static void ComputeBoundingBox30(real3 &bmin, real3 &bmax, int numVertsPerSolid,
                                 const T *vertices, const uint32_t *indices,
                                 const std::vector<IndexKey30> &keys,
                                 unsigned int leftIndex,
                                 unsigned int rightIndex) {
  const real kEPS = REAL_EPSILON * 16;

  bmin[0] = REAL_MAX;
  bmin[1] = REAL_MAX;
  bmin[2] = REAL_MAX;
  bmax[0] = -REAL_MAX;
  bmax[1] = -REAL_MAX;
  bmax[2] = -REAL_MAX;

  if ((rightIndex - leftIndex) == 0) {
    // empty.
    return;
  }

  if (leftIndex >= rightIndex) {
    printf("left = %d, right = %d\n", leftIndex, rightIndex);
    assert(leftIndex < rightIndex);
  }

  size_t i = leftIndex;
  size_t idx = keys[i].index;
  bmin[0] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 0] - kEPS;
  bmin[1] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 1] - kEPS;
  bmin[2] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 2] - kEPS;
  bmax[0] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 0] + kEPS;
  bmax[1] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 1] + kEPS;
  bmax[2] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 2] + kEPS;

  for (i = leftIndex; i < rightIndex; i++) { // for each primitives
    size_t idx = keys[i].index;
    for (int j = 0; j < numVertsPerSolid; j++) {
      size_t fid = indices[numVertsPerSolid * idx + j];
      for (int k = 0; k < 3; k++) { // xyz
        real minval = vertices[3 * fid + k] - kEPS;
        real maxval = vertices[3 * fid + k] + kEPS;
        if (bmin[k] > minval)
          bmin[k] = minval;
        if (bmax[k] < maxval)
          bmax[k] = maxval;
      }
    }
  }
}

#ifdef _OPENMP
template <typename T>
void ComputeBoundingBoxOMP(real3 &bmin, real3 &bmax, int numVertsPerSolid,
                           const T *vertices, const uint32_t *indices,
                           unsigned int *primIndices, unsigned int leftIndex,
                           unsigned int rightIndex) {

  // assert(leftIndex < rightIndex);
  // assert(rightIndex - leftIndex > 0);

  const real kEPS = REAL_EPSILON * 16;

  bmin[0] = (std::numeric_limits<real>::max)();
  bmin[1] = (std::numeric_limits<real>::max)();
  bmin[2] = (std::numeric_limits<real>::max)();
  bmax[0] = -(std::numeric_limits<real>::max)();
  bmax[1] = -(std::numeric_limits<real>::max)();
  bmax[2] = -(std::numeric_limits<real>::max)();

  {
    size_t i = leftIndex;
    size_t idx = primIndices[i];
    bmin[0] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 0] - kEPS;
    bmin[1] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 1] - kEPS;
    bmin[2] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 2] - kEPS;
    bmax[0] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 0] + kEPS;
    bmax[1] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 1] + kEPS;
    bmax[2] = vertices[3 * indices[numVertsPerSolid * idx + 0] + 2] + kEPS;
  }

  // We could use min and max reduction if OpenMP 3.1 ready compiler was
  // available.

  real local_bmin[3] = {bmin[0], bmin[1], bmin[2]};
  real local_bmax[3] = {bmax[0], bmax[1], bmax[2]};

  size_t n = rightIndex - leftIndex;

#pragma omp parallel firstprivate(local_bmin, local_bmax) if (n > (1024 * 124))
  {

#pragma omp for
    for (long long i = leftIndex; i < rightIndex; i++) {

      size_t idx = primIndices[i];

      for (int k = 0; k < numVertsPerSolid; k++) {
        real minval_x =
            vertices[3 * indices[numVertsPerSolid * idx + k] + 0] - kEPS;
        real minval_y =
            vertices[3 * indices[numVertsPerSolid * idx + k] + 1] - kEPS;
        real minval_z =
            vertices[3 * indices[numVertsPerSolid * idx + k] + 2] - kEPS;

        real maxval_x =
            vertices[3 * indices[numVertsPerSolid * idx + k] + 0] + kEPS;
        real maxval_y =
            vertices[3 * indices[numVertsPerSolid * idx + k] + 1] + kEPS;
        real maxval_z =
            vertices[3 * indices[numVertsPerSolid * idx + k] + 2] + kEPS;

        local_bmin[0] = (std::min)(local_bmin[0], minval_x);
        local_bmin[1] = (std::min)(local_bmin[1], minval_y);
        local_bmin[2] = (std::min)(local_bmin[2], minval_z);

        local_bmax[0] = (std::max)(local_bmax[0], maxval_x);
        local_bmax[1] = (std::max)(local_bmax[1], maxval_y);
        local_bmax[2] = (std::max)(local_bmax[2], maxval_z);
      }
    }

#pragma omp critical
    {
      for (int k = 0; k < 3; k++) { // xyz

        if (local_bmin[k] < bmin[k]) {
          {
            if (local_bmin[k] < bmin[k])
              bmin[k] = local_bmin[k];
          }
        }

        if (local_bmax[k] > bmax[k]) {
          {
            if (local_bmax[k] > bmax[k])
              bmax[k] = local_bmax[k];
          }
        }
      }
    }
  }
}
#endif

inline void InvalidateBoundingBox(real3 &bmin, real3 &bmax) {
  bmin[0] = bmin[1] = bmin[2] = (std::numeric_limits<real>::max)();
  bmax[0] = bmax[1] = bmax[2] = -(std::numeric_limits<real>::max)();
}

inline void MergeBoundingBox(real3 &bmin, real3 &bmax, const real3 &leftBMin,
                             const real3 &leftBMax, const real3 &rightBMin,
                             const real3 &rightBMax) {
  bmin = leftBMin;
  bmax = leftBMax;

  for (int k = 0; k < 3; k++) {
    bmin[k] = (std::min)(bmin[k], rightBMin[k]);
    bmax[k] = (std::max)(bmax[k], rightBMax[k]);
  }
}

inline int CountLeadingZeros32(uint32_t x) {
  // Note: we can use clz(), count leading zeros instruction, or
  //       convert integert to float then caculate 31 - floor(log2(x))
  //       for better performance
  //
  // If you have popc instruction, clz() can be computed as,
  //
  // function clz(x):
  //       for each y in {1, 2, 4, 8, 16}: x ← x | (x >> y)
  //       return 32 − popc(x)
  //

  if (x == 0)
    return 32;

  int n = 0;
  if ((x & 0xFFFF0000) == 0) {
    n = n + 16;
    x = x << 16;
  }
  if ((x & 0xFF000000) == 0) {
    n = n + 8;
    x = x << 8;
  }
  if ((x & 0xF0000000) == 0) {
    n = n + 4;
    x = x << 4;
  }
  if ((x & 0xC0000000) == 0) {
    n = n + 2;
    x = x << 2;
  }
  if ((x & 0x80000000) == 0) {
    n = n + 1;
  }
  return n;
}

int GetSplitAxis(uint32_t key) {
  int clz = CountLeadingZeros32(key);

  int n = clz - 2;
  // assert(n >= 0);
  // assert(n < 30);

  // 0 -> x, 1 -> y, 2 -> z, 3 -> x, ...
  return n % 3;
}

template <typename T>
void MakeLeaf32(SolidNode &leaf, int numVersPerSolid, const T *vertices,
                const uint32_t *indices, real3 &bmin, real3 &bmax,
                const std::vector<IndexKey30> &keys, uint32_t leftIndex,
                uint32_t rightIndex) {

  // 1. Compute leaf AABB
  ComputeBoundingBox30(bmin, bmax, numVersPerSolid, vertices, indices, keys,
                       leftIndex, rightIndex);

  // 2. Create leaf node. `n' may be null(create empty leaf node for that case.)
  int n = rightIndex - leftIndex;

  leaf.bmin[0] = bmin[0];
  leaf.bmin[1] = bmin[1];
  leaf.bmin[2] = bmin[2];

  leaf.bmax[0] = bmax[0];
  leaf.bmax[1] = bmax[1];
  leaf.bmax[2] = bmax[2];

  leaf.flag = 1; // leaf
  leaf.data[0] = n;
  leaf.data[1] = (uint32_t)leftIndex;

  return;
}

//
// Build BVH tree from bottom to up manner
//
template <typename T>
size_t BuildTreeRecursive32(std::vector<SolidNode> &nodes, real3 &bmin,
                            real3 &bmax, int numVersPerSolid,
                            const std::vector<IndexKey30> &keys,
                            const std::vector<NodeInfo32> &nodeInfos,
                            const T *vertices, const uint32_t *indices,
                            uint32_t rootIndex, uint32_t leftIndex,
                            uint32_t rightIndex, bool isLeaf, int depth) {
  InvalidateBoundingBox(bmin, bmax);

  uint32_t n = rightIndex - leftIndex;

  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= MAX_LEAF_ELEMENTS) || (depth > MAX_TREE_DEPTH_32BIT)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;
    // if (leftIndex == rightIndex) { // this would be OK. 1 tri in 1 leaf case.
    //  endIndex++;
    //}

    SolidNode leaf;
    MakeLeaf32(leaf, numVersPerSolid, vertices, indices, bmin, bmax, keys,
               leftIndex, endIndex);

    size_t offset = nodes.size();
    nodes.push_back(leaf); // need atomic update.

    return offset;
  }

  //
  // Intermediate node. We already know split position.
  //
  uint32_t midIndex = nodeInfos[rootIndex].childIndex;

  bool isLeftLeaf =
      (nodeInfos[rootIndex].leftType == NODE_TYPE_LEAF) ? true : false;
  bool isRightLeaf =
      (nodeInfos[rootIndex].rightType == NODE_TYPE_LEAF) ? true : false;

  SolidNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  size_t offset = nodes.size();
  nodes.push_back(node);

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  if (midIndex < leftIndex) {
    printf("rootIndex = %d, mid = %d, left = %d, leftLeaf = %d\n", rootIndex,
           midIndex, leftIndex, isLeftLeaf);
    assert(leftIndex <= midIndex);
  }

  if (midIndex > rightIndex) {
    printf("rootIndex = %d, mid = %d, right = %d, rightLeaf = %d\n", rootIndex,
           midIndex, rightIndex, isRightLeaf);
    assert(midIndex <= rightIndex);
  }

  size_t leftChildIndex = BuildTreeRecursive32(
      nodes, leftBMin, leftBMax, numVersPerSolid, keys, nodeInfos, vertices,
      indices, midIndex, leftIndex, midIndex, isLeftLeaf, depth + 1);
  size_t rightChildIndex = BuildTreeRecursive32(
      nodes, rightBMin, rightBMax, numVersPerSolid, keys, nodeInfos, vertices,
      indices, midIndex + 1, midIndex + 1, rightIndex, isRightLeaf, depth + 1);

  MergeBoundingBox(bmin, bmax, leftBMin, leftBMax, rightBMin, rightBMax);

  // printf("[%d] -> (%d, %d)\n", offset, leftChildIndex, rightChildIndex);

  node.data[0] = leftChildIndex;
  node.data[1] = rightChildIndex;

  node.bmin[0] = bmin[0];
  node.bmin[1] = bmin[1];
  node.bmin[2] = bmin[2];

  node.bmax[0] = bmax[0];
  node.bmax[1] = bmax[1];
  node.bmax[2] = bmax[2];

  nodes[offset] = node;

  return offset;
}

#if defined(_OPENMP) && !defined(_MSC_VER)
template <typename T>
size_t BuildTreeRecursive32OMP(std::vector<SolidNode> &nodes, real3 &bmin,
                               real3 &bmax, int numVersPerSolid,
                               const std::vector<IndexKey30> &keys,
                               const std::vector<NodeInfo32> &nodeInfos,
                               const T *vertices, const uint32_t *indices,
                               uint32_t rootIndex, uint32_t leftIndex,
                               uint32_t rightIndex, bool isLeaf, int depth) {
  InvalidateBoundingBox(bmin, bmax);

  uint32_t n = rightIndex - leftIndex;

  // printf("[%d] rootIndex = %d, range (%d - %d), leaf = %d, n = %d\n", depth,
  // rootIndex, leftIndex, rightIndex, isLeaf, n);

  if (isLeaf || (n <= MAX_LEAF_ELEMENTS) || (depth > MAX_TREE_DEPTH_32BIT)) {
    // if (isLeaf || (n <= 0)) {
    // printf("  makeLeaf: range (%d - %d)\n", leftIndex, rightIndex);

    uint32_t endIndex = rightIndex + 1;

    SolidNode leaf;
    MakeLeaf32(leaf, numVersPerSolid, vertices, indices, bmin, bmax, keys,
               leftIndex, endIndex);

    // @{ critical section }
    size_t offset;
#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(leaf); // need atomic update.
    }

    return offset;
  }

  //
  // Intermediate node. We already know split position.
  //
  uint32_t midIndex = nodeInfos[rootIndex].childIndex;
  bool isLeftLeaf =
      (nodeInfos[rootIndex].leftType == NODE_TYPE_LEAF) ? true : false;
  bool isRightLeaf =
      (nodeInfos[rootIndex].rightType == NODE_TYPE_LEAF) ? true : false;

  SolidNode node;
  node.axis = GetSplitAxis(keys[rootIndex].code);
  node.flag = 0; // 0 = branch

  size_t offset;

  // Recurively split tree.
  real3 leftBMin, leftBMax;
  real3 rightBMin, rightBMax;

  size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

  if (depth > 6) {

    // Enough number of tasks was launched. Switch to sequential code.
    std::vector<SolidNode> sub_nodes;

    leftChildIndex =
        BuildTreeRecursive32(sub_nodes, leftBMin, leftBMax, numVersPerSolid,
                             keys, nodeInfos, vertices, indices, midIndex,
                             leftIndex, midIndex, isLeftLeaf, depth + 1);

    rightChildIndex =
        BuildTreeRecursive32(sub_nodes, rightBMin, rightBMax, numVersPerSolid,
                             keys, nodeInfos, vertices, indices, midIndex + 1,
                             midIndex + 1, rightIndex, isRightLeaf, depth + 1);

#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(node);

      // printf("offset = %d\n", offset);

      // add sub nodes
      for (size_t i = 0; i < sub_nodes.size(); i++) {
        if (sub_nodes[i].flag == 0) { // branch node
          // printf("sub[%d] = %d, %d\n", i, sub_nodes[i].data[0],
          // sub_nodes[i].data[1]);
          sub_nodes[i].data[0] += offset + 1;
          sub_nodes[i].data[1] += offset + 1;
          assert(sub_nodes[i].data[0] < nodes.size() + sub_nodes.size());
          assert(sub_nodes[i].data[1] < nodes.size() + sub_nodes.size());
          // printf("-> offt = %d, sublen = %d, sub[%d] = %d, %d\n", offset,
          // sub_nodes.size(), i, sub_nodes[i].data[0], sub_nodes[i].data[1]);
        }
      }

      nodes.insert(nodes.end(), sub_nodes.begin(), sub_nodes.end());
      leftChildIndex += offset + 1;
      rightChildIndex += offset + 1;
      assert(leftChildIndex < nodes.size() + sub_nodes.size());
      assert(rightChildIndex < nodes.size() + sub_nodes.size());
    }

  } else {

#pragma omp critical
    {
      offset = nodes.size();
      nodes.push_back(node);
    }

#pragma omp task shared(leftChildIndex, rightChildIndex, nodes, leftBMin,      \
                        leftBMax, keys, nodeInfos, vertices,                   \
                        indices) firstprivate(midIndex) if (depth < 10)
    leftChildIndex = BuildTreeRecursive32OMP(
        nodes, leftBMin, leftBMax, numVersPerSolid, keys, nodeInfos, vertices,
        indices, midIndex, leftIndex, midIndex, isLeftLeaf, depth + 1);

#pragma omp task shared(leftIndex, rightChildIndex, nodes, rightBMin,          \
                        rightBMax, keys, nodeInfos, vertices,                  \
                        indices) firstprivate(midIndex) if (depth < 10)
    rightChildIndex = BuildTreeRecursive32OMP(
        nodes, rightBMin, rightBMax, numVersPerSolid, keys, nodeInfos, vertices,
        indices, midIndex + 1, midIndex + 1, rightIndex, isRightLeaf,
        depth + 1);

#pragma omp taskwait
  }

  assert(leftChildIndex != (size_t)(-1));
  assert(rightChildIndex != (size_t)(-1));

  MergeBoundingBox(bmin, bmax, leftBMin, leftBMax, rightBMin, rightBMax);

  // printf("[%d] -> (%d, %d)\n", offset, leftChildIndex, rightChildIndex);

  node.data[0] = leftChildIndex;
  node.data[1] = rightChildIndex;

  node.bmin[0] = bmin[0];
  node.bmin[1] = bmin[1];
  node.bmin[2] = bmin[2];

  node.bmax[0] = bmax[0];
  node.bmax[1] = bmax[1];
  node.bmax[2] = bmax[2];

#pragma omp critical
  { nodes[offset] = node; }

  return offset;
}
#endif

} // namespace

//
// --
//

size_t SolidAccel::BuildTree(const Solid *solids, unsigned int leftIdx,
                             unsigned int rightIdx, int depth) {
  assert(leftIdx <= rightIdx);

  debug("d: %d, l: %d, r: %d\n", depth, leftIdx, rightIdx);

  assert(solids);
  assert(solids->indices);

  size_t offset = nodes_.size();

  if (buildStats_.maxTreeDepth < depth) {
    buildStats_.maxTreeDepth = depth;
  }

  real3 bmin, bmax;
  if (solids->isDoublePrecisionPos) {
    ComputeBoundingBox<double>(bmin, bmax, solids->numVertsPerSolid,
                               solids->dvertices, solids->indices,
                               &indices_.at(0), leftIdx, rightIdx);
  } else {
    ComputeBoundingBox<float>(bmin, bmax, solids->numVertsPerSolid,
                              solids->vertices, solids->indices,
                              &indices_.at(0), leftIdx, rightIdx);
  }

  size_t n = rightIdx - leftIdx;
  if ((n < options_.minLeafPrimitives) || (depth >= options_.maxTreeDepth)) {
    // Create leaf node.
    SolidNode leaf;

    debug("leaf.bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    debug("leaf.bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

    leaf.bmin[0] = bmin[0];
    leaf.bmin[1] = bmin[1];
    leaf.bmin[2] = bmin[2];

    leaf.bmax[0] = bmax[0];
    leaf.bmax[1] = bmax[1];
    leaf.bmax[2] = bmax[2];

    assert(leftIdx < (std::numeric_limits<unsigned int>::max)());

    leaf.flag = 1; // leaf
    leaf.data[0] = n;
    leaf.data[1] = (unsigned int)leftIdx;

    nodes_.push_back(leaf);

    buildStats_.numLeafNodes++;

    return offset;
  }

  //
  // Create branch node.
  //

  //
  // Compute SAH and find best split axis and position
  //
  int minCutAxis = 0;
  real cutPos[3] = {0.0, 0.0, 0.0};

  BinBuffer bins(options_.binSize);
  debug("binSize = %d\n", options_.binSize);
  ContributeBinBuffer(&bins, bmin, bmax, solids, &indices_.at(0), leftIdx,
                      rightIdx);
  FindCutFromBinBuffer(cutPos, minCutAxis, &bins, bmin, bmax, n,
                       options_.costTaabb);

  debug("depth: %d, cutPos: (%f, %f, %f), cutAxis: %d\n", depth, cutPos[0],
        cutPos[1], cutPos[2], minCutAxis);

  // Try all 3 axis until good cut position avaiable.
  unsigned int midIdx;
  int cutAxis = minCutAxis;
  for (int axisTry = 0; axisTry < 3; axisTry++) {

    unsigned int *begin = &indices_[leftIdx];
    unsigned int *end = &indices_[rightIdx - 1] + 1;
    unsigned int *mid = 0;

    // try minCutAxis first.
    cutAxis = (minCutAxis + axisTry) % 3;

    //
    // Split at (cutAxis, cutPos)
    // indices_ will be modified.
    //
    if (solids->isDoublePrecisionPos) {
      mid = std::partition(begin, end,
                           SAHPredD(cutAxis, cutPos[cutAxis], solids));
    } else {
      mid =
          std::partition(begin, end, SAHPred(cutAxis, cutPos[cutAxis], solids));
    }

    midIdx = leftIdx + (mid - begin);
    if ((midIdx == leftIdx) || (midIdx == rightIdx)) {

      // Can't split well.
      // Switch to object median(which may create unoptimized tree, but stable)
      midIdx = leftIdx + (n >> 1);

      // Try another axis if there's axis to try.

    } else {

      // Found good cut. exit loop.
      break;
    }
  }

  debug("cutAxis = %d, midIdx = %d\n", cutAxis, midIdx);

  SolidNode node;
  node.axis = cutAxis;
  node.flag = 0; // 0 = branch
  nodes_.push_back(node);

  // Recurively split tree.
  unsigned int leftChildIndex = BuildTree(solids, leftIdx, midIdx, depth + 1);
  unsigned int rightChildIndex = BuildTree(solids, midIdx, rightIdx, depth + 1);

  nodes_[offset].data[0] = leftChildIndex;
  nodes_[offset].data[1] = rightChildIndex;

  nodes_[offset].bmin[0] = bmin[0];
  nodes_[offset].bmin[1] = bmin[1];
  nodes_[offset].bmin[2] = bmin[2];

  nodes_[offset].bmax[0] = bmax[0];
  nodes_[offset].bmax[1] = bmax[1];
  nodes_[offset].bmax[2] = bmax[2];

  buildStats_.numBranchNodes++;

  return offset;
}

bool SolidAccel::Build(const Solid *solids, const SolidBuildOptions &options) {
  options_ = options;
  buildStats_ = SolidBuildStatistics();

  assert(options_.binSize > 1);

  assert(solids);

  size_t n = solids->numSolids;
  trace("[SolidAccel] Input # of solids    = %lu\n", solids->numSolids);

  //
  // 1. Create primtive indices(this will be permutated in BuildTree)
  //
  indices_.resize(n);
  for (size_t i = 0; i < n; i++) {
    indices_[i] = i;
  }

  //
  // 2. Build tree
  //
  BuildTree(solids, 0, n, 0);

  // Tree will be null if input triangle count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    real3 bmin(&nodes_[0].bmin[0]);
    real3 bmax(&nodes_[0].bmax[0]);
    trace("[SolidAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
    trace("[SolidAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);
  }

  trace("[SolidAccel] # of nodes = %lu\n", nodes_.size());

  // Store pointer for later use.
  solids_ = solids;

  return true;
}

bool SolidAccel::Build32(const Solid *solids,
                         const SolidBuildOptions &options) {

  options_ = options;
  buildStats_ = SolidBuildStatistics();

  assert(options_.binSize > 1);

  assert(solids);

  size_t n = solids->numSolids;

  if (n < 1024) {
    // Use non-optimized BVH builder.
    return Build(solids, options);
  }

  trace("[LSGL] [SolidAccel2] Input # of solids    = %lu\n", solids->numSolids);

  //
  // 1. Create indices(this will be permutated in BuildTree)
  //
  timerutil t;
  t.start();
  indices_.resize(n);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (long long i = 0; i < n; i++) {
    indices_[i] = i;
  }
  t.end();
  trace("[LSGL] [1:indexing] %d msec\n", (int)t.msec());

  {
    timerutil t;
    t.start();

    real3 bmin, bmax;

    if (solids->isDoublePrecisionPos) {

#ifdef _OPENMP
      ComputeBoundingBoxOMP(bmin, bmax, solids->numVertsPerSolid,
                            solids->dvertices, solids->indices, &indices_.at(0),
                            0, n);
#else
      ComputeBoundingBox(bmin, bmax, solids->numVertsPerSolid,
                         solids->dvertices, solids->indices, &indices_.at(0), 0,
                         n);
#endif

    } else {

#ifdef _OPENMP
      ComputeBoundingBoxOMP(bmin, bmax, solids->numVertsPerSolid,
                            solids->vertices, solids->indices, &indices_.at(0),
                            0, n);
#else
      ComputeBoundingBox(bmin, bmax, solids->numVertsPerSolid, solids->vertices,
                         solids->indices, &indices_.at(0), 0, n);
#endif
    }

    t.end();
    trace("[LSGL] [2:scene bbox calculation] %d msec\n", (int)t.msec());

    trace("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
    trace("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

    std::vector<uint32_t> codes(n);
    assert(sizeof(real) == sizeof(float));

    {
      timerutil t;
      t.start();

      if (solids->isDoublePrecisionPos) {
        CalculateMortonCodesSolidDouble30(
            &codes.at(0), solids->numVertsPerSolid, solids->dvertices,
            solids->indices, bmin, bmax, 0, n);
      } else {
        CalculateMortonCodesSolidFloat30(&codes.at(0), solids->numVertsPerSolid,
                                         solids->vertices, solids->indices,
                                         bmin, bmax, 0, n);
      }
      t.end();
      trace("[LSGL] [3:morton calculation] %d msec\n", (int)t.msec());

      // for (size_t i = 0; i < n; i++) {
      //  printf("code[%d] = %d\n", i, codes[i]);
      //}
    }

    std::vector<IndexKey30> keys(n);

#ifdef _OPENMP
#pragma omp parallel for if (n > (1024 * 1024))
#endif
    for (long long i = 0; i < n; i++) {
      keys[i].index = indices_[i];
      keys[i].code = codes[i];
    }

    { // sort
      timerutil t;
      t.start();
      std::vector<IndexKey30> temp(n);

#ifdef _OPENMP
#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10
      // OpenMP version doesn't work well. Use sequential method for now.
      RadixSort30(&keys.at(0), &keys.at(0) + n);
#else
      RadixSort30OMP(&keys.at(0), &keys.at(0) + n);
#pragma omp barrier
#endif
      t.end();
      trace("[LSGL] [4:radix sort] %d msec\n", (int)t.msec());

#else
      RadixSort30(&keys.at(0), &keys.at(0) + n);
#endif

      // Check
      // for (size_t i = 0; i < n-1; i++) {
      // std::string b = BitString32(keys[i].code);
      // printf("[%08d] i = %010d, c = %010d(%s)\n", i, keys[i].index,
      // keys[i].code, b.c_str());
      // assert(keys[i].code <= keys[i+1].code);
      //}
    }

    std::vector<NodeInfo32> nodeInfos(n - 1);
    {
      timerutil t;
      t.start();

#ifdef _OPENMP
//#pragma omp parallel forif (n > (1024 * 1024))
#pragma omp parallel for schedule(dynamic)
#endif
      for (long long i = 0; i < n - 1; i++) {
        nodeInfos[i] = ConstructBinaryRadixTree30(&keys.at(0), i, n);
        // printf("I[%d].index   = %d\n", i, nodeInfos[i].index);
        // printf("I[%d].leftTy  = %d\n", i, nodeInfos[i].leftType);
        // printf("I[%d].rightTy = %d\n", i, nodeInfos[i].rightType);
      }

      t.end();
      trace("[LSGL] [5:Construct binary radix tree: %d msec\n", (int)t.msec());
    }

    {
      timerutil t;
      t.start();

      nodes_.clear();

      // Explicitly create root node and reserve storage here.
      SolidNode rootNode;
      nodes_.push_back(rootNode);

      bool isLeftLeaf =
          (nodeInfos[0].leftType == NODE_TYPE_LEAF) ? true : false;
      bool isRightLeaf =
          (nodeInfos[0].rightType == NODE_TYPE_LEAF) ? true : false;
      uint32_t midIndex = nodeInfos[0].childIndex;

      real3 leftBMin, leftBMax;
      real3 rightBMin, rightBMax;

      // printf("root: midIndex = %d, range (%d, %d), flag = %d/%d\n", midIndex,
      // 0, n-1, isLeftLeaf, isRightLeaf);

      size_t leftChildIndex = (size_t)(-1), rightChildIndex = (size_t)(-1);

#if defined(_OPENMP) && !defined(_MSC_VER)

      if (solids->isDoublePrecisionPos) {

#pragma omp parallel shared(leftChildIndex, rightChildIndex, keys, nodeInfos)
        {
#pragma omp single
          {
            leftChildIndex = BuildTreeRecursive32OMP(
                nodes_, leftBMin, leftBMax, solids->numVertsPerSolid, keys,
                nodeInfos, solids->dvertices, solids->indices, midIndex, 0,
                midIndex, isLeftLeaf, 0);
          }

#pragma omp single
          {
            rightChildIndex = BuildTreeRecursive32OMP(
                nodes_, rightBMin, rightBMax, solids->numVertsPerSolid, keys,
                nodeInfos, solids->dvertices, solids->indices, midIndex + 1,
                midIndex + 1, n - 1, isRightLeaf, 0);
          }
        }
#pragma omp barrier

      } else {

#pragma omp parallel shared(leftChildIndex, rightChildIndex, keys, nodeInfos)
        {
#pragma omp single
          {
            leftChildIndex = BuildTreeRecursive32OMP(
                nodes_, leftBMin, leftBMax, solids->numVertsPerSolid, keys,
                nodeInfos, solids->vertices, solids->indices, midIndex, 0,
                midIndex, isLeftLeaf, 0);
          }

#pragma omp single
          {
            rightChildIndex = BuildTreeRecursive32OMP(
                nodes_, rightBMin, rightBMax, solids->numVertsPerSolid, keys,
                nodeInfos, solids->vertices, solids->indices, midIndex + 1,
                midIndex + 1, n - 1, isRightLeaf, 0);
          }
        }
#pragma omp barrier
      }

      assert(leftChildIndex != (size_t)(-1));
      assert(rightChildIndex != (size_t)(-1));

#else

      if (solids->isDoublePrecisionPos) {

        leftChildIndex = BuildTreeRecursive32(
            nodes_, leftBMin, leftBMax, solids->numVertsPerSolid, keys,
            nodeInfos, solids->dvertices, solids->indices, midIndex, 0,
            midIndex, isLeftLeaf, 0);
        rightChildIndex = BuildTreeRecursive32(
            nodes_, rightBMin, rightBMax, solids->numVertsPerSolid, keys,
            nodeInfos, solids->dvertices, solids->indices, midIndex + 1,
            midIndex + 1, n - 1, isRightLeaf, 0);

      } else {

        leftChildIndex = BuildTreeRecursive32(
            nodes_, leftBMin, leftBMax, solids->numVertsPerSolid, keys,
            nodeInfos, solids->vertices, solids->indices, midIndex, 0, midIndex,
            isLeftLeaf, 0);
        rightChildIndex = BuildTreeRecursive32(
            nodes_, rightBMin, rightBMax, solids->numVertsPerSolid, keys,
            nodeInfos, solids->vertices, solids->indices, midIndex + 1,
            midIndex + 1, n - 1, isRightLeaf, 0);
      }

      assert(leftChildIndex != (size_t)(-1));
      assert(rightChildIndex != (size_t)(-1));

#endif

      // printf("  leftbmin = %f, %f, %f\n", leftBMin[0], leftBMin[1],
      // leftBMin[2]);
      // printf("  leftbmax = %f, %f, %f\n", leftBMax[0], leftBMax[1],
      // leftBMax[2]);
      // printf("  rightbmin = %f, %f, %f\n", rightBMin[0], rightBMin[1],
      // rightBMin[2]);
      // printf("  rightbmax = %f, %f, %f\n", rightBMax[0], rightBMax[1],
      // rightBMax[2]);
      // printf("  leftIndex = %d, rightIndex = %d\n", leftChildIndex,
      // rightChildIndex);
      // printf("  isLeaf = %d, %d\n", isLeftLeaf, isRightLeaf);

      real3 rootBMin, rootBMax;

      MergeBoundingBox(rootBMin, rootBMax, leftBMin, leftBMax, rightBMin,
                       rightBMax);

      rootNode.axis = GetSplitAxis(keys[0].code);
      rootNode.flag = 0; // = branch

      rootNode.data[0] = leftChildIndex;
      rootNode.data[1] = rightChildIndex;

      rootNode.bmin[0] = rootBMin[0];
      rootNode.bmin[1] = rootBMin[1];
      rootNode.bmin[2] = rootBMin[2];

      rootNode.bmax[0] = rootBMax[0];
      rootNode.bmax[1] = rootBMax[1];
      rootNode.bmax[2] = rootBMax[2];

      // @atomic update
      nodes_[0] = rootNode;

      t.end();
      // printf("root: midIndex = %d, range (%d, %d), flag = %d/%d\n", midIndex,
      // 0, n-1, isLeftLeaf, isRightLeaf);
      // printf("node.total = %d", nodes_.size());
      // printf("[%d] -> (%d, %d)\n", 0, leftChildIndex, rightChildIndex);
      trace("[LSGL] [6:Construct final AABB tree: %d msec\n", (int)t.msec());

      trace("[LSGL]   bmin = %f, %f, %f\n", rootBMin[0], rootBMin[1],
            rootBMin[2]);
      trace("[LSGL]   bmax = %f, %f, %f\n", rootBMax[0], rootBMax[1],
            rootBMax[2]);
    }

    {
      // Store sorted indices.
      assert(indices_.size() == keys.size());
      for (size_t i = 0; i < keys.size(); i++) {
        indices_[i] = keys[i].index;
      }
    }
  }

  // Tree will be null if input triangle count == 0.
  if (!nodes_.empty()) {
    // 0 = root node.
    real3 bmin(&nodes_[0].bmin[0]);
    real3 bmax(&nodes_[0].bmax[0]);
    trace("[SolidAccel] bound min = (%f, %f, %f)\n", bmin[0], bmin[1], bmin[2]);
    trace("[SolidAccel] bound max = (%f, %f, %f)\n", bmax[0], bmax[1], bmax[2]);
  }

  trace("[SolidAccel] # of nodes = %lu\n", nodes_.size());

  // Store pointer for later use.
  solids_ = solids;

  return true;
}

bool SolidAccel::Dump(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "[SolidAccel] Cannot write a file: %s\n", filename);
    return false;
  }

  unsigned long long numNodes = nodes_.size();
  assert(nodes_.size() > 0);

  unsigned long long numIndices = indices_.size();

  int r;
  r = fwrite(&numNodes, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&nodes_.at(0), sizeof(SolidNode), numNodes, fp);
  assert(r == numNodes);

  r = fwrite(&numIndices, sizeof(unsigned long long), 1, fp);
  assert(r == 1);

  r = fwrite(&indices_.at(0), sizeof(unsigned int), numIndices, fp);
  assert(r == numIndices);

  fclose(fp);

  return true;
}

namespace {

const int kMaxStackDepth = 512;

inline bool IntersectRayAABB(real &tminOut, // [out]
                             real &tmaxOut, // [out]
                             real maxT, real const bmin[3], real const bmax[3],
                             real3 rayOrg, real3 rayInvDir, int rayDirSign[3]) {
  real tmin, tmax;

  const real min_x = rayDirSign[0] ? bmax[0] : bmin[0];
  const real min_y = rayDirSign[1] ? bmax[1] : bmin[1];
  const real min_z = rayDirSign[2] ? bmax[2] : bmin[2];
  const real max_x = rayDirSign[0] ? bmin[0] : bmax[0];
  const real max_y = rayDirSign[1] ? bmin[1] : bmax[1];
  const real max_z = rayDirSign[2] ? bmin[2] : bmax[2];

  // X
  const double tmin_x = (min_x - rayOrg[0]) * rayInvDir[0];
  const double tmax_x = (max_x - rayOrg[0]) * rayInvDir[0];

  // Y
  const double tmin_y = (min_y - rayOrg[1]) * rayInvDir[1];
  const double tmax_y = (max_y - rayOrg[1]) * rayInvDir[1];

  tmin = (tmin_x > tmin_y) ? tmin_x : tmin_y;
  tmax = (tmax_x < tmax_y) ? tmax_x : tmax_y;

  // Z
  const double tmin_z = (min_z - rayOrg[2]) * rayInvDir[2];
  const double tmax_z = (max_z - rayOrg[2]) * rayInvDir[2];

  tmin = (tmin > tmin_z) ? tmin : tmin_z;
  tmax = (tmax < tmax_z) ? tmax : tmax_z;

  //
  // Hit include (tmin == tmax) edge case(hit 2D plane).
  //
  if ((tmax > 0.0) && (tmin <= tmax) && (tmin <= maxT)) {

    // printf("tmin, tmax = %f, %f\n", tmin, tmax);
    tminOut = tmin;
    tmaxOut = tmax;

    return true;
  }

  return false; // no hit
}

bool TestLeafNode(Intersection &isect, // [inout]
                  const SolidNode &node,
                  const std::vector<unsigned int> &indices, const Solid *solids,
                  const Ray &ray) {

  bool hit = false;

  unsigned int numSolids = node.data[0];
  unsigned int offset = node.data[1];

  real t = isect.t; // current hit distance

  double3 rayOrg;
  rayOrg[0] = ray.origin()[0];
  rayOrg[1] = ray.origin()[1];
  rayOrg[2] = ray.origin()[2];

  double3 rayDir;
  rayDir[0] = ray.direction()[0];
  rayDir[1] = ray.direction()[1];
  rayDir[2] = ray.direction()[2];

  for (unsigned int i = 0; i < numSolids; i++) {

    int primIdx = indices[i + offset];

    // self-intersection check
    if (primIdx == ray.prev_prim_id) {
      continue;
    }

    int numVertsPerSolid = solids->numVertsPerSolid;
    double3 vtx[16]; // assume numVertsPerSolid < 16

    for (int j = 0; j < numVertsPerSolid; j++) {
      int f = solids->indices[numVertsPerSolid * primIdx + j];

      if (solids->isDoublePrecisionPos) {
        vtx[j][0] = solids->dvertices[3 * f + 0];
        vtx[j][1] = solids->dvertices[3 * f + 1];
        vtx[j][2] = solids->dvertices[3 * f + 2];
      } else {
        vtx[j][0] = solids->vertices[3 * f + 0];
        vtx[j][1] = solids->vertices[3 * f + 1];
        vtx[j][2] = solids->vertices[3 * f + 2];
      }
    }

    Intersection isects[2]; // isects[0].t < isects[1].t

    bool hit = false;

    if (numVertsPerSolid == 5) {
      hit = IntersectPyramidD(rayOrg, rayDir, vtx, isects);
    } else if (numVertsPerSolid == 6) {
      hit = IntersectPrismD(rayOrg, rayDir, vtx, isects);
    } else if (numVertsPerSolid == 8) {
      hit = IntersectHexaD(rayOrg, rayDir, vtx, isects);
    }

    if (hit) {
      // @todo { Record bot hentering and leaving point. }
      if ((isects[0].t >= 0.0) && (isects[0].t < t)) {
        // Update isect state.
        isect.t = isects[0].t;
        isect.u = isects[0].u;
        isect.v = isects[0].v;
        isect.prim_id = primIdx;
        isect.subface_id = 0; // fixme
        isect.position = isects[0].position;
        isect.normal = isects[0].normal;
        float p[3] = {isect.position[0], isect.position[1], isect.position[2]};
        float d[8];
        interpolate(d, p, numVertsPerSolid, vtx, 0);
        isect.d0 = d[0];
        isect.d1 = d[1];
        isect.d2 = d[2];
        isect.d3 = d[3];
        isect.d4 = d[4];
        isect.d5 = d[5];
        isect.d6 = d[6];
        isect.d7 = d[7];
        isect.f0 = numVertsPerSolid * primIdx + 0;
        isect.f1 = numVertsPerSolid * primIdx + 1;
        isect.f2 = numVertsPerSolid * primIdx + 2;
        isect.f3 = numVertsPerSolid * primIdx + 3;
        isect.f4 = numVertsPerSolid * primIdx + 4;
        isect.f5 = numVertsPerSolid * primIdx + 5;
        isect.f6 = numVertsPerSolid * primIdx + 6;
        isect.f7 = numVertsPerSolid * primIdx + 7;
        t = isect.t;
        hit = true;
      }
      if ((isects[1].t >= 0.0) && (isects[1].t < t)) {
        // Update isect state.
        isect.t = isects[1].t;
        isect.u = isects[1].u;
        isect.v = isects[1].v;
        isect.prim_id = primIdx;
        isect.subface_id = 0; // fixme
        isect.position = isects[1].position;
        isect.normal = isects[1].normal;
        float p[3] = {isect.position[0], isect.position[1], isect.position[2]};
        float d[8];
        interpolate(d, p, numVertsPerSolid, vtx, 0);
        isect.d0 = d[0];
        isect.d1 = d[1];
        isect.d2 = d[2];
        isect.d3 = d[3];
        isect.d4 = d[4];
        isect.d5 = d[5];
        isect.d6 = d[6];
        isect.d7 = d[7];
        isect.f0 = numVertsPerSolid * primIdx + 0;
        isect.f1 = numVertsPerSolid * primIdx + 1;
        isect.f2 = numVertsPerSolid * primIdx + 2;
        isect.f3 = numVertsPerSolid * primIdx + 3;
        isect.f4 = numVertsPerSolid * primIdx + 4;
        isect.f5 = numVertsPerSolid * primIdx + 5;
        isect.f6 = numVertsPerSolid * primIdx + 6;
        isect.f7 = numVertsPerSolid * primIdx + 7;
        t = isect.t;
        hit = true;
      }
    }
  }

  return hit;
}

void BuildIntersection(Intersection &isect, const Solid *solids, Ray &ray) {

  // position and normal are already computed.

  isect.geometric = isect.normal;
}

} // namespace

bool SolidAccel::Traverse(Intersection &isect, Ray &ray) const {
  real hitT = REAL_MAX; // far = no hit.

#if ENABLE_TRAVERSAL_STATISTICS
  // @todo { multi-thread safe. }
  traversalStats_.numRays++;
#endif

  int nodeStackIndex = 0;
  int nodeStack[kMaxStackDepth];
  nodeStack[0] = 0;

  // Init isect info as no hit
  isect.t = hitT;
  isect.u = 0.0;
  isect.v = 0.0;
  isect.prim_id = (unsigned int)(-1);

  int dirSign[3];
  dirSign[0] = ray.direction()[0] < 0.0 ? 1 : 0;
  dirSign[1] = ray.direction()[1] < 0.0 ? 1 : 0;
  dirSign[2] = ray.direction()[2] < 0.0 ? 1 : 0;

  // @fixme { Check edge case; i.e., 1/0 }
  real3 rayInvDir;
  rayInvDir[0] = 1.0 / ray.direction()[0];
  rayInvDir[1] = 1.0 / ray.direction()[1];
  rayInvDir[2] = 1.0 / ray.direction()[2];

  real3 rayOrg;
  rayOrg[0] = ray.origin()[0];
  rayOrg[1] = ray.origin()[1];
  rayOrg[2] = ray.origin()[2];

  real minT, maxT;
  while (nodeStackIndex >= 0) {
    int index = nodeStack[nodeStackIndex];
    const SolidNode &node = nodes_[index];

    nodeStackIndex--;

    bool hit = IntersectRayAABB(minT, maxT, hitT, node.bmin, node.bmax, rayOrg,
                                rayInvDir, dirSign);

    if (node.flag == 0) { // branch node

      if (hit) {

#if ENABLE_TRAVERSAL_STATISTICS
        // @todo { multi-thread safe. }
        traversalStats_.numNodeTraversals++;
#endif

        int orderNear = dirSign[node.axis];
        int orderFar = 1 - orderNear;

        // Traverse near first.
        nodeStack[++nodeStackIndex] = node.data[orderFar];
        nodeStack[++nodeStackIndex] = node.data[orderNear];
      }

    } else { // leaf node

#if ENABLE_TRAVERSAL_STATISTICS
      // @todo { multi-thread safe. }
      traversalStats_.numLeafTests++;
#endif

      if (hit) {

#if ENABLE_TRAVERSAL_STATISTICS
        // @todo { multi-thread safe. }
        traversalStats_.numPrimIsectTests += node.data[0];
#endif

        if (TestLeafNode(isect, node, indices_, solids_, ray)) {
          hitT = isect.t;
        }
      }
    }
  }

  assert(nodeStackIndex < kMaxStackDepth);

  if (isect.t < REAL_MAX) {
    BuildIntersection(isect, solids_, ray);
    return true;
  }

  return false;
}

void SolidAccel::BoundingBox(double bmin[3], double bmax[3]) const {
  if (nodes_.empty()) {
    bmin[0] = (std::numeric_limits<double>::max)();
    bmin[1] = (std::numeric_limits<double>::max)();
    bmin[2] = (std::numeric_limits<double>::max)();
    bmax[0] = -(std::numeric_limits<double>::max)();
    bmax[1] = -(std::numeric_limits<double>::max)();
    bmax[2] = -(std::numeric_limits<double>::max)();
  } else {
    bmin[0] = nodes_[0].bmin[0];
    bmin[1] = nodes_[0].bmin[1];
    bmin[2] = nodes_[0].bmin[2];
    bmax[0] = nodes_[0].bmax[0];
    bmax[1] = nodes_[0].bmax[1];
    bmax[2] = nodes_[0].bmax[2];
  }
}

void SolidAccel::ResetTraversalStatistics() const {
  traversalStats_ = SolidTraversalStatistics();
}

void SolidAccel::ReportTraversalStatistics() const {
#if ENABLE_TRAVERSAL_STATISTICS
  double numRays = traversalStats_.numRays;
  printf("[LSGL] SolidAccel | # of rays              : %lld\n",
         traversalStats_.numRays);
  printf("[LSGL] SolidAccel | # of leaf node tests   : %lld(avg %f)\n",
         traversalStats_.numLeafTests, traversalStats_.numLeafTests / numRays);
  printf("[LSGL] SolidAccel | # of node traversals   : %lld(avg %f)\n",
         traversalStats_.numNodeTraversals,
         traversalStats_.numNodeTraversals / numRays);
  printf("[LSGL] SolidAccel | # of prim isect tests  : %lld(avg %f)\n",
         traversalStats_.numPrimIsectTests,
         traversalStats_.numPrimIsectTests / numRays);
#endif
}
