/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>

#if defined(LSGL_OPTIMIZE_GLSL)
#include <emmintrin.h>
#include "../glsl/glsl_simdclass.h"
#endif

#include "render_texture.h"

#ifndef M_PI
#define M_PI 3.141592
#endif

using namespace lsgl::render;

namespace {

inline double dlerp(double t, double a, double b) {
  return (1.0 - t) * a + t * b;
}
}

#if defined(LSGL_OPTIMIZE_GLSL) && defined(__sparc__) && defined(__HPC_ACE__)
#define FORCEINLINE __attribute__((always_inline))

namespace {

FORCEINLINE __vec4_d
    vclamp(const __vec4_d &v, const __vec4_d &low, const __vec4_d &high) {
  return vmin(vmax(v, low), high);
}

FORCEINLINE __vec4_d
    vlerp4(const __vec4_d &t, const __vec4_d &a, const __vec4_d &b) {
  const __vec4_d vone(1.0);
  return (vone - t) * a + t * b;
}

FORCEINLINE __vec2_d
    vlerp2(const __vec2_d &t, const __vec2_d a, const __vec2_d &b) {
  const __vec2_d vone(1.0);
  return (vone - t) * a + t * b;
}
}

#elif defined(LSGL_OPTIMIZE_GLSL) &&                                           \
    (defined(__SSE2__) || (_M_IX86_FP >= 2)) // _M_IX86_FP 2 = VS /arch:SSE2

#define FORCEINLINE __attribute__((always_inline))

FORCEINLINE __vec4_d
    vclamp(const __vec4_d &v, const __vec4_d &low, const __vec4_d &high) {
  return vmin(vmax(v, low), high);
}

FORCEINLINE __vec4_f
    vclampf(const __vec4_f &v, const __vec4_f &low, const __vec4_f &high) {
  return vmin(vmax(v, low), high);
}

FORCEINLINE __vec4_d
    vlerp4(const __vec4_d &t, const __vec4_d &a, const __vec4_d &b) {
  const __vec4_d vone(1.0);
  return (vone - t) * a + t * b;
}

FORCEINLINE __vec4_f
    vlerp4f(const __vec4_f &t, const __vec4_f &a, const __vec4_f &b) {
  const __vec4_f vone(1.0f);
  return (vone - t) * a + t * b;
}

FORCEINLINE __vec2_d
    vlerp2(const __vec2_d &t, const __vec2_d &a, const __vec2_d &b) {
  const __vec2_d vone(1.0);
  return (vone - t) * a + t * b;
}

#endif

namespace {

inline float Lerp(float t, float a, float b) { return (1.f - t) * a + t * b; }

inline int Clamp(int v, int low, int high) {
  return std::min(std::max(v, low), high);
}

// inline bool myisnan(float a) {
//  volatile float d = a;
//  return d != d;
//}

inline int fasterfloorf(const float x) {
  if (x >= 0) {
    return (int)x;
  }

  int y = (int)x;
  if (std::abs(x - y) <= std::numeric_limits<float>::epsilon()) {
    // Do nothing.
  } else {
    y = y - 1;
  }

  return y;
}

inline void FilterByte(float *rgba, const unsigned char *image, int i00,
                       int i10, int i01, int i11, float w[4], // weight
                       int stride) {
  unsigned char texel[4][4];

  const float inv = 1.0f / 255.0f;
  if (stride == 4) {

    for (int i = 0; i < 4; i++) {
      texel[0][i] = image[i00 + i];
      texel[1][i] = image[i10 + i];
      texel[2][i] = image[i01 + i];
      texel[3][i] = image[i11 + i];
    }

    for (int i = 0; i < 4; i++) {
      rgba[i] = w[0] * texel[0][i] + w[1] * texel[1][i] + w[2] * texel[2][i] +
                w[3] * texel[3][i];
      // normalize.
      rgba[i] *= inv;
    }

  } else {

    for (int i = 0; i < stride; i++) {
      texel[0][i] = image[i00 + i];
      texel[1][i] = image[i10 + i];
      texel[2][i] = image[i01 + i];
      texel[3][i] = image[i11 + i];
    }

    for (int i = 0; i < stride; i++) {
      rgba[i] = w[0] * texel[0][i] + w[1] * texel[1][i] + w[2] * texel[2][i] +
                w[3] * texel[3][i];
      // normalize.
      rgba[i] *= inv;
    }
  }

  if (stride < 4) {
    rgba[3] = 1.0; // In OpenGL, alpha = 1.0 by default
  }

  if (stride == 1) {
    // splat
    rgba[1] = rgba[2] = rgba[3] = rgba[0];
  }
}

inline void FilterFloat(float *rgba, const float *image, int i00, int i10,
                        int i01, int i11, float w[4], // weight
                        int stride) {
  float texel[4][4];

  if (stride == 4) {

    for (int i = 0; i < 4; i++) {
      texel[0][i] = image[i00 + i];
      texel[1][i] = image[i10 + i];
      texel[2][i] = image[i01 + i];
      texel[3][i] = image[i11 + i];
    }

    for (int i = 0; i < 4; i++) {
      rgba[i] = w[0] * texel[0][i] + w[1] * texel[1][i] + w[2] * texel[2][i] +
                w[3] * texel[3][i];
    }

  } else {

    for (int i = 0; i < stride; i++) {
      texel[0][i] = image[i00 + i];
      texel[1][i] = image[i10 + i];
      texel[2][i] = image[i01 + i];
      texel[3][i] = image[i11 + i];
    }

    for (int i = 0; i < stride; i++) {
      rgba[i] = w[0] * texel[0][i] + w[1] * texel[1][i] + w[2] * texel[2][i] +
                w[3] * texel[3][i];
    }
  }

  if (stride < 4) {
    rgba[3] = 1.0; // In OpenGL, alpha = 1.0 by default
  }

  if (stride == 1) {
    // splat
    rgba[1] = rgba[2] = rgba[3] = rgba[0];
  }
}

inline void FilterDouble(float *rgba, const double *image, int i00, int i10,
                         int i01, int i11, float w[4], // weight
                         int stride) {
  double texel[4][4];

  if (stride == 4) {

    for (int i = 0; i < 4; i++) {
      texel[0][i] = image[i00 + i];
      texel[1][i] = image[i10 + i];
      texel[2][i] = image[i01 + i];
      texel[3][i] = image[i11 + i];
    }

    for (int i = 0; i < 4; i++) {
      rgba[i] = w[0] * texel[0][i] + w[1] * texel[1][i] + w[2] * texel[2][i] +
                w[3] * texel[3][i];
    }

  } else {

    for (int i = 0; i < stride; i++) {
      texel[0][i] = image[i00 + i];
      texel[1][i] = image[i10 + i];
      texel[2][i] = image[i01 + i];
      texel[3][i] = image[i11 + i];
    }

    for (int i = 0; i < stride; i++) {
      rgba[i] = w[0] * texel[0][i] + w[1] * texel[1][i] + w[2] * texel[2][i] +
                w[3] * texel[3][i];
    }
  }

  if (stride < 4) {
    rgba[3] = 1.0; // In OpenGL, alpha = 1.0 by default
  }

  if (stride == 1) {
    // splat
    rgba[1] = rgba[2] = rgba[3] = rgba[0];
  }
}

} // namespace

void Texture2D::fetch(float *rgba, float u, float v, bool minFiltering, bool magFiltering) const {
  // @todo { Support wrap mode. }

  float sx = fasterfloorf(u);
  float sy = fasterfloorf(v);

  float uu = u - sx;
  float vv = v - sy;

  // clamp
  uu = std::max(uu, 0.0f);
  uu = std::min(uu, 1.0f);
  vv = std::max(vv, 0.0f);
  vv = std::min(vv, 1.0f);

  float px = (m_width - 1) * uu;
  float py = (m_height - 1) * vv;

  if (magFiltering == false && minFiltering == false) {

    int x0 = (int)px;
    int y0 = (int)py;

    int stride = m_components;

    int i00 = stride * (y0 * m_width + x0);

    if (m_format == FORMAT_BYTE) {
      const unsigned char *image = m_image;
      for (int i = 0; i < stride; i++) { 
        rgba[i] = image[i00 + i] / 255.0f;
      }

      if (stride < 4) {
        rgba[3] = 1.0; // In OpenGL, alpha = 1.0 by default
      }

      if (stride == 1) {
        // splat
        rgba[1] = rgba[2] = rgba[3] = rgba[0];
      }
    } else if (m_format == FORMAT_FLOAT32) {
      const float *image = reinterpret_cast<const float *>(m_image);
      for (int i = 0; i < stride; i++) { 
        rgba[i] = image[i00 + i];
      }

      if (stride < 4) {
        rgba[3] = 1.0; // In OpenGL, alpha = 1.0 by default
      }

      if (stride == 1) {
        // splat
        rgba[1] = rgba[2] = rgba[3] = rgba[0];
      }
    } else if (m_format == FORMAT_FLOAT64) {
      const double *image = reinterpret_cast<const double *>(m_image);
      for (int i = 0; i < stride; i++) { 
        rgba[i] = image[i00 + i];
      }

      if (stride < 4) {
        rgba[3] = 1.0; // In OpenGL, alpha = 1.0 by default
      }

      if (stride == 1) {
        // splat
        rgba[1] = rgba[2] = rgba[3] = rgba[0];
      }
    } else { // unknown
    }

  } else {

    int x0 = (int)px;
    int y0 = (int)py;
    int x1 = ((x0 + 1) >= m_width) ? (m_width - 1) : (x0 + 1);
    int y1 = ((y0 + 1) >= m_height) ? (m_height - 1) : (y0 + 1);

    float dx = px - (float)x0;
    float dy = py - (float)y0;

    float w[4];

    w[0] = (1.0f - dx) * (1.0 - dy);
    w[1] = (1.0f - dx) * (dy);
    w[2] = (dx) * (1.0 - dy);
    w[3] = (dx) * (dy);

    int stride = m_components;

    int i00 = stride * (y0 * m_width + x0);
    int i01 = stride * (y0 * m_width + x1);
    int i10 = stride * (y1 * m_width + x0);
    int i11 = stride * (y1 * m_width + x1);

    if (m_format == FORMAT_BYTE) {
      FilterByte(rgba, m_image, i00, i10, i01, i11, w, stride);
    } else if (m_format == FORMAT_FLOAT32) {
      FilterFloat(rgba, reinterpret_cast<const float *>(m_image), i00, i10, i01,
                  i11, w, stride);
    } else if (m_format == FORMAT_FLOAT64) {
      FilterDouble(rgba, reinterpret_cast<const double *>(m_image), i00, i10, i01,
                   i11, w, stride);
    } else { // unknown
    }
  }
}

void Texture2D::fetchD(float *rgba0, float *rgba1, float *rgba2, float u,
                       float v, bool minFiltering, bool magFiltering) const {
  // @todo { optimize! }

  // fetch (i, j)
  fetch(rgba0, u, v, minFiltering, magFiltering);

  // fetch (i+1, j)
  float u1 = u + m_invWidth;
  fetch(rgba1, u1, v, minFiltering, magFiltering);

  // fetch (i, j+1)
  float v1 = v + m_invHeight;
  fetch(rgba2, u, v1, minFiltering, magFiltering);
}

inline float D(int x, int y, int z, int comp, int idx, int dim[3],
               const float density[]) {
  x = Clamp(x, 0, dim[0] - 1);
  y = Clamp(y, 0, dim[1] - 1);
  z = Clamp(z, 0, dim[2] - 1);
  return density[comp * (z * dim[0] * dim[1] + y * dim[0] + x) + idx];
}

inline float ADDR3D_UB(int x, int y, int z, int comp, int idx, int dim[3],
                       const unsigned char density[]) {
  return density[comp * (z * dim[0] * dim[1] + y * dim[0] + x) + idx] / 255.5f;
}

inline float ADDR3D(int x, int y, int z, int comp, int idx, int dim[3],
                    const float density[]) {
  return density[comp * (z * dim[0] * dim[1] + y * dim[0] + x) + idx];
}

inline double ADDR3DD(int x, int y, int z, int comp, int idx, int dim[3],
                      const double density[]) {
  return density[comp * (z * dim[0] * dim[1] + y * dim[0] + x) + idx];
}

inline float ADDR3DF(double x, double y, double z, double comp, double idx,
                     double dim[3], const float density[]) {
  return density[(
      unsigned long long)(comp * (z * dim[0] * dim[1] + y * dim[0] + x) + idx)];
}

#if 0
void Texture3D::fetch(float *rgba, float u, float v, float r) const {

  rgba[0] = rgba[1] = rgba[2] = rgba[3] = 0.0f;

#if 0
  // NaN check.
  if (myisnan(u) || myisnan(v) || myisnan(r)) {
    int stride = m_components;
    for (int i = 0; i < stride; i++) {
      rgba[i] = 0.0f;
    }
    return;
  }
#endif

  // @fixme { REPEAT only }
  float u01 = u - fasterfloorf(u);
  float v01 = v - fasterfloorf(v);
  float r01 = r - fasterfloorf(r);

//#if 0
#if defined(LSGL_OPTIMIZE_GLSL) && defined(__sparc__) && defined(__HPC_ACE__)
  // SIMD optimized code path

  int dim[3];
  dim[0] = width();
  dim[1] = height();
  dim[2] = depth();

  __vec4_d vzero(0.0);
  __vec4_d vone(1.0);
  __vec4_d vdim(width(), height(), depth(), 0.0);
  __vec4_d vcoord(u01, v01, r01, 0.0);

  __vec4_d vox = vcoord * (vdim - vone);

  const float *density = image();

  const int vx0 = (int)(vox.u.v[0]);
  const int vy0 = (int)(vox.u.v[1]);
  const int vz0 = (int)(vox.u.v[2]);
  const int vx = Clamp(vx0, 0, dim[0] - 1);
  const int vy = Clamp(vy0, 0, dim[1] - 1);
  const int vz = Clamp(vz0, 0, dim[2] - 1);
  const int vx1 = std::min(vx0 + 1, dim[0] - 1);
  const int vy1 = std::min(vy0 + 1, dim[1] - 1);
  const int vz1 = std::min(vz0 + 1, dim[2] - 1);
  const __vec4_d vxyz0(vx0, vy0, vz0, 0.0);
  const __vec4_d vxyz = vclamp(vxyz0, vzero, vdim - vone);

  double dx = vox.u.v[0] - vxyz.u.v[0];
  double dy = vox.u.v[1] - vxyz.u.v[1];
  double dz = vox.u.v[2] - vxyz.u.v[2];

  int comps = components();

  double vals[4] = { 0.0, 0.0, 0.0, 0.0 }; // up to 4 components.
  for (int i = 0; i < comps; i++) {

// 1 lerp: 1 sub, 1 add, 2 mul : 4 op
// 7 lerps: 7 sub, 7 add, 14 mul: 28 ops
#if 0

    // Trilinearly interpolate density values to compute local density
    double d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
                      D(vx + 1, vy, vz, comps, i, dim, density));
    double d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
                      D(vx + 1, vy + 1, vz, comps, i, dim, density));
    double d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
                      D(vx + 1, vy, vz + 1, comps, i, dim, density));
    double d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
                      D(vx + 1, vy + 1, vz + 1, comps, i, dim, density));
#endif
//double d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
//                  D(vx + 1, vy, vz, comps, i, dim, density));
//double d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
//                  D(vx + 1, vy + 1, vz, comps, i, dim, density));
//double d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
//                  D(vx + 1, vy, vz + 1, comps, i, dim, density));
//double d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
//                  D(vx + 1, vy + 1, vz + 1, comps, i, dim, density));
#if 0
    const float D000 = ADDR3D(vx, vy, vz, comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy, vz, comps, i, dim, density);
    const float D010 = ADDR3D(vx, vy1, vz, comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz, comps, i, dim, density);
    const float D001 = ADDR3D(vx, vy, vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy, vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx, vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D(vx1, vy1, vz1, comps, i, dim, density);
#else
    const __vec4_d addr_vx(vx, vx1, vx, vx1);
    const __vec4_d addr_vy(vy, vy, vy1, vy1);
    const __vec4_d addr_vz0(vz);
    const __vec4_d addr_vz1(vz1);
    const __vec4_d addr_comp(comps);
    const __vec4_d addr_i(i);
    const __vec4_d addr_dim0(dim[0]);
    const __vec4_d addr_dim1(dim[1]);
    const __vec4_d addr_dim2(dim[2]);

    const __vec4_d addr0 = addr_comp * (addr_vz0 * addr_dim0 * addr_dim1 +
                                        addr_vy * addr_dim0 + addr_vx) + addr_i;
    const __vec4_d addr1 = addr_comp * (addr_vz1 * addr_dim0 * addr_dim1 +
                                        addr_vy * addr_dim0 + addr_vx) + addr_i;
    const float D000 = density[(unsigned long long) addr0.u.v[0]];
    const float D100 = density[(unsigned long long) addr0.u.v[1]];
    const float D010 = density[(unsigned long long) addr0.u.v[2]];
    const float D110 = density[(unsigned long long) addr0.u.v[3]];
    const float D001 = density[(unsigned long long) addr1.u.v[0]];
    const float D101 = density[(unsigned long long) addr1.u.v[1]];
    const float D011 = density[(unsigned long long) addr1.u.v[2]];
    const float D111 = density[(unsigned long long) addr1.u.v[3]];
#endif

    // (d00, d01, d10, d11)
    const __vec4_d vsrc0(D000, D001, D010, D011);
    const __vec4_d vsrc1(D100, D101, D110, D111);

    __vec4_d vd4 = vlerp4(__vec4_d(dx), vsrc0, vsrc1);

    //double d0 = Lerp(dy, d00, d10);
    //double d1 = Lerp(dy, d01, d11);
    const __vec2_d vd01 = vlerp2(__vec2_d(dy), vd4.u.v0, vd4.u.v1);
    //const __vec2_d vd01 = __vec2_d(dy) + vd4.u.v0 + vd4.u.v1;

    double dret = dlerp(dz, vd01.u.v[0], vd01.u.v[1]);
    //printf("dret = %f\n", dret);
    vals[i] = dret;
  }

  if (comps == 1) { // scalar
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }
#elif defined(LSGL_OPTIMIZE_GLSL) && (defined(__SSE2__) || (_M_IX86_FP >= 2))
// SIMD optimized code path

#if 1 // double ver
  int dim[3];
  dim[0] = width();
  dim[1] = height();
  dim[2] = depth();

  __vec4_d vzero(0.0);
  __vec4_d vone(1.0);
  __vec4_d vdim(width(), height(), depth(), 0.0);
  __vec4_d vcoord(u01, v01, r01, 0.0);

  __vec4_d vox = vcoord * (vdim - vone);

  const float *density = image();

  const int vx0 = (int)(vox.u.v[0]);
  const int vy0 = (int)(vox.u.v[1]);
  const int vz0 = (int)(vox.u.v[2]);
  const int vx = Clamp(vx0, 0, dim[0] - 1);
  const int vy = Clamp(vy0, 0, dim[1] - 1);
  const int vz = Clamp(vz0, 0, dim[2] - 1);
  const int vx1 = std::min(vx0 + 1, dim[0] - 1);
  const int vy1 = std::min(vy0 + 1, dim[1] - 1);
  const int vz1 = std::min(vz0 + 1, dim[2] - 1);
  const __vec4_d vxyz0(vx0, vy0, vz0, 0.0);
  const __vec4_d vxyz = vclamp(vxyz0, vzero, vdim - vone);

  double dx = vox.u.v[0] - vxyz.u.v[0];
  double dy = vox.u.v[1] - vxyz.u.v[1];
  double dz = vox.u.v[2] - vxyz.u.v[2];

  int comps = components();

  double vals[4] = { 0.0, 0.0, 0.0, 0.0 }; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const float D000 = ADDR3D(vx, vy, vz, comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy, vz, comps, i, dim, density);
    const float D010 = ADDR3D(vx, vy1, vz, comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz, comps, i, dim, density);
    const float D001 = ADDR3D(vx, vy, vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy, vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx, vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D(vx1, vy1, vz1, comps, i, dim, density);

    // (d00, d01, d10, d11)
    const __vec4_d vsrc0(D000, D001, D010, D011);
    const __vec4_d vsrc1(D100, D101, D110, D111);

    __vec4_d vd4 = vlerp4(__vec4_d(dx), vsrc0, vsrc1);
    const __vec2_d vd01 = vlerp2(__vec2_d(dy), vd4.u.v0, vd4.u.v1);
    //const __vec2_d vd01 = __vec2_d(dy) + vd4.u.v0 + vd4.u.v1;

    double dret = dlerp(dz, vd01.u.v[0], vd01.u.v[1]);
    //printf("dret = %f\n", dret);
    vals[i] = dret;
  }

  if (comps == 1) { // scalar
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }
#else // float ver

  int dim[3];
  dim[0] = width();
  dim[1] = height();
  dim[2] = depth();

  __vec4_f vzero(0.0f);
  __vec4_f vone(1.0f);
  __vec4_f vdim(width(), height(), depth(), 0.0f);
  __vec4_f vcoord(u01, v01, r01, 0.0f);

  __vec4_f vox = vcoord * (vdim - vone);

  const float *density = image();

  const int vx0 = (int)(vox.u.v[0]);
  const int vy0 = (int)(vox.u.v[1]);
  const int vz0 = (int)(vox.u.v[2]);
  const int vx = Clamp(vx0, 0, dim[0] - 1);
  const int vy = Clamp(vy0, 0, dim[1] - 1);
  const int vz = Clamp(vz0, 0, dim[2] - 1);
  const int vx1 = std::min(vx0 + 1, dim[0] - 1);
  const int vy1 = std::min(vy0 + 1, dim[1] - 1);
  const int vz1 = std::min(vz0 + 1, dim[2] - 1);
  const __vec4_f vxyz0(vx0, vy0, vz0, 0.0);
  const __vec4_f vxyz = vclamp(vxyz0, vzero, vdim - vone);

  double dx = vox.u.v[0] - vxyz.u.v[0];
  double dy = vox.u.v[1] - vxyz.u.v[1];
  double dz = vox.u.v[2] - vxyz.u.v[2];

  int comps = components();

  double vals[4] = { 0.0, 0.0, 0.0, 0.0 }; // up to 4 components.
  for (int i = 0; i < comps; i++) {

// 1 lerp: 1 sub, 1 add, 2 mul : 4 op
// 7 lerps: 7 sub, 7 add, 14 mul: 28 ops
#if 0

    // Trilinearly interpolate density values to compute local density
    double d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
                      D(vx + 1, vy, vz, comps, i, dim, density));
    double d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
                      D(vx + 1, vy + 1, vz, comps, i, dim, density));
    double d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
                      D(vx + 1, vy, vz + 1, comps, i, dim, density));
    double d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
                      D(vx + 1, vy + 1, vz + 1, comps, i, dim, density));
#endif
    //double d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
    //                  D(vx + 1, vy, vz, comps, i, dim, density));
    //double d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
    //                  D(vx + 1, vy + 1, vz, comps, i, dim, density));
    //double d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
    //                  D(vx + 1, vy, vz + 1, comps, i, dim, density));
    //double d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
    //                  D(vx + 1, vy + 1, vz + 1, comps, i, dim, density));
    const float D000 = ADDR3D(vx, vy, vz, comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy, vz, comps, i, dim, density);
    const float D010 = ADDR3D(vx, vy1, vz, comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz, comps, i, dim, density);
    const float D001 = ADDR3D(vx, vy, vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy, vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx, vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D(vx1, vy1, vz1, comps, i, dim, density);

    // (d00, d01, d10, d11)
    const __vec4_f vsrc0(D000, D001, D010, D011);
    const __vec4_f vsrc1(D100, D101, D110, D111);

    __vec4_f vd4 = vlerp4f(__vec4_f(dx), vsrc0, vsrc1);

    //double d0 = Lerp(dy, d00, d10);
    //double d1 = Lerp(dy, d01, d11);
    //const __vec2_d vd01 = vlerp2(__vec2_d(dy), vd4.u.v0, vd4.u.v1);
    //const __vec2_d vd01 = __vec2_d(dy) + vd4.u.v0 + vd4.u.v1;
    //double dret= dlerp(dz, vd01.u.v[0], vd01.u.v[1]);
    float d0 = Lerp(dy, vd4.u.v[0], vd4.u.v[2]);
    float d1 = Lerp(dy, vd4.u.v[1], vd4.u.v[3]);
    float dret = Lerp(dz, d0, d1);

    //printf("dret = %f\n", dret);
    vals[i] = dret;
  }

  if (comps == 1) { // scalar
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }

#endif

#else

  int dim[3];
  dim[0] = width();
  dim[1] = height();
  dim[2] = depth();

  float vox[3];
  //vox[0] = u01 * dim[0] - .5f;
  //vox[1] = v01 * dim[1] - .5f;
  //vox[2] = r01 * dim[2] - .5f;
  vox[0] = u01 * (dim[0] - 1.0f);
  vox[1] = v01 * (dim[1] - 1.0f);
  vox[2] = r01 * (dim[2] - 1.0f);

  const float *density = image();

  int vx = (int)(vox[0]), vy = (int)(vox[1]), vz = (int)(vox[2]);
  vx = Clamp(vx, 0, dim[0] - 1);
  vy = Clamp(vy, 0, dim[1] - 1);
  vz = Clamp(vz, 0, dim[2] - 1);
  int vx1 = std::min(vx + 1, dim[0] - 1);
  int vy1 = std::min(vy + 1, dim[1] - 1);
  int vz1 = std::min(vz + 1, dim[2] - 1);

  float dx = vox[0] - vx, dy = vox[1] - vy, dz = vox[2] - vz;
  //assert(dx <= 1.0);
  //assert(dy <= 1.0);
  //assert(dz <= 1.0);
  //assert(dx >= 0.0f);
  //assert(dy >= 0.0f);
  //assert(dz >= 0.0f);

  int comps = components();
  assert(comps >= 0);
  assert(comps <= 4);

  float vals[4] = { 0.0f, 0.0f, 0.0f, 0.0f }; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const float D000 = ADDR3D(vx, vy, vz, comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy, vz, comps, i, dim, density);
    const float D010 = ADDR3D(vx, vy1, vz, comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz, comps, i, dim, density);
    const float D001 = ADDR3D(vx, vy, vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy, vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx, vy1, vz1, comps, i, dim, density);
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
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }
#endif

}
#endif

void FilterTexture3DByte(float *rgba, const Texture3D *tex, float u, float v,
                         float r) {

  rgba[0] = rgba[1] = rgba[2] = rgba[3] = 0.0f;

#if 0 // disable for faster performance.
      // NaN check.
  if (myisnan(u) || myisnan(v) || myisnan(r)) {
    int stride = m_components;
    for (int i = 0; i < stride; i++) {
      rgba[i] = 0.0f;
    }
    return;
  }
#endif

  // @fixme { REPEAT only }
  float u01 = u - fasterfloorf(u);
  float v01 = v - fasterfloorf(v);
  float r01 = r - fasterfloorf(r);

  int dim[3];
  dim[0] = tex->width;
  dim[1] = tex->height;
  dim[2] = tex->depth;

  float vox[3];
  // vox[0] = u01 * dim[0] - .5f;
  // vox[1] = v01 * dim[1] - .5f;
  // vox[2] = r01 * dim[2] - .5f;
  vox[0] = u01 * (dim[0] - 1.0f);
  vox[1] = v01 * (dim[1] - 1.0f);
  vox[2] = r01 * (dim[2] - 1.0f);

  const unsigned char *density = reinterpret_cast<unsigned char *>(tex->image);

  int vx = (int)(vox[0]), vy = (int)(vox[1]), vz = (int)(vox[2]);
  vx = Clamp(vx, 0, dim[0] - 1);
  vy = Clamp(vy, 0, dim[1] - 1);
  vz = Clamp(vz, 0, dim[2] - 1);
  int vx1 = std::min(vx + 1, dim[0] - 1);
  int vy1 = std::min(vy + 1, dim[1] - 1);
  int vz1 = std::min(vz + 1, dim[2] - 1);

  float dx = vox[0] - vx, dy = vox[1] - vy, dz = vox[2] - vz;
  // assert(dx <= 1.0);
  // assert(dy <= 1.0);
  // assert(dz <= 1.0);
  // assert(dx >= 0.0f);
  // assert(dy >= 0.0f);
  // assert(dz >= 0.0f);

  int comps = tex->components;

  float vals[4] = {0.0f, 0.0f, 0.0f, 0.0f}; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const float D000 = ADDR3D_UB(vx, vy, vz, comps, i, dim, density);
    const float D100 = ADDR3D_UB(vx1, vy, vz, comps, i, dim, density);
    const float D010 = ADDR3D_UB(vx, vy1, vz, comps, i, dim, density);
    const float D110 = ADDR3D_UB(vx1, vy1, vz, comps, i, dim, density);
    const float D001 = ADDR3D_UB(vx, vy, vz1, comps, i, dim, density);
    const float D101 = ADDR3D_UB(vx1, vy, vz1, comps, i, dim, density);
    const float D011 = ADDR3D_UB(vx, vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D_UB(vx1, vy1, vz1, comps, i, dim, density);

    // Trilinearly interpolate density values to compute local density
    // float d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
    //                 D(vx + 1, vy, vz, comps, i, dim, density));
    // float d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
    //                 D(vx + 1, vy + 1, vz, comps, i, dim, density));
    // float d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
    //                 D(vx + 1, vy, vz + 1, comps, i, dim, density));
    // float d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
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
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }
}

void FilterTexture3DFloat(float *rgba, const Texture3D *tex, float u, float v,
                          float r) {

  rgba[0] = rgba[1] = rgba[2] = rgba[3] = 0.0f;

#if 0
  // NaN check.
  if (myisnan(u) || myisnan(v) || myisnan(r)) {
    int stride = m_components;
    for (int i = 0; i < stride; i++) {
      rgba[i] = 0.0f;
    }
    return;
  }
#endif

  // @fixme { REPEAT only }
  float u01 = u - fasterfloorf(u);
  float v01 = v - fasterfloorf(v);
  float r01 = r - fasterfloorf(r);

  int dim[3];
  dim[0] = tex->width;
  dim[1] = tex->height;
  dim[2] = tex->depth;

  float vox[3];
  // vox[0] = u01 * dim[0] - .5f;
  // vox[1] = v01 * dim[1] - .5f;
  // vox[2] = r01 * dim[2] - .5f;
  vox[0] = u01 * (dim[0] - 1.0f);
  vox[1] = v01 * (dim[1] - 1.0f);
  vox[2] = r01 * (dim[2] - 1.0f);

  const float *density = reinterpret_cast<float *>(tex->image);

  int vx = (int)(vox[0]), vy = (int)(vox[1]), vz = (int)(vox[2]);
  vx = Clamp(vx, 0, dim[0] - 1);
  vy = Clamp(vy, 0, dim[1] - 1);
  vz = Clamp(vz, 0, dim[2] - 1);
  int vx1 = std::min(vx + 1, dim[0] - 1);
  int vy1 = std::min(vy + 1, dim[1] - 1);
  int vz1 = std::min(vz + 1, dim[2] - 1);

  float dx = vox[0] - vx, dy = vox[1] - vy, dz = vox[2] - vz;
  // assert(dx <= 1.0);
  // assert(dy <= 1.0);
  // assert(dz <= 1.0);
  // assert(dx >= 0.0f);
  // assert(dy >= 0.0f);
  // assert(dz >= 0.0f);

  int comps = tex->components;

  float vals[4] = {0.0f, 0.0f, 0.0f, 0.0f}; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const float D000 = ADDR3D(vx, vy, vz, comps, i, dim, density);
    const float D100 = ADDR3D(vx1, vy, vz, comps, i, dim, density);
    const float D010 = ADDR3D(vx, vy1, vz, comps, i, dim, density);
    const float D110 = ADDR3D(vx1, vy1, vz, comps, i, dim, density);
    const float D001 = ADDR3D(vx, vy, vz1, comps, i, dim, density);
    const float D101 = ADDR3D(vx1, vy, vz1, comps, i, dim, density);
    const float D011 = ADDR3D(vx, vy1, vz1, comps, i, dim, density);
    const float D111 = ADDR3D(vx1, vy1, vz1, comps, i, dim, density);

    // Trilinearly interpolate density values to compute local density
    // float d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
    //                 D(vx + 1, vy, vz, comps, i, dim, density));
    // float d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
    //                 D(vx + 1, vy + 1, vz, comps, i, dim, density));
    // float d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
    //                 D(vx + 1, vy, vz + 1, comps, i, dim, density));
    // float d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
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
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }
}

void FilterTexture3DDouble(float *rgba, const Texture3D *tex, float u, float v,
                           float r) {

  rgba[0] = rgba[1] = rgba[2] = rgba[3] = 0.0f;

#if 0
  // NaN check.
  if (myisnan(u) || myisnan(v) || myisnan(r)) {
    int stride = m_components;
    for (int i = 0; i < stride; i++) {
      rgba[i] = 0.0f;
    }
    return;
  }
#endif

  // @fixme { REPEAT only }
  float u01 = u - fasterfloorf(u);
  float v01 = v - fasterfloorf(v);
  float r01 = r - fasterfloorf(r);

  int dim[3];
  dim[0] = tex->width;
  dim[1] = tex->height;
  dim[2] = tex->depth;

  float vox[3];
  // vox[0] = u01 * dim[0] - .5f;
  // vox[1] = v01 * dim[1] - .5f;
  // vox[2] = r01 * dim[2] - .5f;
  vox[0] = u01 * (dim[0] - 1.0f);
  vox[1] = v01 * (dim[1] - 1.0f);
  vox[2] = r01 * (dim[2] - 1.0f);

  const double *density = reinterpret_cast<double *>(tex->image);

  int vx = (int)(vox[0]), vy = (int)(vox[1]), vz = (int)(vox[2]);
  vx = Clamp(vx, 0, dim[0] - 1);
  vy = Clamp(vy, 0, dim[1] - 1);
  vz = Clamp(vz, 0, dim[2] - 1);
  int vx1 = std::min(vx + 1, dim[0] - 1);
  int vy1 = std::min(vy + 1, dim[1] - 1);
  int vz1 = std::min(vz + 1, dim[2] - 1);

  double dx = vox[0] - vx, dy = vox[1] - vy, dz = vox[2] - vz;
  // assert(dx <= 1.0);
  // assert(dy <= 1.0);
  // assert(dz <= 1.0);
  // assert(dx >= 0.0f);
  // assert(dy >= 0.0f);
  // assert(dz >= 0.0f);

  int comps = tex->components;

  double vals[4] = {0.0f, 0.0f, 0.0f, 0.0f}; // up to 4 components.
  for (int i = 0; i < comps; i++) {

    const double D000 = ADDR3DD(vx, vy, vz, comps, i, dim, density);
    const double D100 = ADDR3DD(vx1, vy, vz, comps, i, dim, density);
    const double D010 = ADDR3DD(vx, vy1, vz, comps, i, dim, density);
    const double D110 = ADDR3DD(vx1, vy1, vz, comps, i, dim, density);
    const double D001 = ADDR3DD(vx, vy, vz1, comps, i, dim, density);
    const double D101 = ADDR3DD(vx1, vy, vz1, comps, i, dim, density);
    const double D011 = ADDR3DD(vx, vy1, vz1, comps, i, dim, density);
    const double D111 = ADDR3DD(vx1, vy1, vz1, comps, i, dim, density);

    // Trilinearly interpolate density values to compute local density
    // float d00 = Lerp(dx, D(vx, vy, vz, comps, i, dim, density),
    //                 D(vx + 1, vy, vz, comps, i, dim, density));
    // float d10 = Lerp(dx, D(vx, vy + 1, vz, comps, i, dim, density),
    //                 D(vx + 1, vy + 1, vz, comps, i, dim, density));
    // float d01 = Lerp(dx, D(vx, vy, vz + 1, comps, i, dim, density),
    //                 D(vx + 1, vy, vz + 1, comps, i, dim, density));
    // float d11 = Lerp(dx, D(vx, vy + 1, vz + 1, comps, i, dim, density),
    //                 D(vx + 1, vy + 1, vz + 1, comps, i, dim, density));
    double d00 = dlerp(dx, D000, D100);
    double d10 = dlerp(dx, D010, D110);
    double d01 = dlerp(dx, D001, D101);
    double d11 = dlerp(dx, D011, D111);
    double d0 = dlerp(dy, d00, d10);
    double d1 = dlerp(dy, d01, d11);

    double d = dlerp(dz, d0, d1);
    vals[i] = d;
  }

  if (comps == 1) { // scalar
    rgba[0] = vals[0];
    rgba[1] = vals[0];
    rgba[2] = vals[0];
    rgba[3] = vals[0];
  } else {
    rgba[0] = vals[0];
    rgba[1] = vals[1];
    rgba[2] = vals[2];
    rgba[3] = vals[3];
  }
}
