/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __RENDER_COMMON_H__
#define __RENDER_COMMON_H__

#include <cmath>

#define RENDER_USE_DOUBLE_PRECISION (0)

namespace lsgl {
namespace render {

#if RENDER_USE_DOUBLE_PRECISION
typedef double real; // double don't work at this time
#else
typedef float real;
#endif

struct real3 {
  real3() {}
  real3(real xx, real yy, real zz) {
    x = xx;
    y = yy;
    z = zz;
  }
  real3(real *p) {
    x = p[0];
    y = p[1];
    z = p[2];
  }

  real3 operator*(real f) const { return real3(x * f, y * f, z * f); }
  real3 operator-(const real3 &f2) const {
    return real3(x - f2.x, y - f2.y, z - f2.z);
  }
  real3 operator*(const real3 &f2) const {
    return real3(x * f2.x, y * f2.y, z * f2.z);
  }
  real3 operator+(const real3 &f2) const {
    return real3(x + f2.x, y + f2.y, z + f2.z);
  }
  real3 operator/(const real3 &f2) const {
    return real3(x / f2.x, y / f2.y, z / f2.z);
  }
  real operator[](int i) const { return (&x)[i]; }
  real &operator[](int i) { return (&x)[i]; }

  real length() const { return sqrt(x * x + y * y + z * z); }

  real3 &normalize() {

    double length2 = x * x + y * y + z * z;
    double length = 0.0;

    if (fabs(length2) > 1.0e-30) {
      length = 1.0 / sqrt(length2);
    }

    x *= length;
    y *= length;
    z *= length;

    return *this;
  }

  real x, y, z;
  // real pad;  // for alignment
};

inline real3 operator*(real f, const real3 &v) {
  return real3(v.x * f, v.y * f, v.z * f);
}

struct double3 {
  double3() {}
  double3(double xx, double yy, double zz) {
    x = xx;
    y = yy;
    z = zz;
  }
  double3(double *p) {
    x = p[0];
    y = p[1];
    z = p[2];
  }

  double3 operator*(double f) const { return double3(x * f, y * f, z * f); }
  double3 operator-(const double3 &f2) const {
    return double3(x - f2.x, y - f2.y, z - f2.z);
  }
  double3 operator*(const double3 &f2) const {
    return double3(x * f2.x, y * f2.y, z * f2.z);
  }
  double3 operator+(const double3 &f2) const {
    return double3(x + f2.x, y + f2.y, z + f2.z);
  }
  double3 operator/(const double3 &f2) const {
    return double3(x / f2.x, y / f2.y, z / f2.z);
  }
  double operator[](int i) const { return (&x)[i]; }
  double &operator[](int i) { return (&x)[i]; }

  double x, y, z;
};

inline real dot(const real3 &lhs, const real3 &rhs) {
  return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]) + (lhs[2] * rhs[2]);
}

inline real3 cross(const real3 &lhs, const real3 &rhs) {
  return real3(lhs[1] * rhs[2] - lhs[2] * rhs[1], // xyzzy
               lhs[2] * rhs[0] - lhs[0] * rhs[2], // yzxxz
               lhs[0] * rhs[1] - lhs[1] * rhs[0]  // zxyyx
               );
}

} // namespace render
} // namespace lsgl

#endif // __RENDER_COMMON_H__
