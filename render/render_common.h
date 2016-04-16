/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __RENDER_COMMON_H__
#define __RENDER_COMMON_H__

#include <cmath>
#include <cfloat>

#define RENDER_USE_DOUBLE_PRECISION (0)

namespace lsgl {
namespace render {

#if RENDER_USE_DOUBLE_PRECISION
typedef double real; // double don't work at this time
#define REAL_EPSILON  DBL_EPSILON
#define REAL_MAX      DBL_MAX
#else
typedef float real;
#define REAL_EPSILON  FLT_EPSILON
#define REAL_MAX      FLT_MAX
#endif

struct real3 {
  real3() {}
  real3(real xx, real yy, real zz) {
    v[0] = xx;
    v[1] = yy;
    v[2] = zz;
  }
  explicit real3(real *p) {
    v[0] = p[0];
    v[1] = p[1];
    v[2] = p[2];
  }

  real3 operator*(real f) const { return real3(x() * f, y() * f, z() * f); }
  real3 operator-(const real3 &f2) const {
    return real3(x() - f2.x(), y() - f2.y(), z() - f2.z());
  }
  real3 operator*(const real3 &f2) const {
    return real3(x() * f2.x(), y() * f2.y(), z() * f2.z());
  }
  real3 operator+(const real3 &f2) const {
    return real3(x() + f2.x(), y() + f2.y(), z() + f2.z());
  }
  real3 operator/(const real3 &f2) const {
    return real3(x() / f2.x(), y() / f2.y(), z() / f2.z());
  }
  real3 operator/(real f) const {
    return real3(x() / f, y() / f, z() / f);
  }
  real operator[](int i) const {
    return v[i];
  }
  real &operator[](int i) {
    return v[i];
  }

  real length() const { return sqrt(x() * x() + y() * y() + z() * z()); }

  real3 &normalize() {

    double length2 = x() * x() + y() * y() + z() * z();
    double length = 0.0;

    if (fabs(length2) > 1.0e-30) {
      length = 1.0 / sqrt(length2);
    }

    v[0] *= length;
    v[1] *= length;
    v[2] *= length;

    return *this;
  }

  inline real x() const {
    return v[0];
  }

  inline real y() const {
    return v[1];
  }

  inline real z() const {
    return v[2];
  }

  real v[3];
  //real pad;  // for alignment
};

inline real3 operator*(real f, const real3 &v) {
  return real3(v.x() * f, v.y() * f, v.z() * f);
}

struct double3 {
  double3() {}
  double3(double xx, double yy, double zz) {
    v[0] = xx;
    v[1] = yy;
    v[2] = zz;
  }
  explicit double3(double *p) {
    v[0] = p[0];
    v[1] = p[1];
    v[2] = p[2];
  }

  double3 operator*(double f) const { return double3(x() * f, y() * f, z() * f); }
  double3 operator-(const double3 &f2) const {
    return double3(x() - f2.x(), y() - f2.y(), z() - f2.z());
  }
  double3 operator*(const double3 &f2) const {
    return double3(x() * f2.x(), y() * f2.y(), z() * f2.z());
  }
  double3 operator+(const double3 &f2) const {
    return double3(x() + f2.x(), y() + f2.y(), z() + f2.z());
  }
  double3 operator/(const double3 &f2) const {
    return double3(x() / f2.x(), y() / f2.y(), z() / f2.z());
  }
  double3 operator/(double f) const {
    return double3(x() / f, y() / f, z() / f);
  }
  double operator[](int i) const { return v[i]; }
  double &operator[](int i) { return v[i]; }

  inline double x() const {
    return v[0];
  }

  inline double y() const {
    return v[1];
  }

  inline double z() const {
    return v[2];
  }


  double v[3];
};

inline double3 operator*(double f, const double3 &v) {
  return double3(v.x() * f, v.y() * f, v.z() * f);
}

inline real dot(const real3 &lhs, const real3 &rhs) {
  return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]) + (lhs[2] * rhs[2]);
}

inline real3 cross(const real3 &lhs, const real3 &rhs) {
  return real3(lhs[1] * rhs[2] - lhs[2] * rhs[1], // xyzzy
               lhs[2] * rhs[0] - lhs[0] * rhs[2], // yzxxz
               lhs[0] * rhs[1] - lhs[1] * rhs[0]  // zxyyx
               );
}

inline real3 normalize(const real3 &v) {
  real3 vn = v;
  vn.normalize();
  return vn;
}

} // namespace render
} // namespace lsgl

#endif // __RENDER_COMMON_H__
