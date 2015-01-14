/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_VECTOR3_H__
#define __LSGL_VECTOR3_H__

#include <cstddef> //std::size_t
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cassert>

#include "simd_util.h"

#ifdef _MSC_VER
#define FORCEINLINE __forceinline
#else
#define FORCEINLINE __attribute__((always_inline)) inline
#endif

namespace lsgl {
namespace render {

template <class T> class vector3t {
public:
  //-----------------------------------------------
  // type defines
  typedef T value_type;
  typedef const T &const_reference;
  typedef vector3t<T> this_type;

  typedef std::size_t size_type;

  typedef T *iterator;
  typedef const float *const_iterator;

public:
  static const size_type c_size = 3; // container size

public:
  // constructors and destructor
  FORCEINLINE vector3t() {}

  FORCEINLINE vector3t(value_type _x, value_type _y, value_type _z)
      : m0(_x), m1(_y), m2(_z) {}

  FORCEINLINE vector3t(const this_type &rhs) : m0(rhs.m0), m1(rhs.m1), m2(rhs.m2) {}

  explicit vector3t(const value_type rhs[c_size])
      : m0(rhs[0]), m1(rhs[1]), m2(rhs[2]) {}

  template <class X> FORCEINLINE vector3t(const vector3t<X> &rhs) {
    m0 = rhs[0];
    m1 = rhs[1];
    m2 = rhs[2];
  }

  ~vector3t() {}

  //-----------------------------------------------
  // inserters

  FORCEINLINE this_type &operator=(const this_type &rhs) {

    m0 = rhs.m0;
    m1 = rhs.m1;
    m2 = rhs.m2;

    return *this;
  }

  template <class X> FORCEINLINE this_type &operator=(const vector3t<X> &rhs) {

    m0 = rhs[0];
    m1 = rhs[1];
    m2 = rhs[2];

    return *this;
  }

  //-----------------------------------------------
  // operators

  FORCEINLINE this_type &negate() {
    m0 = -m0;
    m1 = -m1;
    m2 = -m2;
    return *this;
  }

#define DECLARE_OP_EQUAL(OP)                                                   \
  FORCEINLINE this_type &operator OP(const this_type &rhs) {                   \
    m0 OP rhs.m0;                                                              \
    m1 OP rhs.m1;                                                              \
    m2 OP rhs.m2;                                                              \
    return *this;                                                              \
  }

  DECLARE_OP_EQUAL(+= )
  DECLARE_OP_EQUAL(-= )
  DECLARE_OP_EQUAL(*= )
  DECLARE_OP_EQUAL(/= )

#undef DECLARE_OP_EQUAL

  FORCEINLINE this_type &operator*=(T rhs) {
    m0 *= rhs;
    m1 *= rhs;
    m2 *= rhs;

    return *this;
  }
  FORCEINLINE this_type &operator/=(T rhs) {
    m0 /= rhs;
    m1 /= rhs;
    m2 /= rhs;

    return *this;
  }
  //--------------------------------

  FORCEINLINE value_type &operator[](size_type i) { return element()[i]; }

  FORCEINLINE value_type operator[](size_type i) const { return element()[i]; }

  FORCEINLINE value_type &at(size_type i) {
    // if(c_size<=i){throw std::out_of_range(debug());}
    return element()[i];
  }

  FORCEINLINE const value_type &at(size_type i) const {
    // if(c_size<=i){throw std::out_of_range(debug());}
    return element()[i];
  }

  //-----------------------------------------------
  // utilities
  FORCEINLINE value_type length() const {
    return std::sqrt(sqr_length());
  }
  FORCEINLINE value_type sqr_length() const { return ((m0 * m0) + (m1 * m1) + (m2 * m2)); }

  FORCEINLINE this_type &normalize() {

    value_type length = sqr_length(); //||V||^2
                                      // if (length == T()) return *this;

    length = value_type(1) / sqrt(length);
    m0 *= length;
    m1 *= length;
    m2 *= length;

    return *this;
  }

  FORCEINLINE this_type &fast_normalize() {
    value_type length = sqr_length(); //||V||^2
    float4 v = vset1_f4(length);
    float4 vinv_len = vfastrsqrt_f4(v);
    float buf[4];
    vstoreu_f4(buf, vinv_len);
    float inv_len = buf[0];
    m0 *= inv_len;
    m1 *= inv_len;
    m2 *= inv_len;

    return *this;
  }

private:
  value_type m0, m1, m2;

  T *element() { return &m0; }
  const T *element() const { return &m0; }
};

template <class T> FORCEINLINE const vector3t<T> operator+(const vector3t<T> &rhs) {
  return rhs;
}
template <class T> FORCEINLINE vector3t<T> operator-(const vector3t<T> &rhs) {
  return vector3t<T>(rhs).negate();
}

#define DECLARE_OPERATOR(OP)                                                   \
  template <class T>                                                           \
  inline vector3t<T> operator OP(const vector3t<T> &lhs,                       \
                                 const vector3t<T> &rhs) {                     \
    return vector3t<T>(lhs) OP## = rhs;                                        \
  }
DECLARE_OPERATOR(+)
DECLARE_OPERATOR(-)
DECLARE_OPERATOR(*)
DECLARE_OPERATOR(/ )

#undef DECLARE_OPERATOR

template <class T>
inline vector3t<T> operator*(double lhs, const vector3t<T> &rhs) {
  return vector3t<T>(rhs) *= (T)lhs;
}

template <class T>
inline vector3t<T> operator*(const vector3t<T> &lhs, double rhs) {
  return vector3t<T>(lhs) *= rhs;
}

template <class T>
inline vector3t<T> operator/(const vector3t<T> &lhs, double rhs) {
  return vector3t<T>(lhs) /= rhs;
}

template <class T> inline T length(const vector3t<T> &rhs) {
  return rhs.length();
}

template <class T> inline T sqr_length(const vector3t<T> &rhs) {
  return rhs.sqr_length();
}

template <class T> inline vector3t<T> normalize(const vector3t<T> &rhs) {
  return vector3t<T>(rhs).normalize();
}

template <class T>
inline T dot(const vector3t<T> &lhs, const vector3t<T> &rhs) {
  return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]) + (lhs[2] * rhs[2]);
}

template <class T>
inline vector3t<T> cross(const vector3t<T> &lhs, const vector3t<T> &rhs) {
  return vector3t<T>(lhs[1] * rhs[2] - lhs[2] * rhs[1], // xyzzy
                     lhs[2] * rhs[0] - lhs[0] * rhs[2], // yzxxz
                     lhs[0] * rhs[1] - lhs[1] * rhs[0]  // zxyyx
                     );
}

template <class T>
inline bool operator==(const vector3t<T> &lhs, const vector3t<T> &rhs) {
  return (lhs[0] == rhs[0]) && (lhs[1] == rhs[1]) && (lhs[2] == rhs[2]);
}

template <class T>
inline bool operator!=(const vector3t<T> &lhs, const vector3t<T> &rhs) {
  return !(lhs == rhs);
}

typedef vector3t<float> vector3;
typedef vector3t<float> vector3f;
typedef vector3t<double> vector3d;

} // namespace render
} // namespace lsgl

#endif // __LSGL_VECTOR3_H__
