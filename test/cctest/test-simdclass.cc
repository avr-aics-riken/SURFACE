#include "gtest/gtest.h"

#include <cmath>

// SPARC HPC

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10

#include <emmintrin.h>

// Simple approximate floor.
int inline fasterfloor( const double x ) {
  // assume x is positive
  return (int)x;
}

TEST(SIMDClassTest, DFloor) {
  EXPECT_DOUBLE_EQ(floor(2.8), fasterfloor(2.8));
  EXPECT_DOUBLE_EQ(floor(2.0), fasterfloor(2.0));
  EXPECT_DOUBLE_EQ(floor(1.0), fasterfloor(1.0));
  EXPECT_DOUBLE_EQ(floor(0.0), fasterfloor(0.0));
  EXPECT_DOUBLE_EQ(floor(0.5), fasterfloor(0.5));
}

#endif

// x86 SSE
#if defined(__SSE2__) || (_M_IX86_FP >= 2) 

#include "glsl_simdclass.h"

namespace {

// Simple approximate floor.
int inline fasterfloor( const double x ) {
  if (x >= 0) {
    return (int)x;
  }

  int y = (int)x;
  if (std::abs(x - y) <= std::numeric_limits<double>::epsilon()) {
    // Do nothing.
  } else {
    y = y - 1;
  }

  return y;
}

int inline fasterfloor( const float x ) {
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

} // namespace

TEST(SIMDClassTest, FSet) {
  float a = 0.1;
  float b = 0.2;
  float c = 0.4;
  float d = 0.124;

  __vec4_f v(a, b, c, d);

  EXPECT_FLOAT_EQ(a, v.u.v[0]);
  EXPECT_FLOAT_EQ(b, v.u.v[1]);
  EXPECT_FLOAT_EQ(c, v.u.v[2]);
  EXPECT_FLOAT_EQ(d, v.u.v[3]);
}

TEST(SIMDClassTest, FAdd) {
  float a = 0.1;
  float b = 0.2;
  float c = 0.4;
  float d = 0.124;
  float e = 0.2;
  float f = 0.4;
  float g = 0.7;
  float h = 0.326;

  __vec4_f v0(a, b, c, d);
  __vec4_f v1(e, f, g, h);
  __vec4_f v2 = v0 + v1;

  EXPECT_FLOAT_EQ(a+e, v2.u.v[0]);
  EXPECT_FLOAT_EQ(b+f, v2.u.v[1]);
  EXPECT_FLOAT_EQ(c+g, v2.u.v[2]);
  EXPECT_FLOAT_EQ(d+h, v2.u.v[3]);
}

TEST(SIMDClassTest, FSub) {
  float a = 0.1;
  float b = 0.2;
  float c = 0.4;
  float d = 0.124;
  float e = 0.2;
  float f = 0.4;
  float g = 0.7;
  float h = 0.326;

  __vec4_f v0(a, b, c, d);
  __vec4_f v1(e, f, g, h);
  __vec4_f v2 = v0 - v1;

  EXPECT_FLOAT_EQ(a-e, v2.u.v[0]);
  EXPECT_FLOAT_EQ(b-f, v2.u.v[1]);
  EXPECT_FLOAT_EQ(c-g, v2.u.v[2]);
  EXPECT_FLOAT_EQ(d-h, v2.u.v[3]);
}

TEST(SIMDClassTest, FMul) {
  float a = 0.1;
  float b = 0.2;
  float c = 0.4;
  float d = 0.124;
  float e = 0.2;
  float f = 0.4;
  float g = 0.7;
  float h = 0.326;

  __vec4_f v0(a, b, c, d);
  __vec4_f v1(e, f, g, h);
  __vec4_f v2 = v0 * v1;

  EXPECT_FLOAT_EQ(a*e, v2.u.v[0]);
  EXPECT_FLOAT_EQ(b*f, v2.u.v[1]);
  EXPECT_FLOAT_EQ(c*g, v2.u.v[2]);
  EXPECT_FLOAT_EQ(d*h, v2.u.v[3]);
}

TEST(SIMDClassTest, FDiv) {
  float a = 0.1;
  float b = 0.2;
  float c = 0.4;
  float d = 0.124;
  float e = 0.2;
  float f = 0.4;
  float g = 0.7;
  float h = 0.326;

  __vec4_f v0(a, b, c, d);
  __vec4_f v1(e, f, g, h);
  __vec4_f v2 = v0 / v1;

  EXPECT_FLOAT_EQ(a/e, v2.u.v[0]);
  EXPECT_FLOAT_EQ(b/f, v2.u.v[1]);
  EXPECT_FLOAT_EQ(c/g, v2.u.v[2]);
  EXPECT_FLOAT_EQ(d/h, v2.u.v[3]);
}
 

TEST(SIMDClassTest, FLerp) {
  float a = 0.1;
  float b = 0.2;
  float c = 0.4;
  float d = 0.124;
  float e = 0.2;
  float f = 0.4;
  float g = 0.7;
  float h = 0.326;
  float t = 0.48;

  __vec4_f v0(a, b, c, d);
  __vec4_f v1(e, f, g, h);
  __vec4_f v2 = (1.0 - t) * v0 + t * v1;

  EXPECT_FLOAT_EQ((1.0-t)*a+t*e, v2.u.v[0]);
  EXPECT_FLOAT_EQ((1.0-t)*b+t*f, v2.u.v[1]);
  EXPECT_FLOAT_EQ((1.0-t)*c+t*g, v2.u.v[2]);
  EXPECT_FLOAT_EQ((1.0-t)*d+t*h, v2.u.v[3]);
}

TEST(SIMDClassTest, FFloor) {
  EXPECT_FLOAT_EQ(floorf(2.8), fasterfloor(2.8));
  EXPECT_FLOAT_EQ(floorf(-2.8), fasterfloor(-2.8));
  EXPECT_FLOAT_EQ(floorf(2.0), fasterfloor(2.0));
  EXPECT_FLOAT_EQ(floorf(1.0), fasterfloor(1.0));
  EXPECT_FLOAT_EQ(floorf(0.0), fasterfloor(0.0));
  EXPECT_FLOAT_EQ(floorf(-0.0), fasterfloor(-0.0));
  EXPECT_FLOAT_EQ(floorf(-1.0), fasterfloor(-1.0));
  EXPECT_FLOAT_EQ(floorf(-2.0), fasterfloor(-2.0));
  EXPECT_FLOAT_EQ(floorf(0.5), fasterfloor(0.5));
  EXPECT_FLOAT_EQ(floorf(-0.5), fasterfloor(-0.5));
}

TEST(SIMDClassTest, FRepeat) {
  for (float u = -10.0; u <= 10.0; u+= 0.001) {
    EXPECT_FLOAT_EQ(u - floorf(u), u - fasterfloor(u));
  }
}

TEST(SIMDClassTest, DSet) {
  double a = 0.1;
  double b = 0.2;
  double c = 0.4;
  double d = 0.124;

  __vec4_d v(a, b, c, d);

  EXPECT_DOUBLE_EQ(a, v.u.v[0]);
  EXPECT_DOUBLE_EQ(b, v.u.v[1]);
  EXPECT_DOUBLE_EQ(c, v.u.v[2]);
  EXPECT_DOUBLE_EQ(d, v.u.v[3]);
}

TEST(SIMDClassTest, DAdd) {
  double a = 0.1;
  double b = 0.2;
  double c = 0.4;
  double d = 0.124;
  double e = 0.2;
  double f = 0.4;
  double g = 0.7;
  double h = 0.326;

  __vec4_d v0(a, b, c, d);
  __vec4_d v1(e, f, g, h);
  __vec4_d v2 = v0 + v1;

  EXPECT_DOUBLE_EQ(a+e, v2.u.v[0]);
  EXPECT_DOUBLE_EQ(b+f, v2.u.v[1]);
  EXPECT_DOUBLE_EQ(c+g, v2.u.v[2]);
  EXPECT_DOUBLE_EQ(d+h, v2.u.v[3]);
}

TEST(SIMDClassTest, DSub) {
  double a = 0.1;
  double b = 0.2;
  double c = 0.4;
  double d = 0.124;
  double e = 0.2;
  double f = 0.4;
  double g = 0.7;
  double h = 0.326;

  __vec4_d v0(a, b, c, d);
  __vec4_d v1(e, f, g, h);
  __vec4_d v2 = v0 - v1;

  EXPECT_DOUBLE_EQ(a-e, v2.u.v[0]);
  EXPECT_DOUBLE_EQ(b-f, v2.u.v[1]);
  EXPECT_DOUBLE_EQ(c-g, v2.u.v[2]);
  EXPECT_DOUBLE_EQ(d-h, v2.u.v[3]);
}

TEST(SIMDClassTest, DMul) {
  double a = 0.1;
  double b = 0.2;
  double c = 0.4;
  double d = 0.124;
  double e = 0.2;
  double f = 0.4;
  double g = 0.7;
  double h = 0.326;

  __vec4_d v0(a, b, c, d);
  __vec4_d v1(e, f, g, h);
  __vec4_d v2 = v0 * v1;

  EXPECT_DOUBLE_EQ(a*e, v2.u.v[0]);
  EXPECT_DOUBLE_EQ(b*f, v2.u.v[1]);
  EXPECT_DOUBLE_EQ(c*g, v2.u.v[2]);
  EXPECT_DOUBLE_EQ(d*h, v2.u.v[3]);
}

TEST(SIMDClassTest, DDiv) {
  double a = 0.1;
  double b = 0.2;
  double c = 0.4;
  double d = 0.124;
  double e = 0.2;
  double f = 0.4;
  double g = 0.7;
  double h = 0.326;

  __vec4_d v0(a, b, c, d);
  __vec4_d v1(e, f, g, h);
  __vec4_d v2 = v0 / v1;

  EXPECT_DOUBLE_EQ(a/e, v2.u.v[0]);
  EXPECT_DOUBLE_EQ(b/f, v2.u.v[1]);
  EXPECT_DOUBLE_EQ(c/g, v2.u.v[2]);
  EXPECT_DOUBLE_EQ(d/h, v2.u.v[3]);
}
 

TEST(SIMDClassTest, DLerp) {
  double a = 0.1;
  double b = 0.2;
  double c = 0.4;
  double d = 0.124;
  double e = 0.2;
  double f = 0.4;
  double g = 0.7;
  double h = 0.326;
  double t = 0.48;

  __vec4_d v0(a, b, c, d);
  __vec4_d v1(e, f, g, h);
  __vec4_d v2 = (1.0 - t) * v0 + t * v1;

  EXPECT_DOUBLE_EQ((1.0-t)*a+t*e, v2.u.v[0]);
  EXPECT_DOUBLE_EQ((1.0-t)*b+t*f, v2.u.v[1]);
  EXPECT_DOUBLE_EQ((1.0-t)*c+t*g, v2.u.v[2]);
  EXPECT_DOUBLE_EQ((1.0-t)*d+t*h, v2.u.v[3]);
}

TEST(SIMDClassTest, DFloor) {
  EXPECT_DOUBLE_EQ(floor(2.8), fasterfloor(2.8));
  EXPECT_DOUBLE_EQ(floor(-2.8), fasterfloor(-2.8));
  EXPECT_DOUBLE_EQ(floor(2.0), fasterfloor(2.0));
  EXPECT_DOUBLE_EQ(floor(1.0), fasterfloor(1.0));
  EXPECT_DOUBLE_EQ(floor(0.0), fasterfloor(0.0));
  EXPECT_DOUBLE_EQ(floor(-0.0), fasterfloor(-0.0));
  EXPECT_DOUBLE_EQ(floor(-1.0), fasterfloor(-1.0));
  EXPECT_DOUBLE_EQ(floor(-2.0), fasterfloor(-2.0));
  EXPECT_DOUBLE_EQ(floor(0.5), fasterfloor(0.5));
  EXPECT_DOUBLE_EQ(floor(-0.5), fasterfloor(-0.5));
}

TEST(SIMDClassTest, DRepeat) {
  for (double u = -10.0; u <= 10.0; u+= 0.001) {
    EXPECT_DOUBLE_EQ(u - floor(u), u - fasterfloor(u));
  }
}

#endif
