/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_SIMD_H__
#define __LSGL_SIMD_H__

#ifdef __cplusplus
extern "C" {
#endif
#if defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2
#include <xmmintrin.h>
#include <emmintrin.h>
#if defined(__SSSE3__)
#include <tmmintrin.h> // SSSE3
#endif
#if defined(__SSE4_1__)
#include <smmintrin.h> // _mm_blendv_ps
#endif
#if defined(__SSE4_2__)
#include <nmmintrin.h> // SSE4.2
#endif
#if defined(__AVX__)
#include <immintrin.h>
#endif
#else
// non non SSE
#endif

#ifdef __cplusplus
}
#endif

#include <float.h>

// ----------------------------------------------------------------------------
//
// Defines
//
// ----------------------------------------------------------------------------

#if defined(_MSC_VER)
#define ATTRIB_ALIGN_16
#else // Assume gcc
#define ATTRIB_ALIGN_16 __attribute__((aligned(16)))
#endif

// ----------------------------------------------------------------------------
//
// SIMD data types
//
// ----------------------------------------------------------------------------

#if defined(__SSE2__) || (_M_IX86_FP >= 2)
typedef __m128 float4;
typedef __m128i int4;
#if defined(__AVX__)
typedef __m256 float8;
typedef __m256d double4;
#endif
#else
typedef struct { float v[4]; } float4;

typedef struct { int v[4]; } int4;

typedef struct { unsigned int v[4]; } uint4;
#endif

// ----------------------------------------------------------------------------
//
// SIMD operations
//
// ----------------------------------------------------------------------------

#if defined(__SSE2__) || (_M_IX86_FP >= 2) // sse

// Float operations

#define vadd_f4(a, b) _mm_add_ps((a), (b))
#define vsub_f4(a, b) _mm_sub_ps((a), (b))
#define vmul_f4(a, b) _mm_mul_ps((a), (b))
#define vcmplt_f4(a, b) _mm_cmplt_ps((a), (b))
#define vcmple_f4(a, b) _mm_cmple_ps((a), (b))
#define vcmpgt_f4(a, b) _mm_cmpgt_ps((a), (b))
#define vcmpge_f4(a, b) _mm_cmpge_ps((a), (b))
#define vmax_f4(a, b) _mm_max_ps((a), (b))
#define vmin_f4(a, b) _mm_min_ps((a), (b))
#define vand_f4(a, b) _mm_and_ps((a), (b))
#define vor_f4(a, b) _mm_or_ps((a), (b))

#define vdot_f4(ax, ay, az, bx, by, bz)                                        \
  (_mm_add_ps(_mm_mul_ps((ax), (bx)),                                          \
              _mm_add_ps(_mm_mul_ps((ay), (by)), _mm_mul_ps((az), (bz)))))

#define vcross_f4(a, b, c, d)                                                  \
  (_mm_sub_ps(_mm_mul_ps((a), (c)), _mm_mul_ps((b), (d))))

#if defined(__SSE4_1__)
#define vsel_f4(a, b, m) _mm_blendv_ps((b), (a), (m))
#else
#define vsel_f4(a, b, m)                                                       \
  _mm_or_ps(_mm_and_ps((m), (a)), _mm_andnot_ps((m), (b)))
#endif

#define vset1_f4(a) _mm_set1_ps((a))
#define vset4_f4(a, b, c, d) _mm_set_ps((a), (b), (c), (d))
#define vsetr4_f4(a, b, c, d) _mm_setr_ps((a), (b), (c), (d))
#define vzero_f4() _mm_setzero_ps()
#define vone_f4() _mm_set1_ps(1.0f)
#define vshuffle_f4(a, b, i) _mm_shuffle_ps((a), (b), (i))
#define vloadu_f4(p) _mm_loadu_ps((p))
#define vstoreu_f4(p, v) _mm_storeu_ps((p), (v))

#define vconvtoi4_f4(v) _mm_cvttps_epi32((v))

#define vtranspose4_f4(x0, x1, x2, x3) _MM_TRANSPOSE4_PS(x0, x1, x2, x3)

#if defined(_MSC_VER)
#define vcasttoi4_f4(v) _mm_castps_si128((v))
#define vcasttof4_i4(v) _mm_castsi128_ps((v))
#else
#define vcasttoi4_f4(v) (__m128i)((v))
#define vcasttof4_i4(v) (__m128)((v))
#endif

static inline float4 vabs_f4(const float4 &v) {
  // fabs(x) = (asUInt(x) << 1) >> 1)
  return vcasttof4_i4(_mm_srli_epi32(_mm_slli_epi32(vcasttoi4_f4(v), 1), 1));
}

// v must be positive value.
static inline float4 vfastinv_f4(const float4 &v) {
  // Approximated division by reciprocal estimate + one round of Newton-Raphson
  // arg = 0 -> NaN
  // if (arg > -eps) arg = -eps
  // if (arg <  eps) arg =  eps
  const float4 vzero = vset1_f4(0.0f);
  const float4 vone = vset1_f4(1.0f);
  const float4 vnegone = vset1_f4(-1.0f);
  const float4 veps = vset1_f4(1.0e-6f);
  const float4 vabs = vabs_f4(v);
  const float4 filter_vabs = vmax_f4(veps, vabs);
  const float4 vsignmask = vcmpgt_f4(v, vzero);
  const float4 vsign = vsel_f4(vone, vnegone, vsignmask);
  const float4 filter_v = vmul_f4(vsign, filter_vabs);

  const float4 rcp = _mm_rcp_ps(filter_v);
  return _mm_sub_ps(_mm_add_ps(rcp, rcp),
                    _mm_mul_ps(_mm_mul_ps(rcp, rcp), filter_v));
}

static inline float4 vfastrsqrt_f4(const float4 &v) {
  const float4 _three = vset1_f4(3.0f);
  const float4 _half4 = vset1_f4(0.5f);
  const float4 approx = _mm_rsqrt_ps(v);
  const float4 muls = _mm_mul_ps(_mm_mul_ps(v, approx), approx);
  return _mm_mul_ps(_mm_mul_ps(_half4, approx), _mm_sub_ps(_three, muls));
}

// Assume v >= 0.0
static inline float4 vfastsqrt_f4(const float4 &v) {
  return _mm_mul_ps(v, vfastsqrt_f4(v));
}

// [xyzxyzxyzxyz] -> [xxxx, yyyy, zzzz]
static inline void vaos2soa_xyz_f4(float4 &x, float4 &y, float4 &z,
                                   const float *p) {
  __m128 x0y0z0x1 = _mm_load_ps(p + 0);
  __m128 y1z1x2y2 = _mm_load_ps(p + 4);
  __m128 z2x3y3z3 = _mm_load_ps(p + 8);
  __m128 x2y2x3y3 = _mm_shuffle_ps(y1z1x2y2, z2x3y3z3, _MM_SHUFFLE(2, 1, 3, 2));
  __m128 y0z0y1z1 = _mm_shuffle_ps(x0y0z0x1, y1z1x2y2, _MM_SHUFFLE(1, 0, 2, 1));
  x = _mm_shuffle_ps(x0y0z0x1, x2y2x3y3, _MM_SHUFFLE(2, 0, 3, 0)); // x0x1x2x3
  y = _mm_shuffle_ps(y0z0y1z1, x2y2x3y3, _MM_SHUFFLE(3, 1, 2, 0)); // y0y1y2y3
  z = _mm_shuffle_ps(y0z0y1z1, z2x3y3z3, _MM_SHUFFLE(3, 0, 3, 1)); // z0z1z2z3
}

// [xxxx,yyyy,zzzz] -> [xyzxyzxyzxyz]
static inline void vsoa2aos_xyz_f4(float *p, const float4 &x, const float4 &y,
                                   const float4 &z) {
  __m128 x0x2y0y2 = _mm_shuffle_ps(x, y, _MM_SHUFFLE(2, 0, 2, 0));
  __m128 y1y3z1z3 = _mm_shuffle_ps(y, z, _MM_SHUFFLE(3, 1, 3, 1));
  __m128 z0z2x1x3 = _mm_shuffle_ps(z, x, _MM_SHUFFLE(3, 1, 2, 0));

  __m128 rx0y0z0x1 =
      _mm_shuffle_ps(x0x2y0y2, z0z2x1x3, _MM_SHUFFLE(2, 0, 2, 0));
  __m128 ry1z1x2y2 =
      _mm_shuffle_ps(y1y3z1z3, x0x2y0y2, _MM_SHUFFLE(3, 1, 2, 0));
  __m128 rz2x3y3z3 =
      _mm_shuffle_ps(z0z2x1x3, y1y3z1z3, _MM_SHUFFLE(3, 1, 3, 1));

  _mm_store_ps(p + 0, rx0y0z0x1);
  _mm_store_ps(p + 4, ry1z1x2y2);
  _mm_store_ps(p + 8, rz2x3y3z3);
}

// example: mask=[1, 0, 0, 1], v=[1,2,3,4] -> ret = [1,4,x,x]
static inline float4 vgather_left_f4(const float4 &v, const int mask) {
#if defined(__SSE4_2__)

  __m128i maskTable[16] = {
      // [0, 0, 0, 0] -> (x,x,x,x)
      // [0, 0, 0, 1] -> (x,x,x,0)
      // [0, 0, 1, 0] -> (x,x,x,1)
      // [0, 0, 1, 1] -> (x,x,1,0)
      _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x80808080),
      _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x03020100),
      _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x07060504),
      _mm_set_epi32(0x80808080, 0x80808080, 0x07060504, 0x03020100),

      // [0, 1, 0, 0] -> (x,x,x,2)
      // [0, 1, 0, 1] -> (x,x,2,0)
      // [0, 1, 1, 0] -> (x,x,2,1)
      // [0, 1, 1, 1] -> (x,2,1,0)
      _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x0b0a0908),
      _mm_set_epi32(0x80808080, 0x80808080, 0x0b0a0908, 0x03020100),
      _mm_set_epi32(0x80808080, 0x80808080, 0x0b0a0908, 0x07060504),
      _mm_set_epi32(0x80808080, 0x0b0a0908, 0x07060504, 0x03020100),

      // [1, 0, 0, 0] -> (x,x,x,3)
      // [1, 0, 0, 1] -> (x,x,3,0)
      // [1, 0, 1, 0] -> (x,x,3,1)
      // [1, 0, 1, 1] -> (x,3,1,0)
      _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x0f0e0d0c),
      _mm_set_epi32(0x80808080, 0x80808080, 0x0f0e0d0c, 0x03020100),
      _mm_set_epi32(0x80808080, 0x80808080, 0x0f0e0d0c, 0x07060504),
      _mm_set_epi32(0x80808080, 0x0f0e0d0c, 0x07060504, 0x03020100),

      // [1, 1, 0, 0] -> (x,x,3,2)
      // [1, 1, 0, 1] -> (x,3,2,0)
      // [1, 1, 1, 0] -> (x,3,2,1)
      // [1, 1, 1, 1] -> (3,2,1,0)
      _mm_set_epi32(0x80808080, 0x80808080, 0x0f0e0d0c, 0x0b0a0908),
      _mm_set_epi32(0x80808080, 0x0f0e0d0c, 0x0b0a0908, 0x03020100),
      _mm_set_epi32(0x80808080, 0x0f0e0d0c, 0x0b0a0908, 0x07060504),
      _mm_set_epi32(0x0f0e0d0c, 0x0b0a0908, 0x07060504, 0x03020100),
  };

  return vcasttof4_i4(_mm_shuffle_epi8(vcasttoi4_f4(v), maskTable[mask]));
#else
  // @todo
  return vset1_f4(0.0f);
#endif
}

// example: mask=[1, 0, 0, 1], v=[1,2,3,4] -> ret = [x,x,1,4]
static inline float4 vgather_right_f4(const float4 &v, const int mask) {
#if defined(__SSE4_2__)
  __m128i maskTable[16] = {
      // [0, 0, 0, 0] -> (x,x,x,x)
      // [0, 0, 0, 1] -> (0,x,x,x)
      // [0, 0, 1, 0] -> (1,x,x,x)
      // [0, 0, 1, 1] -> (1,0,x,x)
      _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x80808080),
      _mm_set_epi32(0x03020100, 0x80808080, 0x80808080, 0x80808080),
      _mm_set_epi32(0x07060504, 0x80808080, 0x80808080, 0x80808080),
      _mm_set_epi32(0x07060504, 0x03020100, 0x80808080, 0x80808080),

      // [0, 1, 0, 0] -> (2,x,x,x)
      // [0, 1, 0, 1] -> (2,0,x,x)
      // [0, 1, 1, 0] -> (2,1,x,x)
      // [0, 1, 1, 1] -> (2,1,0,x)
      _mm_set_epi32(0x0b0a0908, 0x80808080, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0b0a0908, 0x03020100, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0b0a0908, 0x07060504, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0b0a0908, 0x07060504, 0x03020100, 0x80808080),

      // [1, 0, 0, 0] -> (3,x,x,x)
      // [1, 0, 0, 1] -> (3,0,x,x)
      // [1, 0, 1, 0] -> (3,1,x,x)
      // [1, 0, 1, 1] -> (3,1,0,x)
      _mm_set_epi32(0x0f0e0d0c, 0x80808080, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0f0e0d0c, 0x03020100, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0f0e0d0c, 0x07060504, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0f0e0d0c, 0x07060504, 0x03020100, 0x80808080),

      // [1, 1, 0, 0] -> (3,2,x,x)
      // [1, 1, 0, 1] -> (3,2,0,x)
      // [1, 1, 1, 0] -> (3,2,1,x)
      // [1, 1, 1, 1] -> (3,2,1,0)
      _mm_set_epi32(0x0f0e0d0c, 0x0b0a0908, 0x80808080, 0x80808080),
      _mm_set_epi32(0x0f0e0d0c, 0x0b0a0908, 0x03020100, 0x80808080),
      _mm_set_epi32(0x0f0e0d0c, 0x0b0a0908, 0x07060504, 0x80808080),
      _mm_set_epi32(0x0f0e0d0c, 0x0b0a0908, 0x07060504, 0x03020100),
  };

  return vcasttof4_i4(_mm_shuffle_epi8(vcasttoi4_f4(v), maskTable[mask]));
#else
  // @todo
  return vset1_f4(0.0f);
#endif
}

static inline int popcount_i4(unsigned int i) {
#if defined(__SSE4_2__)
  return _mm_popcnt_u32(i); // Assume i is calculated by movemask_ps(4bit mask)
#else
  return (i & 0x1) + ((i >> 1) & 0x1) + ((i >> 2) & 0x1) + ((i >> 3) & 0x1);
#endif
}

/*
 * Fast SSE2 implementation of special math functions.
 */

#define POLY0(x, c0) _mm_set1_ps(c0)
#define POLY1(x, c0, c1)                                                       \
  _mm_add_ps(_mm_mul_ps(POLY0(x, c1), x), _mm_set1_ps(c0))
#define POLY2(x, c0, c1, c2)                                                   \
  _mm_add_ps(_mm_mul_ps(POLY1(x, c1, c2), x), _mm_set1_ps(c0))
#define POLY3(x, c0, c1, c2, c3)                                               \
  _mm_add_ps(_mm_mul_ps(POLY2(x, c1, c2, c3), x), _mm_set1_ps(c0))
#define POLY4(x, c0, c1, c2, c3, c4)                                           \
  _mm_add_ps(_mm_mul_ps(POLY3(x, c1, c2, c3, c4), x), _mm_set1_ps(c0))
#define POLY5(x, c0, c1, c2, c3, c4, c5)                                       \
  _mm_add_ps(_mm_mul_ps(POLY4(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0))

#define EXP_POLY_DEGREE 3
#define LOG_POLY_DEGREE 5

/**
 * See http://www.devmaster.net/forums/showthread.php?p=43580
 */
static inline float4 vexp2_f4(float4 x) {
  __m128i ipart;
  __m128 fpart, expipart, expfpart;

  x = _mm_min_ps(x, _mm_set1_ps(129.00000f));
  x = _mm_max_ps(x, _mm_set1_ps(-126.99999f));

  /* ipart = int(x - 0.5) */
  ipart = _mm_cvtps_epi32(_mm_sub_ps(x, _mm_set1_ps(0.5f)));

  /* fpart = x - ipart */
  fpart = _mm_sub_ps(x, _mm_cvtepi32_ps(ipart));

  /* expipart = (float) (1 << ipart) */
  expipart = _mm_castsi128_ps(
      _mm_slli_epi32(_mm_add_epi32(ipart, _mm_set1_epi32(127)), 23));

/* minimax polynomial fit of 2**x, in range [-0.5, 0.5[ */
#if EXP_POLY_DEGREE == 5
  expfpart = POLY5(fpart, 9.9999994e-1f, 6.9315308e-1f, 2.4015361e-1f,
                   5.5826318e-2f, 8.9893397e-3f, 1.8775767e-3f);
#elif EXP_POLY_DEGREE == 4
  expfpart = POLY4(fpart, 1.0000026f, 6.9300383e-1f, 2.4144275e-1f,
                   5.2011464e-2f, 1.3534167e-2f);
#elif EXP_POLY_DEGREE == 3
  expfpart =
      POLY3(fpart, 9.9992520e-1f, 6.9583356e-1f, 2.2606716e-1f, 7.8024521e-2f);
#elif EXP_POLY_DEGREE == 2
  expfpart = POLY2(fpart, 1.0017247f, 6.5763628e-1f, 3.3718944e-1f);
#else
#error
#endif

  return _mm_mul_ps(expipart, expfpart);
}

/**
 * See http://www.devmaster.net/forums/showthread.php?p=43580
 */
static inline float4 vlog2_f4(float4 x) {
  __m128i expmask = _mm_set1_epi32(0x7f800000);
  __m128i mantmask = _mm_set1_epi32(0x007fffff);
  __m128 one = _mm_set1_ps(1.0f);

  __m128i i = _mm_castps_si128(x);

  /* exp = (float) exponent(x) */
  __m128 exp = _mm_cvtepi32_ps(_mm_sub_epi32(
      _mm_srli_epi32(_mm_and_si128(i, expmask), 23), _mm_set1_epi32(127)));

  /* mant = (float) mantissa(x) */
  __m128 mant = _mm_or_ps(_mm_castsi128_ps(_mm_and_si128(i, mantmask)), one);

  __m128 logmant;

/* Minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[
 * These coefficients can be generate with
 * http://www.boost.org/doc/libs/1_36_0/libs/math/doc/sf_and_dist/html/math_toolkit/toolkit/internals2/minimax.html
 */
#if LOG_POLY_DEGREE == 6
  logmant = POLY5(mant, 3.11578814719469302614f, -3.32419399085241980044f,
                  2.59883907202499966007f, -1.23152682416275988241f,
                  0.318212422185251071475f, -0.0344359067839062357313f);
#elif LOG_POLY_DEGREE == 5
  logmant = POLY4(mant, 2.8882704548164776201f, -2.52074962577807006663f,
                  1.48116647521213171641f, -0.465725644288844778798f,
                  0.0596515482674574969533f);
#elif LOG_POLY_DEGREE == 4
  logmant = POLY3(mant, 2.61761038894603480148f, -1.75647175389045657003f,
                  0.688243882994381274313f, -0.107254423828329604454f);
#elif LOG_POLY_DEGREE == 3
  logmant = POLY2(mant, 2.28330284476918490682f, -1.04913055217340124191f,
                  0.204446009836232697516f);
#else
#error
#endif

  /* This effectively increases the polynomial degree by one, but ensures that
   * log2(1) == 0*/
  logmant = _mm_mul_ps(logmant, _mm_sub_ps(mant, one));

  return _mm_add_ps(logmant, exp);
}

static inline float4 vpow_f4(float4 x, float4 y) {
  return vexp2_f4(_mm_mul_ps(vlog2_f4(x), y));
}

// float8 operations
#if defined(__AVX__)
#define vadd_f8(a, b) _mm256_add_ps((a), (b))
#define vsub_f8(a, b) _mm256_sub_ps((a), (b))
#define vmul_f8(a, b) _mm256_mul_ps((a), (b))
#define vcmplt_f8(a, b) _mm256_cmplt_ps((a), (b))
#define vcmple_f8(a, b) _mm256_cmple_ps((a), (b))
#define vcmpgt_f8(a, b) _mm256_cmpgt_ps((a), (b))
#define vcmpge_f8(a, b) _mm256_cmpge_ps((a), (b))
#define vmax_f8(a, b) _mm256_max_ps((a), (b))
#define vmin_f8(a, b) _mm256_min_ps((a), (b))
#define vand_f8(a, b) _mm256_and_ps((a), (b))
#define vor_f8(a, b) _mm256_or_ps((a), (b))

#define vdot_f8(ax, ay, az, bx, by, bz)                                        \
  (_mm256_add_ps(                                                              \
      _mm256_mul_ps((ax), (bx)),                                               \
      _mm256_add_ps(_mm256_mul_ps((ay), (by)), _mm256_mul_ps((az), (bz)))))

#define vcross_f8(a, b, c, d)                                                  \
  (_mm256_sub_ps(_mm256_mul_ps((a), (c)), _mm256_mul_ps((b), (d))))

#define vsel_f8(a, b, m) _mm256_blendv_ps((a), (b), (m))

#define vset1_f8(a) _mm256_set1_ps((a))
#define vset4_f8(a, b, c, d, e, f, g, h)                                       \
  _mm256_set_ps((a), (b), (c), (d), (e), (f), (g), (h))
#define vsetr4_f8(a, b, c, d, e, f, g, h)                                      \
  _mm256_setr_ps((a), (b), (c), (d), (e), (f), (g), (h))
#define vzero_f8() _mm256_setzero_ps()
#define vone_f8() _mm256_set1_ps(1.0f)
#define vshuffle_f8(a, b, i) _mm256_shuffle_ps((a), (b), (i))
#define vloadu_f8(p) _mm256_loadu_ps((p))
#define vstoreu_f8(p, v) _mm256_storeu_ps((p), (v))

#if defined(_MSC_VER)
#define vcasttoi8_f8(v) _mm256_castps_si256((v))
#define vcasttof8_i8(v) _mm256_castsi128_ps((v))
#else
#define vcasttoi8_f8(v) (__m256i)((v))
#define vcasttof8_i8(v) (__m256)((v))
#endif

// [xyzxyzxyzxyzxyzxyzxyzxyz] -> [xxxxxxxx, yyyyyyyy, zzzzzzz]
static inline void vaos2soa_xyz_f8(float8 &x, float8 &y, float8 &z,
                                   const float *p) {
  // See details at
  // http://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx/
  __m128 *m = (__m128 *)p;
  __m256 m03;
  __m256 m14;
  __m256 m25;
  m03 = _mm256_castps128_ps256(m[0]); // load lower halves
  m14 = _mm256_castps128_ps256(m[1]);
  m25 = _mm256_castps128_ps256(m[2]);
  m03 = _mm256_insertf128_ps(m03, m[3], 1); // load upper halves
  m14 = _mm256_insertf128_ps(m14, m[4], 1);
  m25 = _mm256_insertf128_ps(m25, m[5], 1);

  __m256 xy =
      _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE(2, 1, 3, 2)); // upper x's and y's
  __m256 yz =
      _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE(1, 0, 2, 1)); // lower y's and z's
  x = _mm256_shuffle_ps(m03, xy, _MM_SHUFFLE(2, 0, 3, 0));
  y = _mm256_shuffle_ps(yz, xy, _MM_SHUFFLE(3, 1, 2, 0));
  z = _mm256_shuffle_ps(yz, m25, _MM_SHUFFLE(3, 0, 3, 1));
}

// [xxxxxxxx, yyyyyyyy, zzzzzzz] -> [xyzxyzxyzxyzxyzxyzxyzxyz]
static inline void vsoa2aos_xyz_f8(float *p, const float8 &x, const float8 &y,
                                   const float8 &z) {
  // See details at
  // http://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx/
  __m128 *m = (__m128 *)p;

  __m256 rxy = _mm256_shuffle_ps(x, y, _MM_SHUFFLE(2, 0, 2, 0));
  __m256 ryz = _mm256_shuffle_ps(y, z, _MM_SHUFFLE(3, 1, 3, 1));
  __m256 rzx = _mm256_shuffle_ps(z, x, _MM_SHUFFLE(3, 1, 2, 0));

  __m256 r03 = _mm256_shuffle_ps(rxy, rzx, _MM_SHUFFLE(2, 0, 2, 0));
  __m256 r14 = _mm256_shuffle_ps(ryz, rxy, _MM_SHUFFLE(3, 1, 2, 0));
  __m256 r25 = _mm256_shuffle_ps(rzx, ryz, _MM_SHUFFLE(3, 1, 3, 1));

  m[0] = _mm256_castps256_ps128(r03);
  m[1] = _mm256_castps256_ps128(r14);
  m[2] = _mm256_castps256_ps128(r25);
  m[3] = _mm256_extractf128_ps(r03, 1);
  m[4] = _mm256_extractf128_ps(r14, 1);
  m[5] = _mm256_extractf128_ps(r25, 1);
}

#endif // __AVX__

// Integer operations

#define vset1_i4(a) _mm_set1_epi32((a))
#define vset_i4(a, b, c, d) _mm_set_epi32((a), (b), (c), (d))
#define vsetr_i4(a, b, c, d) _mm_setr_epi32((a), (b), (c), (d))
#define vloadu_i4(p) _mm_loadu_si128((__m128i *)(p))
#define vstoreu_i4(p, a) _mm_storeu_si128((__m128i *)(p), (a))
#define vsel_i4(a, b, m)                                                       \
  _mm_or_si128(_mm_and_si128((m), (a)), _mm_andnot_si128((m), (b)))
#define vcmpneq_i4(a, b)                                                       \
  _mm_andnot_si128(_mm_cmpeq_epi32((a), (b)), _mm_set1_epi32(0xFFFFFFFF))

#else // non SSE

static float4 vset1_f4(float a) {
  float4 c;
  c.v[0] = a;
  c.v[1] = a;
  c.v[2] = a;
  c.v[3] = a;
  return c;
}

static float4 vzero_f4() {
  float4 c;
  c.v[0] = 0.0f;
  c.v[1] = 0.0f;
  c.v[2] = 0.0f;
  c.v[3] = 0.0f;
  return c;
}

static float4 vone_f4() {
  float4 c;
  c.v[0] = 1.0f;
  c.v[1] = 1.0f;
  c.v[2] = 1.0f;
  c.v[3] = 1.0f;
  return c;
}

static float4 vloadu_f4(const float *p) {
  float4 c;
  c.v[0] = p[0];
  c.v[1] = p[1];
  c.v[2] = p[2];
  c.v[3] = p[3];
  return c;
}

static void vstoreu_f4(float *p, const float4 &v) {
  p[0] = v.v[0];
  p[1] = v.v[1];
  p[2] = v.v[2];
  p[3] = v.v[3];
}

static float4 vadd_f4(float4 a, float4 b) {
  float4 c;
  c.v[0] = a.v[0] + b.v[0];
  c.v[1] = a.v[1] + b.v[1];
  c.v[2] = a.v[2] + b.v[2];
  c.v[3] = a.v[3] + b.v[3];
  return c;
}

static float4 vsub_f4(float4 a, float4 b) {
  float4 c;
  c.v[0] = a.v[0] - b.v[0];
  c.v[1] = a.v[1] - b.v[1];
  c.v[2] = a.v[2] - b.v[2];
  c.v[3] = a.v[3] - b.v[3];
  return c;
}

static float4 vmul_f4(float4 a, float4 b) {
  float4 c;
  c.v[0] = a.v[0] * b.v[0];
  c.v[1] = a.v[1] * b.v[1];
  c.v[2] = a.v[2] * b.v[2];
  c.v[3] = a.v[3] * b.v[3];
  return c;
}

//#define vcmplt_f4(a, b) vcltq_f32((a), (b))
//#define vcmple_f4(a, b) vcleq_f32((a), (b))
//#define vcmpgt_f4(a, b) vcgtq_f32((a), (b))
//#define vcmpge_f4(a, b) vcgeq_f32((a), (b))
//#define vmax_f4(a, b) vmaxq_f32((a), (b))
//#define vmin_f4(a, b) vminq_f32((a), (b))
//#define vand_f4(a, b) vandq_f32((a), (b))
//#define vor_f4(a, b) vorq_f32((a), (b))

static float4 vdot_f4(float4 ax, float4 ay, float4 az, float4 bx, float4 by,
                      float4 bz) {
  float4 c;
  c.v[0] = ax.v[0] * bx.v[0] + ay.v[0] * by.v[0] + az.v[0] * bz.v[0];
  c.v[1] = ax.v[1] * bx.v[1] + ay.v[1] * by.v[1] + az.v[1] * bz.v[1];
  c.v[2] = ax.v[2] * bx.v[2] + ay.v[2] * by.v[2] + az.v[2] * bz.v[2];
  c.v[3] = ax.v[3] * bx.v[3] + ay.v[3] * by.v[3] + az.v[3] * bz.v[3];
  return c;
}

static float4 vcross_f4(float4 a, float4 b, float4 c, float4 d) {
  float4 e;
  e.v[0] = a.v[0] * c.v[0] - b.v[0] * d.v[0];
  e.v[1] = a.v[1] * c.v[1] - b.v[1] * d.v[1];
  e.v[2] = a.v[2] * c.v[2] - b.v[2] * d.v[2];
  e.v[3] = a.v[3] * c.v[3] - b.v[3] * d.v[3];
  return e;
}

static float4 vabs_f4(float4 a) {
  float4 c;
  c.v[0] = fabsf(a.v[0]);
  c.v[1] = fabsf(a.v[1]);
  c.v[2] = fabsf(a.v[2]);
  c.v[3] = fabsf(a.v[3]);
  return c;
}

static inline float4 vfastinv_f4(const float4 &v) {
  // @todo { Use fast inv inst }
  float4 inv;

  inv.v[0] = 1.0f / v.v[0];
  inv.v[1] = 1.0f / v.v[1];
  inv.v[2] = 1.0f / v.v[2];
  inv.v[3] = 1.0f / v.v[3];

  return inv;
}

static inline float4 vfastrsqrt_f4(const float4 &v) {
  // @todo { Use fast sqrt inst }
  float4 inv;

  inv.v[0] = 1.0f / sqrt(v.v[0]);
  inv.v[1] = 1.0f / sqrt(v.v[1]);
  inv.v[2] = 1.0f / sqrt(v.v[2]);
  inv.v[3] = 1.0f / sqrt(v.v[3]);

  return inv;
}

#endif // !SSE

#endif // __LSGL_SIMD_H__
