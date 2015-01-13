#ifndef __SIMDCLASS_H__
#define __SIMDCLASS_H__

#if defined(__sparc__) && defined(__HPC_ACE__) // K/FX10

#include <emmintrin.h>

namespace {

#define FORCEINLINE __attribute__((always_inline))
#define POST_ALIGN(x)  __attribute__ ((aligned(x)))

static FORCEINLINE __m128d __d2_div(__m128d a, __m128d x)
{
    // Newton Raphson Reciprocal
    // (2 * Rcp(x)) - (x * Rcp(x) * Rcp(x))]
    // = (2 - (x * Rcp(x)) * Rcp(x)
    const __m128d two = _mm_set_pd(2.0, 2.0);
    const __m128d rcp = _fjsp_rcpa_v2r8(x);
    const __m128d rcp_nr0   = _fjsp_nmsub_v2r8(x, rcp, two);
    const __m128d rcp_nr1   = _mm_mul_pd(rcp_nr0, rcp);

    return _mm_mul_pd(a, rcp_nr1);
}

struct __vec2_d {

  __vec2_d() { }

  __vec2_d(double v0, double v1) {
    u.v0 = _mm_set_pd(v1, v0);
  }

  __vec2_d(__m128d v0) {
    u.v0 = v0;
  }

  __vec2_d(double v0) {
    u.v0 = _mm_set_pd(v0, v0);
  }

  __vec2_d(const __vec2_d &rhs) {
    u.v0 = rhs.u.v0;
  }

  __vec2_d operator=(const __vec2_d& rhs) {
    u.v0 = rhs.u.v0;
    return (*this);
  }

  __vec2_d operator+(const __vec2_d& rhs) const {
    return __vec2_d(_mm_add_pd(u.v0, rhs.u.v0));
  }

  __vec2_d operator-(const __vec2_d& rhs) const {
    return __vec2_d(_mm_sub_pd(u.v0, rhs.u.v0));
  }
  
  __vec2_d operator*(const __vec2_d& rhs) const {
    return __vec2_d(_mm_mul_pd(u.v0, rhs.u.v0));
  }

  __vec2_d operator/(const __vec2_d& rhs) const {
    return __vec2_d(__d2_div(u.v0, rhs.u.v0));
  }
  
  union {
    double v[2];
    struct { __m128d v0; };
  } u;
} POST_ALIGN(32);

inline __vec2_d operator*(double f, const __vec2_d& rhs) {
  return __vec2_d(f) * rhs;
}


inline __vec2_d vmin(const __vec2_d& a, const __vec2_d& b) {
  return __vec2_d(_mm_min_pd(a.u.v0, b.u.v0));
}

inline __vec2_d vmax(const __vec2_d& a, const __vec2_d& b) {
  return __vec2_d(_mm_max_pd(a.u.v0, b.u.v0));
}

struct __vec4_d {

  __vec4_d() { }

  inline __vec4_d(double v0, double v1, double v2, double v3) {
    u.v0 = _mm_set_pd(v1, v0);
    u.v1 = _mm_set_pd(v3, v2);
  }

  inline __vec4_d(const double* p) {
    u.v0 = _mm_load_pd(p);
    u.v1 = _mm_load_pd(p + 2);
  }

  inline __vec4_d(__m128d v0, __m128d v1) {
    u.v0 = v0;
    u.v1 = v1;
  }

  inline __vec4_d(double v0) {
    u.v0 = _mm_set_pd(v0, v0);
    u.v1 = _mm_set_pd(v0, v0);
  }

  inline __vec4_d(const __vec4_d &rhs) {
    u.v0 = rhs.u.v0;
    u.v1 = rhs.u.v1;
  }

  inline __vec4_d operator=(const __vec4_d& rhs) {
    u.v0 = rhs.u.v0;
    u.v1 = rhs.u.v1;
    return (*this);
  }

  inline __vec4_d operator+(const __vec4_d& rhs) const {
    return __vec4_d(_mm_add_pd(u.v0, rhs.u.v0),
                    _mm_add_pd(u.v1, rhs.u.v1));
  }

  inline __vec4_d operator-(const __vec4_d& rhs) const {
    return __vec4_d(_mm_sub_pd(u.v0, rhs.u.v0),
                    _mm_sub_pd(u.v1, rhs.u.v1));
  }
  
  inline __vec4_d operator*(const __vec4_d& rhs) const {
    return __vec4_d(_mm_mul_pd(u.v0, rhs.u.v0),
                    _mm_mul_pd(u.v1, rhs.u.v1));
  }

  inline __vec4_d operator/(const __vec4_d& rhs) const {
    return __vec4_d(__d2_div(u.v0, rhs.u.v0),
                    __d2_div(u.v1, rhs.u.v1));
  }
  
  union {
    double v[4];
    struct { __m128d v0, v1; };
  } u;
} POST_ALIGN(32);

inline __vec4_d operator*(double f, const __vec4_d& rhs) {
  return __vec4_d(f) * rhs;
}


inline __vec4_d vmin(const __vec4_d& a, const __vec4_d& b) {
  return __vec4_d(_mm_min_pd(a.u.v0, b.u.v0), _mm_min_pd(a.u.v1, b.u.v1));
}

inline __vec4_d vmax(const __vec4_d& a, const __vec4_d& b) {
  return __vec4_d(_mm_max_pd(a.u.v0, b.u.v0), _mm_max_pd(a.u.v1, b.u.v1));
}

}


#endif

#if defined(__SSE2__) || (_M_IX86_FP >= 2) // _M_IX86_FP 2 = VS /arch:SSE2

// SIMD class
#include <emmintrin.h>

namespace {

static inline __m128d rcpe(__m128d a)
{ return _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(a))); }


// 1/a with 53bits precision.
static inline __m128d rcp(__m128d a)
{
  __m128d b = rcpe(a);
  b = _mm_sub_pd(_mm_add_pd(b, b), _mm_mul_pd(_mm_mul_pd(a, b), b));
  b = _mm_sub_pd(_mm_add_pd(b, b), _mm_mul_pd(_mm_mul_pd(a, b), b));
  b = _mm_sub_pd(_mm_add_pd(b, b), _mm_mul_pd(_mm_mul_pd(a, b), b));
  return b;
}

// 1.0f / a with 24bits precision..
inline __m128 rcp(__m128 a)
{
  __m128 b  = _mm_rcp_ps(a);
  // (b+b) - a*b*b
  b = _mm_sub_ps(_mm_add_ps(b, b), _mm_mul_ps(_mm_mul_ps(a, b), b));
  b = _mm_sub_ps(_mm_add_ps(b, b), _mm_mul_ps(_mm_mul_ps(a, b), b));
  return b;
}

}

struct __vec4_f {

  __vec4_f() { }

  __vec4_f(float v0, float v1, float v2, float v3) {
    u.v0 = _mm_setr_ps(v0, v1, v2, v3);
  }

  __vec4_f(__m128 v0) {
    u.v0 = v0;
  }

  __vec4_f(float v0) {
    u.v0 = _mm_set1_ps(v0);
  }

  __vec4_f(const __vec4_f &rhs) {
    u.v0 = rhs.u.v0;
  }

  __vec4_f operator=(const __vec4_f& rhs) {
    u.v0 = rhs.u.v0;
    return (*this);
  }

  __vec4_f operator+(const __vec4_f& rhs) const {
    return __vec4_f(_mm_add_ps(u.v0, rhs.u.v0));
  }

  __vec4_f operator-(const __vec4_f& rhs) const {
    return __vec4_f(_mm_sub_ps(u.v0, rhs.u.v0));
  }
  
  __vec4_f operator*(const __vec4_f& rhs) const {
    return __vec4_f(_mm_mul_ps(u.v0, rhs.u.v0));
  }

  __vec4_f operator/(const __vec4_f& rhs) const {
    return __vec4_f(_mm_mul_ps(u.v0, rcp(rhs.u.v0)));
  }
  
  union {
    float v[4];
    struct { __m128 v0; };
  } u;
};

inline __vec4_f operator*(float f, const __vec4_f& rhs) {
  return __vec4_f(f) * rhs;
}

inline __vec4_f vmin(const __vec4_f& a, const __vec4_f& b) {
  return __vec4_f(_mm_min_ps(a.u.v0, b.u.v0));
}

inline __vec4_f vmax(const __vec4_f& a, const __vec4_f& b) {
  return __vec4_f(_mm_max_ps(a.u.v0, b.u.v0));
}

struct __vec2_d {

  __vec2_d() { }

  __vec2_d(double v0, double v1) {
    u.v0 = _mm_setr_pd(v0, v1);
  }

  __vec2_d(__m128d v0) {
    u.v0 = v0;
  }

  __vec2_d(double v0) {
    u.v0 = _mm_set_pd(v0, v0);
  }

  __vec2_d(const __vec2_d &rhs) {
    u.v0 = rhs.u.v0;
  }

  __vec2_d operator=(const __vec2_d& rhs) {
    u.v0 = rhs.u.v0;
    return (*this);
  }

  __vec2_d operator+(const __vec2_d& rhs) const {
    return __vec2_d(_mm_add_pd(u.v0, rhs.u.v0));
  }

  __vec2_d operator-(const __vec2_d& rhs) const {
    return __vec2_d(_mm_sub_pd(u.v0, rhs.u.v0));
  }
  
  __vec2_d operator*(const __vec2_d& rhs) const {
    return __vec2_d(_mm_mul_pd(u.v0, rhs.u.v0));
  }

  __vec2_d operator/(const __vec2_d& rhs) const {
    return __vec2_d(_mm_mul_pd(u.v0, rcp(rhs.u.v0)));
  }
  
  union {
    double v[2];
    struct { __m128d v0; };
  } u;
};

inline __vec2_d operator*(double f, const __vec2_d& rhs) {
  return __vec2_d(f) * rhs;
}


inline __vec2_d vmin(const __vec2_d& a, const __vec2_d& b) {
  return __vec2_d(_mm_min_pd(a.u.v0, b.u.v0));
}

inline __vec2_d vmax(const __vec2_d& a, const __vec2_d& b) {
  return __vec2_d(_mm_max_pd(a.u.v0, b.u.v0));
}

struct __vec4_d {

  __vec4_d() { }

  __vec4_d(double v0, double v1, double v2, double v3) {
    u.v0 = _mm_setr_pd(v0, v1);
    u.v1 = _mm_setr_pd(v2, v3);
  }

  __vec4_d(__m128d v0, __m128d v1) {
    u.v0 = v0;
    u.v1 = v1;
  }

  __vec4_d(double v0) {
    u.v0 = _mm_setr_pd(v0, v0);
    u.v1 = _mm_setr_pd(v0, v0);
  }

  __vec4_d(const __vec4_d &rhs) {
    u.v0 = rhs.u.v0;
    u.v1 = rhs.u.v1;
  }

  __vec4_d operator=(const __vec4_d& rhs) {
    u.v0 = rhs.u.v0;
    u.v1 = rhs.u.v1;
    return (*this);
  }

  __vec4_d operator+(const __vec4_d& rhs) const {
    return __vec4_d(_mm_add_pd(u.v0, rhs.u.v0),
                    _mm_add_pd(u.v1, rhs.u.v1));
  }

  __vec4_d operator-(const __vec4_d& rhs) const {
    return __vec4_d(_mm_sub_pd(u.v0, rhs.u.v0),
                    _mm_sub_pd(u.v1, rhs.u.v1));
  }
  
  __vec4_d operator*(const __vec4_d& rhs) const {
    return __vec4_d(_mm_mul_pd(u.v0, rhs.u.v0),
                    _mm_mul_pd(u.v1, rhs.u.v1));
  }

  __vec4_d operator/(const __vec4_d& rhs) const {
    return __vec4_d(_mm_mul_pd(u.v0, rcp(rhs.u.v0)),
                    _mm_mul_pd(u.v1, rcp(rhs.u.v1)));
  }
  
  union {
    double v[4];
    struct { __m128d v0, v1; };
  } u;
};

inline __vec4_d operator*(double f, const __vec4_d& rhs) {
  return __vec4_d(f) * rhs;
}

inline __vec4_d vmin(const __vec4_d& a, const __vec4_d& b) {
  return __vec4_d(_mm_min_pd(a.u.v0, b.u.v0), _mm_min_pd(a.u.v1, b.u.v1));
}

inline __vec4_d vmax(const __vec4_d& a, const __vec4_d& b) {
  return __vec4_d(_mm_max_pd(a.u.v0, b.u.v0), _mm_max_pd(a.u.v1, b.u.v1));
}

#endif

#endif  // __SIMDCLASS_H__
