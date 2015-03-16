#ifndef __GLSL_MATH_FUNC_H__
#define __GLSL_MATH_FUNC_H__

#include <math.h>
#include <assert.h>

#include "glsl_simdclass.h"

const float PI = 3.141592;

inline float __glsl_radians(float degrees)
{
    return (degrees * PI) / 180;
}

inline vec2 __glsl_radians(vec2 degrees)
{
    vec2 t_vec2;
    t_vec2.v[0] = degrees.v[0] * PI / 180;
    t_vec2.v[1] = degrees.v[1] * PI / 180;
    return t_vec2;
}

inline vec3 __glsl_radians(vec3 degrees)
{
    vec3 t_vec3;
    t_vec3.v[0] = degrees.v[0] * PI / 180;
    t_vec3.v[1] = degrees.v[1] * PI / 180;
    t_vec3.v[2] = degrees.v[2] * PI / 180;
    return t_vec3;
}

inline vec4 __glsl_radians(vec4 degrees)
{
    vec4 t_vec4;
    t_vec4.v[0] = degrees.v[0] * PI / 180;
    t_vec4.v[1] = degrees.v[1] * PI / 180;
    t_vec4.v[2] = degrees.v[2] * PI / 180;
    t_vec4.v[3] = degrees.v[3] * PI / 180;
    return t_vec4;
}

inline float __glsl_degrees(float radians)
{
    return radians * 180 / PI;
}

inline vec2 __glsl_degrees(vec2 radians)
{
    vec2 t_vec2;
    t_vec2.v[0]=radians.v[0] * 180 / PI;
    t_vec2.v[1]=radians.v[1] * 180 / PI;
    return t_vec2;
}

inline vec3 __glsl_degrees(vec3 radians)
{
    vec3 t_vec3;
    t_vec3.v[0]=radians.v[0] * 180 / PI;
    t_vec3.v[1]=radians.v[1] * 180 / PI;
    t_vec3.v[2]=radians.v[2] * 180 / PI;
    return t_vec3;
}

inline vec4 __glsl_degrees(vec4 radians)
{
    vec4 t_vec4;
    t_vec4.v[0]=radians.v[0] * 180 / PI;
    t_vec4.v[1]=radians.v[1] * 180 / PI;
    t_vec4.v[2]=radians.v[2] * 180 / PI;
    t_vec4.v[3]=radians.v[3] * 180 / PI;
    return t_vec4;
}

inline float __glsl_sin(float angle)
{
    return sinf(angle);
}

inline vec2 __glsl_sin(vec2 angle)
{
    vec2 t_vec2;
    t_vec2.v[0]=sinf(angle.v[0]);
    t_vec2.v[1]=sinf(angle.v[1]);
    return t_vec2;
}

inline vec3 __glsl_sin(vec3 angle)
{
    vec3 t_vec3;
    t_vec3.v[0]=sinf(angle.v[0]);
    t_vec3.v[1]=sinf(angle.v[1]);
    t_vec3.v[2]=sinf(angle.v[2]);
    return t_vec3;
}

inline vec4 __glsl_sin(vec4 angle)
{
    vec4 t_vec4;
    t_vec4.v[0]=sinf(angle.v[0]);
    t_vec4.v[1]=sinf(angle.v[1]);
    t_vec4.v[2]=sinf(angle.v[2]);
    t_vec4.v[3]=sinf(angle.v[3]);
    return t_vec4;
}

inline float __glsl_cos(float angle)
{
    return cosf(angle);
}

inline vec2 __glsl_cos(vec2 angle)
{
    vec2 t_vec2;
    t_vec2.v[0]=cosf(angle.v[0]);
    t_vec2.v[1]=cosf(angle.v[1]);
    return t_vec2;
}

inline vec3 __glsl_cos(vec3 angle)
{
    vec3 t_vec3;
    t_vec3.v[0]=cosf(angle.v[0]);
    t_vec3.v[1]=cosf(angle.v[1]);
    t_vec3.v[2]=cosf(angle.v[2]);
    return t_vec3;
}

inline vec4 __glsl_cos(vec4 angle)
{
    vec4 t_vec4;
    t_vec4.v[0]=cosf(angle.v[0]);
    t_vec4.v[1]=cosf(angle.v[1]);
    t_vec4.v[2]=cosf(angle.v[2]);
    t_vec4.v[3]=cosf(angle.v[3]);
    return t_vec4;
}


inline float __glsl_tan(float angle)
{
    return tanf(angle);
}

inline vec2 __glsl_tan(vec2 angle)
{
    vec2 t_vec2;
    t_vec2.v[0]=tanf(angle.v[0]);
    t_vec2.v[1]=tanf(angle.v[1]);
    return t_vec2;
}

inline vec3 __glsl_tan(vec3 angle)
{
    vec3 t_vec3;
    t_vec3.v[0]=tanf(angle.v[0]);
    t_vec3.v[1]=tanf(angle.v[1]);
    t_vec3.v[2]=tanf(angle.v[2]);
    return t_vec3;
}

inline vec4 __glsl_tan(vec4 angle)
{
    vec4 t_vec4;
    t_vec4.v[0]=tanf(angle.v[0]);
    t_vec4.v[1]=tanf(angle.v[1]);
    t_vec4.v[2]=tanf(angle.v[2]);
    t_vec4.v[3]=tanf(angle.v[3]);
    return t_vec4;
}

inline float __glsl_asin(float x)
{
    return asinf(x);
}

inline vec2 __glsl_asin(vec2 angle)
{
    vec2 t_vec2;
    t_vec2.v[0]=asinf(angle.v[0]);
    t_vec2.v[1]=asinf(angle.v[1]);
    return t_vec2;
}

inline vec3 __glsl_asin(vec3 angle)
{
    vec3 t_vec3;
    t_vec3.v[0]=asinf(angle.v[0]);
    t_vec3.v[1]=asinf(angle.v[1]);
    t_vec3.v[2]=asinf(angle.v[2]);
    return t_vec3;
}

inline vec4 __glsl_asin(vec4 angle)
{
    vec4 t_vec4;
    t_vec4.v[0]=asinf(angle.v[0]);
    t_vec4.v[1]=asinf(angle.v[1]);
    t_vec4.v[2]=asinf(angle.v[2]);
    t_vec4.v[3]=asinf(angle.v[3]);
    return t_vec4;
}

inline float __glsl_acos(float x)
{
    return acosf(x);
}

inline vec2 __glsl_acos(vec2 angle)
{
    vec2 t_vec2;
    t_vec2.v[0]=acosf(angle.v[0]);
    t_vec2.v[1]=acosf(angle.v[1]);
    return t_vec2;
}

inline vec3 __glsl_acos(vec3 angle)
{
    vec3 t_vec3;
    t_vec3.v[0]=acosf(angle.v[0]);
    t_vec3.v[1]=acosf(angle.v[1]);
    t_vec3.v[2]=acosf(angle.v[2]);
    return t_vec3;
}

inline vec4 __glsl_acos(vec4 angle)
{
    vec4 t_vec4;
    t_vec4.v[0]=acosf(angle.v[0]);
    t_vec4.v[1]=acosf(angle.v[1]);
    t_vec4.v[2]=acosf(angle.v[2]);
    t_vec4.v[3]=acosf(angle.v[3]);
    return t_vec4;
}

inline float __glsl_atan(float x)
{
    return atanf(x);
}

inline vec2 __glsl_atan(vec2 angle)
{
    vec2 t_vec2;
    t_vec2.v[0]=atanf(angle.v[0]);
    t_vec2.v[1]=atanf(angle.v[1]);
    return t_vec2;
}

inline vec3 __glsl_atan(vec3 angle)
{
    vec3 t_vec3;
    t_vec3.v[0]=atanf(angle.v[0]);
    t_vec3.v[1]=atanf(angle.v[1]);
    t_vec3.v[2]=atanf(angle.v[2]);
    return t_vec3;
}

inline vec4 __glsl_atan(vec4 angle)
{
    vec4 t_vec4;
    t_vec4.v[0]=atanf(angle.v[0]);
    t_vec4.v[1]=atanf(angle.v[1]);
    t_vec4.v[2]=atanf(angle.v[2]);
    t_vec4.v[3]=atanf(angle.v[3]);
    return t_vec4;
}

inline float __glsl_atan(float y,float x)
{
    return atan2f(y,x);
}

inline vec2 __glsl_atan(vec2 x, vec2 y)
{
    vec2 t_vec2;
    t_vec2.v[0]=atan2f(x.v[0], y.v[0]);
    t_vec2.v[1]=atan2f(x.v[1], y.v[1]);
    return t_vec2;
}

inline vec3 __glsl_atan(vec3 x, vec3 y)
{
    vec3 t_vec3;
    t_vec3.v[0]=atan2f(x.v[0], y.v[0]);
    t_vec3.v[1]=atan2f(x.v[1], y.v[1]);
    t_vec3.v[2]=atan2f(x.v[2], y.v[2]);
    return t_vec3;
}

inline vec4 __glsl_atan(vec4 x, vec4 y)
{
    vec4 t_vec4;
    t_vec4.v[0]=atan2f(x.v[0], y.v[0]);
    t_vec4.v[1]=atan2f(x.v[1], y.v[1]);
    t_vec4.v[2]=atan2f(x.v[2], y.v[2]);
    t_vec4.v[3]=atan2f(x.v[3], x.v[3]);
    return t_vec4;
}

inline float __glsl_pow(float x, float y)
{
    return powf(x,y);
}

inline vec2 __glsl_pow(vec2 x, vec2 y)
{
    vec2 t_vec2;
    t_vec2.v[0]=powf(x.v[0], y.v[0]);
    t_vec2.v[1]=powf(x.v[1], y.v[1]);
    return t_vec2;
}

inline vec3 __glsl_pow(vec3 x, vec3 y)
{
    vec3 t_vec3;
    t_vec3.v[0]=powf(x.v[0], y.v[0]);
    t_vec3.v[1]=powf(x.v[1], y.v[1]);
    t_vec3.v[2]=powf(x.v[2], y.v[2]);
    return t_vec3;
}

inline vec4 __glsl_pow(vec4 x, vec4 y)
{
    vec4 t_vec4;
    t_vec4.v[0]=powf(x.v[0], y.v[0]);
    t_vec4.v[0]=powf(x.v[1], y.v[1]);
    t_vec4.v[0]=powf(x.v[2], y.v[2]);
    t_vec4.v[0]=powf(x.v[3], y.v[3]);
    return t_vec4;
}

inline float __glsl_exp(float x)
{
    return expf(x);
}

inline vec2 __glsl_exp(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=expf(x.v[0]);
    t_vec2.v[1]=expf(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_exp(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=expf(x.v[0]);
    t_vec3.v[1]=expf(x.v[1]);
    t_vec3.v[2]=expf(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_exp(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=expf(x.v[0]);
    t_vec4.v[1]=expf(x.v[1]);
    t_vec4.v[2]=expf(x.v[2]);
    t_vec4.v[3]=expf(x.v[3]);
    return t_vec4;
}

inline float __glsl_log(float x)
{
    return logf(x);
}

inline vec2 __glsl_log(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=logf(x.v[0]);
    t_vec2.v[1]=logf(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_log(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=logf(x.v[0]);
    t_vec3.v[1]=logf(x.v[1]);
    t_vec3.v[2]=logf(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_log(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=logf(x.v[0]);
    t_vec4.v[1]=logf(x.v[1]);
    t_vec4.v[2]=logf(x.v[2]);
    t_vec4.v[3]=logf(x.v[3]);
    return t_vec4;
}

inline float __glsl_exp2(float x)
{
    return powf(2.0f, x);
}

inline vec2 __glsl_exp2(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_exp2(x.v[0]);
    t_vec2.v[1]=__glsl_exp2(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_exp2(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_exp2(x.v[0]);
    t_vec3.v[1]=__glsl_exp2(x.v[1]);
    t_vec3.v[2]=__glsl_exp2(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_exp2(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_exp2(x.v[0]);
    t_vec4.v[1]=__glsl_exp2(x.v[1]);
    t_vec4.v[2]=__glsl_exp2(x.v[2]);
    t_vec4.v[3]=__glsl_exp2(x.v[3]);
    return t_vec4;
}

inline float __glsl_log2(float x)
{
    return logf(x) / logf(2.0f);
}

inline vec2 __glsl_log2(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_log2(x.v[0]);
    t_vec2.v[1]=__glsl_log2(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_log2(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_log2(x.v[0]);
    t_vec3.v[1]=__glsl_log2(x.v[1]);
    t_vec3.v[2]=__glsl_log2(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_log2(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_log2(x.v[0]);
    t_vec4.v[1]=__glsl_log2(x.v[1]);
    t_vec4.v[2]=__glsl_log2(x.v[2]);
    t_vec4.v[3]=__glsl_log2(x.v[3]);
    return t_vec4;
}

inline float __glsl_sqrt(float x)
{
    return sqrtf(x);
}

inline vec2 __glsl_sqrt(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=sqrtf(x.v[0]);
    t_vec2.v[1]=sqrtf(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_sqrt(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=sqrtf(x.v[0]);
    t_vec3.v[1]=sqrtf(x.v[1]);
    t_vec3.v[2]=sqrtf(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_sqrt(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=sqrtf(x.v[0]);
    t_vec4.v[1]=sqrtf(x.v[1]);
    t_vec4.v[2]=sqrtf(x.v[2]);
    t_vec4.v[3]=sqrtf(x.v[3]);
    return t_vec4;
}

inline float __glsl_inversesqrt(float x)
{
    return 1/sqrtf(x);
}

inline vec2 __glsl_inversesqrt(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=sqrtf(x.v[0]);
    t_vec2.v[1]=sqrtf(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_inversesqrt(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=sqrtf(x.v[0]);
    t_vec3.v[1]=sqrtf(x.v[1]);
    t_vec3.v[2]=sqrtf(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_inversesqrt_v4(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=sqrtf(x.v[0]);
    t_vec4.v[1]=sqrtf(x.v[1]);
    t_vec4.v[2]=sqrtf(x.v[2]);
    t_vec4.v[3]=sqrtf(x.v[3]);
    return t_vec4;
}

inline float __glsl_abs(float x)
{
    if(x>=0)return x;
    return -x;
}

inline vec2 __glsl_abs(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_abs(x.v[0]);
    t_vec2.v[1]=__glsl_abs(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_abs(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_abs(x.v[0]);
    t_vec3.v[1]=__glsl_abs(x.v[1]);
    t_vec3.v[2]=__glsl_abs(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_abs(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_abs(x.v[0]);
    t_vec4.v[1]=__glsl_abs(x.v[1]);
    t_vec4.v[2]=__glsl_abs(x.v[2]);
    t_vec4.v[3]=__glsl_abs(x.v[3]);
    return t_vec4;
}

inline float __glsl_sign(float x)
{
    if(x>0)
    {
        return 1.0;
    }
    else if(x==0)
    {
        return 0.0;
    }
    else
    {
        return -1.0;
    }
}

inline vec2 __glsl_sign(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_sign(x.v[0]);
    t_vec2.v[1]=__glsl_sign(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_sign(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_sign(x.v[0]);
    t_vec3.v[1]=__glsl_sign(x.v[1]);
    t_vec3.v[2]=__glsl_sign(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_sign(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_sign(x.v[0]);
    t_vec4.v[1]=__glsl_sign(x.v[1]);
    t_vec4.v[2]=__glsl_sign(x.v[2]);
    t_vec4.v[3]=__glsl_sign(x.v[3]);
    return t_vec4;
}

inline float __glsl_floor(float x)
{
    return floorf(x);
}

inline vec2 __glsl_floor(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_floor(x.v[0]);
    t_vec2.v[1]=__glsl_floor(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_floor(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_floor(x.v[0]);
    t_vec3.v[1]=__glsl_floor(x.v[1]);
    t_vec3.v[2]=__glsl_floor(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_floor(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_floor(x.v[0]);
    t_vec4.v[1]=__glsl_floor(x.v[1]);
    t_vec4.v[2]=__glsl_floor(x.v[2]);
    t_vec4.v[3]=__glsl_floor(x.v[3]);
    return t_vec4;
}

inline float __glsl_ceil(float x)
{
    return ceilf(x);
}

inline vec2 __glsl_ceil(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=ceilf(x.v[0]);
    t_vec2.v[1]=ceilf(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_ceil(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=ceilf(x.v[0]);
    t_vec3.v[1]=ceilf(x.v[1]);
    t_vec3.v[2]=ceilf(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_ceil(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=ceilf(x.v[0]);
    t_vec4.v[1]=ceilf(x.v[1]);
    t_vec4.v[2]=ceilf(x.v[2]);
    t_vec4.v[3]=ceilf(x.v[3]);
    return t_vec4;
}

inline float __glsl_fract(float x)
{
    return x-floorf(x);
}

inline vec2 __glsl_fract(vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=x.v[0]-floorf(x.v[0]);
    t_vec2.v[1]=x.v[1]-floorf(x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_fract(vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=x.v[0]-floorf(x.v[0]);
    t_vec3.v[1]=x.v[1]-floorf(x.v[1]);
    t_vec3.v[2]=x.v[2]-floorf(x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_fract(vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=x.v[0]-floorf(x.v[0]);
    t_vec4.v[1]=x.v[1]-floorf(x.v[1]);
    t_vec4.v[2]=x.v[2]-floorf(x.v[2]);
    t_vec4.v[3]=x.v[3]-floorf(x.v[3]);
    return t_vec4;
}

inline float __glsl_mod(float x,float y)
{
    return x-y*floor(x/y);
}

inline vec2 __glsl_mod(vec2 x, vec2 y)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_mod(x.v[0], y.v[0]);
    t_vec2.v[1]=__glsl_mod(x.v[1], y.v[1]);
    return t_vec2;
}

inline vec3 __glsl_mod(vec3 x, vec3 y)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_mod(x.v[0], y.v[0]);
    t_vec3.v[1]=__glsl_mod(x.v[1], y.v[1]);
    t_vec3.v[2]=__glsl_mod(x.v[2], y.v[2]);
    return t_vec3;
}

inline vec4 __glsl_mod(vec4 x, vec4 y)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_mod(x.v[0], y.v[0]);
    t_vec4.v[1]=__glsl_mod(x.v[1], y.v[1]);
    t_vec4.v[2]=__glsl_mod(x.v[2], y.v[2]);
    t_vec4.v[3]=__glsl_mod(x.v[3], y.v[3]);
    return t_vec4;
}

inline float __glsl_min(float x,float y)
{
    return y<x ? y : x;
}

inline vec2 __glsl_min(vec2 x, vec2 y)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_min(x.v[0], y.v[0]);
    t_vec2.v[1]=__glsl_min(x.v[1], y.v[1]);
    return t_vec2;
}

inline vec3 __glsl_min(vec3 x, vec3 y)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_min(x.v[0], y.v[0]);
    t_vec3.v[1]=__glsl_min(x.v[1], y.v[1]);
    t_vec3.v[2]=__glsl_min(x.v[2], y.v[2]);
    return t_vec3;
}

inline vec4 __glsl_min(vec4 x, vec4 y)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_min(x.v[0], y.v[0]);
    t_vec4.v[1]=__glsl_min(x.v[1], y.v[1]);
    t_vec4.v[2]=__glsl_min(x.v[2], y.v[2]);
    t_vec4.v[3]=__glsl_min(x.v[3], y.v[3]);
    return t_vec4;
}

inline float __glsl_max(float x,float y)
{
    return x>y ? x : y;
}

inline vec2 __glsl_max(vec2 x, vec2 y)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_max(x.v[0], y.v[0]);
    t_vec2.v[1]=__glsl_max(x.v[1], y.v[1]);
    return t_vec2;
}

inline vec3 __glsl_max(vec3 x, vec3 y)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_max(x.v[0], y.v[0]);
    t_vec3.v[1]=__glsl_max(x.v[1], y.v[1]);
    t_vec3.v[2]=__glsl_max(x.v[2], y.v[2]);
    return t_vec3;
}

inline vec4 __glsl_max(vec4 x, vec4 y)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_max(x.v[0], y.v[0]);
    t_vec4.v[1]=__glsl_max(x.v[1], y.v[1]);
    t_vec4.v[2]=__glsl_max(x.v[2], y.v[2]);
    t_vec4.v[3]=__glsl_max(x.v[3], y.v[3]);
    return t_vec4;
}

inline float __glsl_clamp(float x, float minVal, float maxVal)
{
    return x < minVal ? minVal : x > maxVal ? maxVal : x;
}

inline vec2 __glsl_clamp(vec2 x, vec2 minVal, vec2 maxVal)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_clamp(x.v[0], minVal.v[0], maxVal.v[0]);
    t_vec2.v[1]=__glsl_clamp(x.v[1], minVal.v[1], maxVal.v[1]);
    return t_vec2;
}

inline vec2 __glsl_clamp(vec2 x, float minVal, float maxVal)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_clamp(x.v[0], minVal, maxVal);
    t_vec2.v[1]=__glsl_clamp(x.v[1], minVal, maxVal);
    return t_vec2;
}

inline vec3 __glsl_clamp(vec3 x, vec3 minVal, vec3 maxVal)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_clamp(x.v[0], minVal.v[0], maxVal.v[0]);
    t_vec3.v[1]=__glsl_clamp(x.v[1], minVal.v[1], maxVal.v[1]);
    t_vec3.v[2]=__glsl_clamp(x.v[2], minVal.v[2], maxVal.v[2]);
    return t_vec3;
}

inline vec3 __glsl_clamp(vec3 x, float minVal, float maxVal)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_clamp(x.v[0], minVal, maxVal);
    t_vec3.v[1]=__glsl_clamp(x.v[1], minVal, maxVal);
    t_vec3.v[2]=__glsl_clamp(x.v[2], minVal, maxVal);
    return t_vec3;
}

inline vec4 __glsl_clamp(vec4 x, vec4 minVal, vec4 maxVal)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_clamp(x.v[0], minVal.v[0], maxVal.v[0]);
    t_vec4.v[1]=__glsl_clamp(x.v[1], minVal.v[1], maxVal.v[1]);
    t_vec4.v[2]=__glsl_clamp(x.v[2], minVal.v[2], maxVal.v[2]);
    t_vec4.v[3]=__glsl_clamp(x.v[3], minVal.v[3], maxVal.v[3]);
    return t_vec4;
}

inline vec4 __glsl_clamp(vec4 x, float minVal, float maxVal)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_clamp(x.v[0], minVal, maxVal);
    t_vec4.v[1]=__glsl_clamp(x.v[1], minVal, maxVal);
    t_vec4.v[2]=__glsl_clamp(x.v[2], minVal, maxVal);
    t_vec4.v[3]=__glsl_clamp(x.v[3], minVal, maxVal);
    return t_vec4;
}

inline float __glsl_mix(float x, float y, float a)
{
    return x*(1-a)+y*a;
}

inline vec2 __glsl_mix(vec2 x, vec2 y, vec2 a)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_mix(x.v[0], y.v[0], a.v[0]);
    t_vec2.v[1]=__glsl_mix(x.v[1], y.v[1], a.v[1]);
    return t_vec2;
}

inline vec2 __glsl_mix(vec2 x, vec2 y, float a)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_mix(x.v[0], y.v[0], a);
    t_vec2.v[1]=__glsl_mix(x.v[1], y.v[1], a);
    return t_vec2;
}

inline vec3 __glsl_mix(vec3 x, vec3 y, vec3 a)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_mix(x.v[0], y.v[0], a.v[0]);
    t_vec3.v[1]=__glsl_mix(x.v[1], y.v[1], a.v[1]);
    t_vec3.v[2]=__glsl_mix(x.v[2], y.v[2], a.v[2]);
    return t_vec3;
}

inline vec3 __glsl_mix(vec3 x, vec3 y, float a)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_mix(x.v[0], y.v[0], a);
    t_vec3.v[1]=__glsl_mix(x.v[1], y.v[1], a);
    t_vec3.v[2]=__glsl_mix(x.v[2], y.v[2], a);
    return t_vec3;
}

inline vec4 __glsl_mix(vec4 x, vec4 y, vec4 a)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_mix(x.v[0], y.v[0], a.v[0]);
    t_vec4.v[1]=__glsl_mix(x.v[1], y.v[1], a.v[1]);
    t_vec4.v[2]=__glsl_mix(x.v[2], y.v[2], a.v[2]);
    t_vec4.v[3]=__glsl_mix(x.v[3], y.v[3], a.v[3]);
    return t_vec4;
}

inline vec4 __glsl_mix(vec4 x, vec4 y, float a)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_mix(x.v[0], y.v[0], a);
    t_vec4.v[1]=__glsl_mix(x.v[1], y.v[1], a);
    t_vec4.v[2]=__glsl_mix(x.v[2], y.v[2], a);
    t_vec4.v[3]=__glsl_mix(x.v[3], y.v[3], a);
    return t_vec4;
}

inline float __glsl_step(float edge, float x)
{
    return x < edge ? 0.0 : 1.0;
}

inline vec2 __glsl_step(vec2 edge, vec2 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_step(edge.v[0], x.v[0]);
    t_vec2.v[1]=__glsl_step(edge.v[1], x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_step(vec3 edge, vec3 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_step(edge.v[0], x.v[0]);
    t_vec3.v[1]=__glsl_step(edge.v[1], x.v[1]);
    t_vec3.v[2]=__glsl_step(edge.v[2], x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_step(vec4 edge, vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_step(edge.v[0], x.v[0]);
    t_vec4.v[1]=__glsl_step(edge.v[1], x.v[1]);
    t_vec4.v[2]=__glsl_step(edge.v[2], x.v[2]);
    t_vec4.v[3]=__glsl_step(edge.v[3], x.v[3]);
    return t_vec4;
}

inline float __glsl_smoothstep(float edge0, float edge1, float x)
{
    float t,ct;
    ct = (x-edge0)/(edge1-edge0);  
    t = ct < 0.0 ? 0.0 : ct > 1.0 ? 1.0 : ct; // clamp to [0, 1]
    return t * t* (3.0-2.0*t);
}

inline vec2 __glsl_smoothstep(vec2 edge0, vec2 edge1, vec4 x)
{
    vec2 t_vec2;
    t_vec2.v[0]=__glsl_smoothstep(edge0.v[0], edge1.v[0], x.v[0]);
    t_vec2.v[1]=__glsl_smoothstep(edge0.v[1], edge1.v[1], x.v[1]);
    return t_vec2;
}

inline vec3 __glsl_smoothstep(vec3 edge0, vec3 edge1, vec4 x)
{
    vec3 t_vec3;
    t_vec3.v[0]=__glsl_smoothstep(edge0.v[0], edge1.v[0], x.v[0]);
    t_vec3.v[1]=__glsl_smoothstep(edge0.v[1], edge1.v[1], x.v[1]);
    t_vec3.v[2]=__glsl_smoothstep(edge0.v[2], edge1.v[2], x.v[2]);
    return t_vec3;
}

inline vec4 __glsl_smoothstep(vec4 edge0, vec4 edge1, vec4 x)
{
    vec4 t_vec4;
    t_vec4.v[0]=__glsl_smoothstep(edge0.v[0], edge1.v[0], x.v[0]);
    t_vec4.v[1]=__glsl_smoothstep(edge0.v[1], edge1.v[1], x.v[1]);
    t_vec4.v[2]=__glsl_smoothstep(edge0.v[2], edge1.v[2], x.v[2]);
    t_vec4.v[3]=__glsl_smoothstep(edge0.v[3], edge1.v[3], x.v[3]);
    return t_vec4;
}

inline float __glsl_length(float x)
{
    return x;
}

inline float __glsl_length(vec2 x)
{
    float len;    
    len=sqrtf(x.v[0] * x.v[0] + x.v[1] * x.v[1]);
    return len;
}

inline float __glsl_length(vec3 x)
{
    float len;    
    len=sqrtf(x.v[0]*x.v[0] + x.v[1]*x.v[1] + x.v[2]*x.v[2]);
    return len;
}

inline float __glsl_length(vec4 x)
{
    float len;    
    len=sqrtf(x.v[0]*x.v[0] + x.v[1]*x.v[1] + x.v[2]*x.v[2] + x.v[3]*x.v[3]);
    return len;
}

inline float __glsl_distance(float p0, float p1)
{
    return fabsf(p0-p1);
}

inline float __glsl_distance(vec2 p0, vec2 p1)
{
    float distance;    
    vec2 t_vec2;
    t_vec2.v[0] = p0.v[0] - p1.v[0];
    t_vec2.v[1] = p0.v[1] - p1.v[1];    
    distance=__glsl_length(t_vec2);
    return distance;
}

inline float __glsl_distance(vec3 p0, vec3 p1)
{
    float distance;
    vec3 t_vec3;
    t_vec3.v[0] = p0.v[0] - p1.v[0];
    t_vec3.v[1] = p0.v[1] - p1.v[1];
    t_vec3.v[2] = p0.v[2] - p1.v[2];
    distance=__glsl_length(t_vec3);
    return distance;
}

inline float __glsl_distance(vec4 p0, vec4 p1)
{
    float distance;
    vec4 t_vec4;
    t_vec4.v[0] = p0.v[0] - p1.v[0];
    t_vec4.v[1] = p0.v[1] - p1.v[1];
    t_vec4.v[2] = p0.v[2] - p1.v[2];
    t_vec4.v[3] = p0.v[3] - p1.v[3];    
    distance=__glsl_length(t_vec4);
    return distance;
}

inline float __glsl_dot(float p0, float p1)
{
    return p0 * p1;
}

inline float __glsl_dot(vec2 p0, vec2 p1)
{
    return p0.v[0] * p1.v[0] + p0.v[1] * p1.v[1];
}

inline float __glsl_dot(vec3 p0, vec3 p1)
{
    return p0.v[0] * p1.v[0] + p0.v[1] * p1.v[1] + p0.v[2] * p1.v[2];
}

inline float __glsl_dot(vec4 p0, vec4 p1)
{
    return p0.v[0] * p1.v[0] + p0.v[1] * p1.v[1] + p0.v[2] * p1.v[2] + p0.v[3] * p1.v[3];
}

inline vec3 __glsl_cross(vec3 x, vec3 y)
{
    vec3 t_vec3;
    t_vec3.v[0] = x.v[1] * y.v[2] - y.v[1] * x.v[2];
    t_vec3.v[1] = x.v[2] * y.v[0] - y.v[2] * x.v[0];
    t_vec3.v[2] = x.v[0] * y.v[1] - y.v[0] * x.v[1];
    return  t_vec3;
}

inline float __glsl_normalize_f(float x)
{
	return x;
	// x/length(x)
}

inline vec2 __glsl_normalize(vec2 x)
{
	float len;
	vec2 t_vec2;
	len=sqrtf(x.v[0] * x.v[0] + x.v[1] * x.v[1]);
	t_vec2.v[0] = x.v[0]/len;
	t_vec2.v[1] = x.v[1]/len;	
	return t_vec2;
}

inline vec3 __glsl_normalize(vec3 x)
{
	float len;
	vec3 t_vec3;
	len=sqrtf(x.v[0] * x.v[0] + x.v[1] * x.v[1] + x.v[2] * x.v[2]);
	t_vec3.v[0] = x.v[0]/len;
	t_vec3.v[1] = x.v[1]/len;
	t_vec3.v[2] = x.v[2]/len;	
	return t_vec3;
}

inline vec4 __glsl_normalize_v4(vec4 x)
{
	float len;
	vec4 t_vec4;
	len=sqrtf(x.v[0]*x.v[0] + x.v[1]*x.v[1] + x.v[2]*x.v[2] + x.v[3]*x.v[3]);
	t_vec4.v[0] = x.v[0]/len;
	t_vec4.v[1] = x.v[1]/len;
	t_vec4.v[2] = x.v[2]/len;
	t_vec4.v[3] = x.v[3]/len;	
	return t_vec4;
}

inline float __glsl_faceforward_f(float N, float I, float Nref)
{
	return Nref * I < 0 ? N : -N;
}

inline vec2 __glsl_faceforward(vec2 N, vec2 I, vec2 Nref)
{
	vec2 t_vec2;
	if(__glsl_dot(Nref,I)<0)
	{
		t_vec2.v[0]=N.v[0];
		t_vec2.v[1]=N.v[1];
	}
	else
	{
		t_vec2.v[0]=-N.v[0];
		t_vec2.v[1]=-N.v[1];
	}
	return t_vec2;
}

inline vec3 __glsl_faceforward(vec3 N, vec3 I, vec3 Nref)
{
	vec3 t_vec3;
	if(__glsl_dot(Nref,I)<0)
	{
		t_vec3.v[0]=N.v[0];
		t_vec3.v[1]=N.v[1];
		t_vec3.v[2]=N.v[2];
	}
	else
	{
		t_vec3.v[0]=-N.v[0];
		t_vec3.v[1]=-N.v[1];
		t_vec3.v[2]=-N.v[2];
	}
	return t_vec3;
}

inline vec4 __glsl_faceforward_v4(vec4 N, vec4 I, vec4 Nref)
{
	vec4 t_vec4;
	if(__glsl_dot(Nref,I)<0)
	{
		t_vec4.v[0]=N.v[0];
		t_vec4.v[1]=N.v[1];
		t_vec4.v[2]=N.v[2];
		t_vec4.v[3]=N.v[3];
	}
	else
	{
		t_vec4.v[0]=-N.v[0];
		t_vec4.v[1]=-N.v[1];
		t_vec4.v[2]=-N.v[2];
		t_vec4.v[3]=-N.v[3];
	}
	return t_vec4;
}

inline float __glsl_reflect_f(float I, float N)
{
	return I-2*(N*I)*N;
}

inline vec2 __glsl_reflect(vec2 I, vec2 N)
{
	float t;
	vec2 t_vec2;
	t=__glsl_dot(N,I);
	t_vec2.v[0] = I.v[0] - 2 * t * N.v[0];
	t_vec2.v[1] = I.v[1] - 2 * t * N.v[1];
	return t_vec2;
}

inline vec3 __glsl_reflect(vec3 I, vec3 N)
{
	float t;
	vec3 t_vec3;
	t=__glsl_dot(N,I);
	t_vec3.v[0] = I.v[0] - 2 * t * N.v[0];
	t_vec3.v[1] = I.v[1] - 2 * t * N.v[1];
	t_vec3.v[2] = I.v[2] - 2 * t * N.v[2];
	return t_vec3;
}

inline vec4 __glsl_reflect(vec4 I, vec4 N)
{
	float t;
	vec4 t_vec4;
	t=__glsl_dot(N,I);
	t_vec4.v[0] = I.v[0] - 2 * t * N.v[0];
	t_vec4.v[1] = I.v[1] - 2 * t * N.v[1];
	t_vec4.v[2] = I.v[2] - 2 * t * N.v[2];
	t_vec4.v[3] = I.v[3] - 2 * t * N.v[3];
	return t_vec4;
}

inline vec3 __glsl_refract(vec3 I, vec3 N, float eta)
{
  float k = 1.0f - eta * eta * (1.0f - __glsl_dot(N, I) * __glsl_dot(N, I));
	vec3 R;
  if (k < 0.0) {
    R.v[0] = 0.0f;
    R.v[1] = 0.0f;
    R.v[2] = 0.0f;
  } else {
    float sqK = sqrtf(k);
	  R.v[0] = eta * I.v[0] - (eta * N.v[0] * I.v[0] + sqK) * N.v[0];
	  R.v[1] = eta * I.v[1] - (eta * N.v[1] * I.v[1] + sqK) * N.v[1];
	  R.v[2] = eta * I.v[2] - (eta * N.v[2] * I.v[2] + sqK) * N.v[2];
  }

  return R;
}

inline mat2 __glsl_matrixCompMult(mat2 x, mat2 y)
{
    assert(0);
	mat2 t_mat2 = x;
	//t_mat2.v[0]1 = x.v[0]1 * y.v[0]1;
	//t_mat2.v[0]2 = x.v[0]2 * y.v[0]2;
	//t_mat2.v[1]1 = x.v[1]1 * y.v[1]1;
	//t_mat2.v[1]2 = x.v[1]2 * y.v[1]2;
	return t_mat2;
}

inline bvec2 __glsl_lessThan(vec2 x, vec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] < y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] < y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_lessThan(vec3 x, vec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] < y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] < y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] < y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_lessThan(vec4 x, vec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] < y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] < y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] < y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] < y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_lessThan(ivec2 x, ivec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] <= y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_lessThan(ivec3 x, ivec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] <= y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] <= y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_lessThan(ivec4 x, ivec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] <= y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] <= y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] <= y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_lessThanEqual(vec2 x, vec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] <= y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_lessThanEqual(vec3 x, vec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] <= y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] <= y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_lessThanEqual(vec4 x, vec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] <= y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] <= y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] <= y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_lessThanEqual(ivec2 x, ivec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] <= y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_lessThanEqual(ivec3 x, ivec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] <= y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] <= y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_lessThanEqual(ivec4 x, ivec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] <= y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] <= y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] <= y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] <= y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_greaterThan(vec2 x, vec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] > y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] > y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_greaterThan(vec3 x, vec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] > y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] > y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] > y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_greaterThan(vec4 x, vec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] > y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] > y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] > y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] > y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_greaterThan(ivec2 x, ivec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] > y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] > y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_greaterThan(ivec3 x, ivec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] > y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] > y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] > y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_greaterThan(ivec4 x, ivec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] > y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] > y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] > y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] > y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_greaterThanEqual(vec2 x, vec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] >= y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] >= y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_greaterThanEqual(vec3 x, vec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] >= y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] >= y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] >= y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_greaterThanEqual(vec4 x, vec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] >= y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] >= y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] >= y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] >= y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_greaterThanEqual(ivec2 x, ivec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] >= y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] >= y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_greaterThanEqual(ivec3 x, ivec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] >= y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] >= y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] >= y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_greaterThanEqual(ivec4 x, ivec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] >= y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] >= y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] >= y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] >= y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_equal(vec2 x, vec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] == y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_equal(vec3 x, vec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] == y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] == y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_equal(vec4 x, vec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] == y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] == y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] == y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_equal(ivec2 x, ivec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] == y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_equal(ivec3 x, ivec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] == y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] == y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_equal(ivec4 x, ivec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] == y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] == y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] == y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_equal(bvec2 x, bvec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] == y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_equal(bvec3 x, bvec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] == y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] == y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_equal(bvec4 x, bvec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] == y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] == y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] == y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] == y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_notEqual(vec2 x, vec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] != y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_notEqual(vec3 x, vec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] != y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] != y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_notEqual(vec4 x, vec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] != y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] != y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] != y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_notEqual(ivec2 x, ivec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] != y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_notEqual(ivec3 x, ivec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] != y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] != y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_notEqual(ivec4 x, ivec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] != y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] != y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] != y.v[3] ? true:false;
	return t_bvec4;
}

inline bvec2 __glsl_notEqual(bvec2 x, bvec2 y)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec2.v[1] = x.v[1] != y.v[1] ? true:false;
	return t_bvec2;
}

inline bvec3 __glsl_notEqual(bvec3 x, bvec3 y)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec3.v[1] = x.v[1] != y.v[1] ? true:false;
	t_bvec3.v[2] = x.v[2] != y.v[2] ? true:false;
	return t_bvec3;
}

inline bvec4 __glsl_notEqual(bvec4 x, bvec4 y)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] != y.v[0] ? true:false;
	t_bvec4.v[1] = x.v[1] != y.v[1] ? true:false;
	t_bvec4.v[2] = x.v[2] != y.v[2] ? true:false;
	t_bvec4.v[3] = x.v[3] != y.v[3] ? true:false;
	return t_bvec4;
}

inline bool __glsl_any(bvec2 x)
{
	return x.v[0] == true ? true : x.v[1] == true ? true : false;
}

inline bool __glsl_any(bvec3 x)
{
	return x.v[0] == true ? true : x.v[1] == true ? true : x.v[2] == true ? true : false;
}

inline bool __glsl_any(bvec4 x)
{
	return x.v[0] == true ? true : x.v[1] == true ? true : x.v[2] == true ? 
		true : x.v[3] == true ? true : false;
}

inline bool __glsl_all(bvec2 x)
{
	return x.v[0] == false ? false : x.v[1] == false ? false : true;
}

inline bool __glsl_all(bvec3 x)
{
	return x.v[0] == false ? false : x.v[1] == false ? false : x.v[2] == false ? 
		false : true;
}

inline bool __glsl_all(bvec4 x)
{
	return x.v[0] == false ? false : x.v[1] == false ? false : x.v[2] == false ? 
		false : x.v[3] == false ? false : true;
}


inline bvec2 __glsl_not(bvec2 x)
{
	bvec2 t_bvec2;
	t_bvec2.v[0] = x.v[0] == true ? false : true; 
	t_bvec2.v[1] = x.v[1] == true ? false : true; 
	return t_bvec2;
}

inline bvec3 __glsl_not(bvec3 x)
{
	bvec3 t_bvec3;
	t_bvec3.v[0] = x.v[0] == true ? false : true; 
	t_bvec3.v[0] = x.v[1] == true ? false : true; 
	t_bvec3.v[0] = x.v[2] == true ? false : true; 
	return t_bvec3;
}

inline bvec4 __glsl_not(bvec4 x)
{
	bvec4 t_bvec4;
	t_bvec4.v[0] = x.v[0] == true ? false : true; 
	t_bvec4.v[1] = x.v[1] == true ? false : true; 
	t_bvec4.v[2] = x.v[2] == true ? false : true; 
	t_bvec4.v[3] = x.v[3] == true ? false : true; 
	return t_bvec4;
}

#endif  // __GLSL_MATH_FUNC_H__
