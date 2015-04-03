/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_CAMERA_HPP__
#define __LSGL_CAMERA_HPP__

#include "render_common.h"

namespace lsgl {
namespace render {

/// Camera class
class Camera {
public:
  Camera(real3 eye, real3 lookat, real3 up, float fov, bool isOrtho = false) {
    m_eye = eye;
    m_lookat = lookat;
    m_up = up;
    m_fov = fov;
    m_isOrtho = isOrtho;
    m_isStereoEnv = false;

    real3 look = m_lookat - m_eye;
    real3 u = cross(m_up, look);
    u.normalize();

    real3 v = cross(look, u);
    v.normalize();

    look.normalize();

    m_uvw[0] = u;
    m_uvw[1] = v;
    m_uvw[2] = look;
  }

  // Create stereo environment camera
  Camera(real3 eye, real3 lookat, real3 up, float zeroParallax,
         float eyeSeparation) {
    m_eye = eye;
    m_lookat = lookat;
    m_up = up;
    m_fov = zeroParallax; // fixme
    m_isOrtho = false;
    m_isStereoEnv = true;
    m_zeroParallax = zeroParallax;
    m_eyeSeparation = eyeSeparation;

    real3 look = m_lookat - m_eye;
    real3 u = cross(m_up, look);
    u.normalize();

    real3 v = cross(look, u);
    v.normalize();

    look.normalize();

    m_uvw[0] = u;
    m_uvw[1] = v;
    m_uvw[2] = look;
  }

  ~Camera() {}

  /// Build camera frame useful for setting up raytracing.
  void buildFrame(real3 &corner, real3 &u, real3 &v, int imageWidth,
                  int imageHeight);

  inline const real3 &getEye() const { return m_eye; }
  inline const real3 &getLookat() const { return m_lookat; }
  inline const real3 &getUp() const { return m_up; }
  inline float getFov() const { return m_fov; }
  inline bool isOrtho() const { return m_isOrtho; }
  inline bool isStereoEnv() const { return m_isStereoEnv; }

  /// Env & stereo camera(longitude-latitude coordinate)
  void generateStereoEnvRay(real3 &org, real3 &dir, float imageU, float imageV,
                            int imageWidth, int imageHeight);

private:
  real3 m_eye;
  real3 m_lookat;
  real3 m_up;
  float m_fov;
  float m_zeroParallax;
  float m_eyeSeparation;
  bool m_isOrtho;
  bool m_isStereoEnv;

  real3 m_uvw[3];
};

} // namespace render
} // namespace lsgl

#endif // __LSGL_CONTEXT_HPP__
