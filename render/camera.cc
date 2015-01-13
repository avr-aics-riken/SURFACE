/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#define _USE_MATH_DEFINES // VS specific. Includes M_PI
#include <cmath>

#include "camera.h"

using namespace lsgl::render;

void Camera::buildFrame(vector3 &corner, vector3 &u, vector3 &v, int imageWidth,
                        int imageHeight) {
  if (isOrtho()) {

    vector3 look = m_lookat - m_eye;
    u = cross(look, m_up);
    u.normalize();

    v = cross(look, u);
    v.normalize();

    look.normalize();

    // Map to [-1, 1]^2
    u[0] /= (float)imageWidth;
    u[1] /= (float)imageWidth;
    u[2] /= (float)imageWidth;
    v[0] /= (float)imageHeight;
    v[1] /= (float)imageHeight;
    v[2] /= (float)imageHeight;
    corner = -(u * imageWidth + v * imageHeight);

    // @fixme { Is this correct calculation for ortho projection? }
    u = 2.0f * u;
    v = 2.0f * v;

  } else {
    // perspecive
    float flen =
        (0.5f * (float)imageHeight / tan(0.5 * (m_fov * M_PI / 180.0)));

    // Assume right-handed coordinate.
    vector3 look = m_lookat - m_eye;
    u = cross(look, m_up);
    u.normalize();

    v = cross(look, u);
    v.normalize();

    look.normalize();
    look = flen * look + m_eye;

    corner = look - 0.5 * (imageWidth * u + imageHeight * v);
  }

}

void Camera::generateStereoEnvRay(vector3& org, vector3 &dir, float imageU, float imageV,
                            int imageWidth, int imageHeight) {

  const bool isLeft = imageV < (imageHeight / 2);
  const float focalLength = m_zeroParallax;
  const float r = m_eyeSeparation;

  // [0, pi]
  float theta = M_PI * fmod(2.0 * imageV / imageHeight, 1.0); // [0, pi]

  // [0, 2pi]
  float phi = 2.0 * M_PI * (imageU / imageWidth);

  // +Y up
  vector3 dirO;
  dirO[0] = sinf(theta) * cosf(phi);
  dirO[1] = cosf(theta);
  dirO[2] = sinf(theta) * sinf(phi);

  vector3 parallax;
  if (isLeft) {
    // positive rotation
    parallax[0] = -dirO[2];
    parallax[1] = 0.0;
    parallax[2] = dirO[0];
  } else {
    // negative rotation
    parallax[0] = dirO[2];
    parallax[1] = 0.0;
    parallax[2] = -dirO[0];
  }

  parallax.normalize();
  parallax = parallax * r;

  org = m_eye + parallax;

  float psi = atan2(r, focalLength);

  if (isLeft) {
    psi = -psi; // negative rotataion
  }

  vector3 d;
  d[0] = dirO[0] * cos(psi) - dirO[2] * sin(psi);
  d[1] = dirO[1];
  d[2] = dirO[0] * sin(psi) + dirO[2] * cos(psi);


  // local -> world
  dir[0] = m_uvw[0][0] * d[0] + m_uvw[1][0] * d[1] + m_uvw[2][0] * d[2];
  dir[1] = m_uvw[0][1] * d[0] + m_uvw[1][1] * d[1] + m_uvw[2][1] * d[2];
  dir[2] = m_uvw[0][2] * d[0] + m_uvw[1][2] * d[1] + m_uvw[2][2] * d[2];

}
