/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_RENDER_TEXTURE_HPP__
#define __LSGL_RENDER_TEXTURE_HPP__

extern "C" {

// 3D Volume texture is defined in "C" linkage for optimization.
#define LSGL_RENDER_TEXTURE3D_FORMAT_INVALID (0)
#define LSGL_RENDER_TEXTURE3D_FORMAT_BYTE (1 << 0)
#define LSGL_RENDER_TEXTURE3D_FORMAT_USHORT (1 << 1)
#define LSGL_RENDER_TEXTURE3D_FORMAT_UINT (1 << 2)
#define LSGL_RENDER_TEXTURE3D_FORMAT_FLOAT (1 << 3)
#define LSGL_RENDER_TEXTURE3D_FORMAT_DOUBLE (1 << 4)

#ifdef __GNUC__
#pragma pack(push, 1)
#define ALIGNMENT __attribute__((packed))
#else
#pragma pack(1)
#define ALIGNMENT
#endif // __GNUC__
typedef struct _Texture3D {
  unsigned char *image; // Needs reinterpret cast based on 'type' parameter
  int width;
  int height;
  int depth;
  int components;
  int data_type; // one of LSGL_RENDER_TEXTURE3D_FORMAT_***
} Texture3D;
#ifdef __GNUC__
#pragma pack(pop)
#else
#pragma pack()
#endif // __GNUC__

// Trilinear textel fetch for byte format 3D texture.
extern void FilterTexture3DByte(float *rgba, const Texture3D *tex, float u,
                                float v, float r);

// Trilinear textel fetch for float format 3D texture.
extern void FilterTexture3DFloat(float *rgba, const Texture3D *tex, float u,
                                 float v, float r);

// Trilinear textel fetch for double format 3D texture.
extern void FilterTexture3DDouble(float *rgba, const Texture3D *tex, float u,
                                  float v, float r);
}

namespace lsgl {
namespace render {

/// Image texture class
class Texture2D {

public:
  typedef enum { FORMAT_BYTE, FORMAT_FLOAT32, FORMAT_FLOAT64, } Format;

  Texture2D(const unsigned char *image, int width, int height, int components,
            Format format) {
    m_width = width;
    m_height = height;
    m_image = image;
    m_invWidth = 1.0f / width;
    m_invHeight = 1.0f / height;
    m_components = components;
    m_format = format;
  }

  ~Texture2D() {}

  int width() const { return m_width; }

  int height() const { return m_height; }

  int components() const { return m_components; }

  Format format() const { return m_format; }

  const unsigned char *image() const { return m_image; }

  // Bilinear textel fetch.
  void fetch(float *rgba, float u, float v) const;

  // Fetch filtered (i, j), (i+1, j), (i, j+1) texel.
  void fetchD(float *rgba0, float *rgba1, float *rgba2, float u, float v) const;

private:
  int m_width;
  int m_height;
  float m_invWidth;
  float m_invHeight;
  int m_components;
  const unsigned char *m_image;
  Format m_format;
};

#if 0
/// 3D volume texture class.
class Texture3D {
public:
  Texture3D(const float *image, int width, int height, int depth,
            int components) {
    m_width = width;
    m_height = height;
    m_depth = depth;
    m_image = image;
    m_components = components;
  }

  ~Texture3D() {
    // Texture image is maintained in application, so don't free it here.
  }

  int width() const { return m_width; }

  int height() const { return m_height; }

  int depth() const { return m_depth; }

  int components() const { return m_components; }

  const float *image() const { return m_image; }

  // Trilinear textel fetch.
  void fetch(float *rgba, float u, float v, float r) const;

private:

  int m_width;
  int m_height;
  int m_depth;
  int m_components;
  const float *m_image;

};
#endif

} // render
} // lsgl

#endif // __LSGL_RENER_TEXTURE_HPP__
