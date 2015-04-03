/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
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
  typedef enum {
    FORMAT_BYTE,
    FORMAT_FLOAT32,
    FORMAT_FLOAT64,
  } Format;

  typedef enum {
    WRAP_REPEAT,
    WRAP_MIRRORED_REPEAT,
    WRAP_CLAMP_TO_EDGE,
    WRAP_CLAMP_TO_BORDER,
    WRAP_MIRROR_CLAMP_TO_BORDER,
  } WrapMode;

  Texture2D(const unsigned char *image, int width, int height, int components,
            Format format, WrapMode wrapMode = WRAP_REPEAT) {
    m_width = width;
    m_height = height;
    m_image = image;
    m_invWidth = 1.0f / width;
    m_invHeight = 1.0f / height;
    m_components = components;
    m_format = format;
    m_wrapMode[0] = wrapMode;
    m_wrapMode[1] = wrapMode;
    m_wrapMode[2] = wrapMode;
  }

  ~Texture2D() {}

  void setWrapS(WrapMode mode) { m_wrapMode[0] = mode; }

  void setWrapT(WrapMode mode) { m_wrapMode[1] = mode; }

  void setWrapR(WrapMode mode) { m_wrapMode[2] = mode; }

  int width() const { return m_width; }

  int height() const { return m_height; }

  int components() const { return m_components; }

  Format format() const { return m_format; }

  const unsigned char *image() const { return m_image; }

  /// Fetch texel color.
  void fetch(float *rgba, float u, float v, bool minFiltering,
             bool magFiltering) const;

  /// Fetch (i, j), (i+1, j), (i, j+1) texels.
  void fetchD(float *rgba0, float *rgba1, float *rgba2, float u, float v,
              bool minFiltering, bool magFiltering) const;

private:
  int m_width;
  int m_height;
  float m_invWidth;
  float m_invHeight;
  int m_components;
  const unsigned char *m_image;
  Format m_format;
  WrapMode m_wrapMode[3]; // S, T and R
};

} // render
} // lsgl

#endif // __LSGL_RENER_TEXTURE_HPP__
