/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

// @todo { Rewrite. }
#ifndef __LSGL_GLES_ACCEL_BUILDER_H__
#define __LSGL_GLES_ACCEL_BUILDER_H__

#include <time.h>
#include <vector>
#include <map>

#include "../render/render_prim_triangle.h"
#include "../render/render_accel_triangle.h"
#include "../render/render_accel_particle.h"
#include "../render/render_accel_line.h"
#include "../render/render_accel_tetra.h"
#include "../render/render_accel_solid.h"
#include "gles_common.h"
#include "GLES2/gl2.h"

namespace lsgl {

//
// Accel builder
//
class LSGLES_EXPORT AccelBuilder {
public:
  typedef ParticleAccel ParticleAccelerator;
  typedef LineAccel LineAccelerator;
  typedef TetraAccel TetraAccelerator;
  typedef SolidAccel SolidAccelerator;  // Pyramid, Prism, Hexa
  typedef TriangleAccel TriangleAccelerator;

  typedef enum {
    PRIMITIVE_INVALID,
    PRIMITIVE_TRIANGLES,
    PRIMITIVE_POINTS,
    PRIMITIVE_LINES,
    PRIMITIVE_TETRAHEDRONS,
    PRIMITIVE_PYRAMIDS,
    PRIMITIVE_PRISMS,
    PRIMITIVE_HEXAHEDRONS
  } PrimitiveType;

  struct ArrayBufInfo {
    GLuint offsetPosition;
    bool operator==(const ArrayBufInfo &info) const {
      return (offsetPosition == info.offsetPosition);
    }
  };

  struct PrimData {
    Mesh mesh;           // for triangle data
    Lines *lines;        // for line data
    Particles *points;   // for point data
    Tetrahedron *tetras; // for tetra data
    Solid *solids;       // for pyramid, prism, hexa data

    DataBuffer<GLuint> indexBuffer;
    DataBuffer<GLfloat> positionBuffer;
    DataBuffer<double> dpositionBuffer;
    DataBuffer<GLfloat> radiusBuffer; // for points
    unsigned char *accel;
    PrimitiveType type;

    inline PrimData()
        : accel(NULL), type(PRIMITIVE_INVALID), lines(NULL), points(NULL),
          tetras(NULL) {}

    inline ~PrimData() {
      // clear pointers inside mesh structure so it doesn't try to free our
      // internal data buffers in it's destructor
      mesh.faces = NULL;
      mesh.vertices = NULL;
      mesh.dvertices = NULL;

      assert(type != PRIMITIVE_INVALID);

      // For some reason, accel shold be already deleted in FreePrim(), not
      // here.
    }

    inline GLuint GetByteSize() {
      // @todo { remove }
      return indexBuffer.GetByteSize() + positionBuffer.GetByteSize() +
             dpositionBuffer.GetByteSize();
    }
  };

  AccelBuilder(GLuint maxCacheSize = (8192 * 1024));
  ~AccelBuilder();

  TriangleAccelerator *
  BuildTriangleAccel(const Buffer *elembuf, const Buffer *arraybuf,
                     bool isDoublePrecisionPos,
                     const VertexAttribute *vertexAttributes,
                     const GLuint *texture2D, GLsizei count, GLuint offset);

  ParticleAccelerator *BuildParticleAccel(
      const Buffer *elembuf, const Buffer *arraybuf, bool isDoublePrecisionPos,
      const VertexAttribute *vertexAttributes, GLsizei count, GLuint offset,
      GLfloat constantWidth, const GLfloat *widthBuf, GLsizei widthBufLen);
  LineAccelerator *BuildLineAccel(const Buffer *elembuf, const Buffer *arraybuf,
                                  bool isDoublePrecisionPos,
                                  const VertexAttribute *vertexAttributes,
                                  GLsizei count, GLuint offset,
                                  GLfloat constantWidth,
                                  const GLfloat *widthBuf, GLsizei widthBufLen,
                                  bool cap);
  TetraAccelerator *BuildTetraAccel(const Buffer *elembuf,
                                    const Buffer *arraybuf,
                                    bool isDoublePrecisionPos,
                                    const VertexAttribute *vertexAttributes,
                                    GLsizei count, GLuint offset);
  SolidAccelerator *BuildSolidAccel(GLenum solidType,
                                    const Buffer *elembuf,
                                    const Buffer *arraybuf,
                                    bool isDoublePrecisionPos,
                                    const VertexAttribute *vertexAttributes,
                                    GLsizei count, GLuint offset);
  static bool AddTriangleData(PrimData *md, TriangleAccelerator *accel);
  static bool AddParticleData(PrimData *md, ParticleAccelerator *accel);
  static bool AddTetraData(PrimData *md, TetraAccelerator *accel);
  static bool AddSolidData(PrimData *md, SolidAccelerator *accel);
  void Invalidate(const Buffer *buf);
  void EndFrame();

  inline GLuint GetCacheSize() const { return maxCacheSize_; }
  inline GLuint GetCacheUsed() const { return cacheSizeUsed_; }
  inline GLuint GetCacheFree() const { return GetCacheSize() - GetCacheUsed(); }
  inline float GetCachePercent() const {
    return ((float)GetCacheUsed() / (float)GetCacheSize()) * 100.0f;
  }

  static clock_t GetBuiltTime(TriangleAccelerator *accel);

private:
  struct CacheEntry {
    const Buffer *elembuf;
    const Buffer *arraybuf;
    bool isDoublePrecisionPos;
    ArrayBufInfo abinfo;
    GLsizei count;
    GLuint offset;
    GLuint hitCount;
    clock_t builtTime;
    bool dynamic;
    PrimData primData;
  };

private:
  template <typename T>
  static void AddData(DataBuffer<T> &dest, const Buffer *src, GLuint count,
                      GLuint offset = 0, GLuint adjust = 0,
                      T * maxValue = NULL);

  // template <typename T>
  // static void AddDataPad(bool enabled, const T *meshBuf, GLuint scalar,
  //                       GLuint numVertices, DataBuffer<T> &dest,
  //                       const Buffer *src, GLuint count, GLuint offset = 0,
  //                       GLuint adjust = 0, T *maxValue = NULL);

  static void AddTriangleData(PrimData *md, const Buffer *elembuf,
                              const Buffer *arraybuf, bool isDoublePrecisionPos,
                              const ArrayBufInfo *abinfo, GLsizei count,
                              GLuint offset);
  static void AddParticleData(PrimData *md, const Buffer *elembuf,
                              const Buffer *arraybuf, bool isDoublePrecisionPos,
                              const ArrayBufInfo *abinfo, GLsizei count,
                              GLuint offset, GLfloat constantWidth,
                              const GLfloat *widthBuf, GLsizei widthBufLen);
  static void AddLineData(PrimData *md, const Buffer *elembuf,
                          const Buffer *arraybuf, const ArrayBufInfo *abinfo,
                          GLsizei count, GLuint offset, GLfloat constantWidth,
                          const GLfloat *widthBuf, GLsizei widthBufLen,
                          bool cap);
  static void AddTetraData(PrimData *md, const Buffer *elembuf,
                           const Buffer *arraybuf, bool isDoublePrecisionPos,
                           const ArrayBufInfo *abinfo, GLsizei count,
                           GLuint offset);
  static void AddSolidData(GLenum solidType, PrimData *md, const Buffer *elembuf,
                           const Buffer *arraybuf, bool isDoublePrecisionPos,
                           const ArrayBufInfo *abinfo, GLsizei count,
                           GLuint offset);
  static void DebugDumpMesh(const Mesh *mesh);

  unsigned char *Locate(const Buffer *elembuf, const Buffer *arraybuf,
                        const ArrayBufInfo *abinfo, GLsizei count,
                        GLuint offset) const;
  void FreePrim(CacheEntry *ce);

  std::vector<CacheEntry *> staticList_, dynamicList_;

  typedef std::map<unsigned char *, CacheEntry *> PrimDataMap;
  static PrimDataMap primDataMap_;

  GLuint maxCacheSize_;
  GLuint cacheSizeUsed_;
};
}

#endif // __LSGL_GLES_ACCEL_BUILDER_H__
