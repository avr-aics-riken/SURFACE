/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
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
#include "gles_common.h"
#include "GLES2/gl2.h"

namespace lsgl {

//
// Mesh builder
//
class AccelBuilder {
public:
  typedef ParticleAccel ParticleAccelerator;
  typedef LineAccel LineAccelerator;
  typedef TetraAccel TetraAccelerator;
  typedef TriangleAccel MeshAccelerator;

  typedef enum {
    PRIMITIVE_INVALID,
    PRIMITIVE_TRIANGLES,
    PRIMITIVE_POINTS,
    PRIMITIVE_LINES,
    PRIMITIVE_TETRAHEDRONS
  } PrimitiveType;

  struct ArrayBufInfo {
    GLuint offsetPosition;
    bool operator==(const ArrayBufInfo &info) const {
      return (offsetPosition == info.offsetPosition);
    }
  };

  struct MeshData {
    // @toto { varying }
    Mesh mesh;
    DataBuffer<GLuint> indexBuffer;
    DataBuffer<GLfloat> positionBuffer;
    DataBuffer<double> dpositionBuffer;
    DataBuffer<GLfloat> radiusBuffer; // for points
    unsigned char *accel;
    PrimitiveType type;

    inline MeshData() : accel(NULL), type(PRIMITIVE_INVALID) {}

    inline ~MeshData() {
      // clear pointers inside mesh structure so it doesn't try to free our
      // internal data buffers in it's destructor
      mesh.faces = NULL;
      mesh.vertices = NULL;
      mesh.dvertices = NULL;

      assert(type != PRIMITIVE_INVALID);

      // For some reason, accel shold be already deleted in FreeMesh(), not
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

  MeshAccelerator *BuildMeshAccel(const Buffer *elembuf, const Buffer *arraybuf,
                                  bool isDoublePrecisionPos,
                                  const VertexAttribute *vertexAttributes,
                                  const GLuint *texture2D, GLsizei count,
                                  GLuint offset);

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
  static bool AddMeshData(MeshData *md, MeshAccelerator *accel);
  static bool AddParticleData(MeshData *md, ParticleAccelerator *accel);
  static bool AddTetraData(MeshData *md, TetraAccelerator *accel);
  void Invalidate(const Buffer *buf);
  void EndFrame();

  inline GLuint GetCacheSize() const { return maxCacheSize_; }
  inline GLuint GetCacheUsed() const { return cacheSizeUsed_; }
  inline GLuint GetCacheFree() const { return GetCacheSize() - GetCacheUsed(); }
  inline float GetCachePercent() const {
    return ((float)GetCacheUsed() / (float)GetCacheSize()) * 100.0f;
  }

  static clock_t GetBuiltTime(MeshAccelerator *accel);

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
    MeshData meshData;
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

  static void AddMeshData(MeshData *md, const Buffer *elembuf,
                          const Buffer *arraybuf, bool isDoublePrecisionPos,
                          const ArrayBufInfo *abinfo, GLsizei count,
                          GLuint offset);
  static void AddParticleData(MeshData *md, const Buffer *elembuf,
                              const Buffer *arraybuf, bool isDoublePrecisionPos,
                              const ArrayBufInfo *abinfo, GLsizei count,
                              GLuint offset, GLfloat constantWidth,
                              const GLfloat *widthBuf, GLsizei widthBufLen);
  static void AddLineData(MeshData *md, const Buffer *elembuf,
                          const Buffer *arraybuf, const ArrayBufInfo *abinfo,
                          GLsizei count, GLuint offset, GLfloat constantWidth,
                          const GLfloat *widthBuf, GLsizei widthBufLen,
                          bool cap);
  static void AddTetraData(MeshData *md, const Buffer *elembuf,
                           const Buffer *arraybuf, bool isDoublePrecisionPos,
                           const ArrayBufInfo *abinfo, GLsizei count,
                           GLuint offset);
  static void DebugDumpMesh(const Mesh *mesh);

  unsigned char *Locate(const Buffer *elembuf, const Buffer *arraybuf,
                        const ArrayBufInfo *abinfo, GLsizei count,
                        GLuint offset) const;
  void FreeMesh(CacheEntry *ce);

  std::vector<CacheEntry *> staticList_, dynamicList_;

  typedef std::map<unsigned char *, CacheEntry *> MeshDataMap;
  static MeshDataMap meshDataMap_;

  GLuint maxCacheSize_;
  GLuint cacheSizeUsed_;
};
}

#endif // __LSGL_GLES_ACCEL_BUILDER_H__
