/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cstring>

#include "gles_accel_builder.h"
#include "../render/render_timerutil.h"

using namespace lsgl;

AccelBuilder::PrimDataMap AccelBuilder::primDataMap_;

AccelBuilder::AccelBuilder(GLuint maxCacheSize)
    : maxCacheSize_(maxCacheSize), cacheSizeUsed_(0) {}

AccelBuilder::~AccelBuilder() {
  // free all dynamic meshes
  EndFrame();

  // free all meshes in our static cache
  for (int t = 0; t < staticList_.size(); t++) {
    FreePrim(staticList_[t]);
  }
  staticList_.clear();
}

AccelBuilder::TriangleAccelerator *AccelBuilder::BuildTriangleAccel(
    const Buffer *elembuf, const Buffer *arraybuf, bool isDoublePrecisionPos,
    const VertexAttribute *vertexAttributes, const GLuint *texture2D,
    GLsizei count, GLuint offset) {
  // ensure the position vertex attribute is enabled, at minimum
  const VertexAttribute *attrPos = &vertexAttributes[kVtxAttrPosition];
  if (attrPos->enabled == false) {
    return NULL;
  }

  // initialize array buffer info
  ArrayBufInfo abinfo = {0};
  abinfo.offsetPosition = (GLubyte *)attrPos->ptr - (GLubyte *)NULL;

  // FIXME: assume non-interleaved XYZ float for position (for now)
  if ((attrPos->size == 3) && (attrPos->type == GL_FLOAT) &&
      (attrPos->normalized == GL_FALSE) &&
      (attrPos->stride == (3 * sizeof(float)))) {
    // OK
  } else if ((attrPos->size == 3) && (attrPos->type == GL_DOUBLE) &&
             (attrPos->normalized == GL_FALSE) &&
             (attrPos->stride == (3 * sizeof(double)))) {
    // OK. Double preicsion
    assert(isDoublePrecisionPos == true);
  } else {
    fprintf(stderr, "[LSGL] Unsupported primitive format.\n");
    return NULL;
  }
  // if ((attrPos->size != 3) || ((attrPos->type != GL_FLOAT) && (attrPos->type
  // != GL_DOUBLE)) ||
  //    (attrPos->normalized != GL_FALSE) ||
  //    (attrPos->stride != (3 * sizeof(float)))) {
  //  assert((attrPos->size == 3) && "vtxattr position size should be 3");
  //  assert((attrPos->type == GL_FLOAT) &&
  //         "vtxattr position type should be GL_FLOAT");
  //  assert((attrPos->normalized == GL_FALSE) &&
  //         "vtxattr position normalized should be GL_FALSE");
  //  assert((attrPos->stride == (3 * sizeof(float))) &&
  //         "vtxattr position stride must be equal to 3*sizeof(float)");
  //  return NULL;
  //}

  // first try to locate this mesh in our static cache
  TriangleAccelerator *accel = reinterpret_cast<TriangleAccelerator *>(
      Locate(elembuf, arraybuf, &abinfo, count, offset));
  if (accel != NULL) {
    return accel;
  }

  // allocate a new cache entry and initialize
  CacheEntry *ce = new CacheEntry;
  ce->elembuf = elembuf;
  ce->arraybuf = arraybuf;
  ce->isDoublePrecisionPos = isDoublePrecisionPos;
  ce->abinfo = abinfo;
  ce->count = count;
  ce->offset = offset;
  ce->hitCount = 1;
  ce->builtTime = clock();

  // build the mesh accelerator
  AddTriangleData(&ce->primData, elembuf, arraybuf, isDoublePrecisionPos,
                  &abinfo, count, offset);

  // if both buffers are marked as static, add this mesh to the static cache
  // list, otherwise add it to the dynamic list instead
  if ((elembuf->IsStatic() == true) && (arraybuf->IsStatic() == true)) {
    staticList_.push_back(ce);
    ce->dynamic = false;
  } else {
    dynamicList_.push_back(ce);
    ce->dynamic = true;
  }

  // update used cache size, and add the mesh accelerator to our map
  cacheSizeUsed_ += ce->primData.GetByteSize();
  primDataMap_[ce->primData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] MeshBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<AccelBuilder::TriangleAccelerator *>(
      ce->primData.accel);
}

AccelBuilder::ParticleAccelerator *AccelBuilder::BuildParticleAccel(
    const Buffer *elembuf, const Buffer *arraybuf, bool isDoublePrecisionPos,
    const VertexAttribute *vertexAttributes, GLsizei count, GLuint offset,
    GLfloat constantWidth, const GLfloat *widthBuf, GLsizei widthBufLen) {
  // ensure the position vertex attribute is enabled, at minimum
  const VertexAttribute *attrPos = &vertexAttributes[kVtxAttrPosition];
  if (attrPos->enabled == false) {
    return NULL;
  }

  // initialize array buffer info
  ArrayBufInfo abinfo = {0};
  abinfo.offsetPosition = (GLubyte *)attrPos->ptr - (GLubyte *)NULL;

  // FIXME: assume non-interleaved XYZ float for position (for now)
  if ((attrPos->size != 3) || (attrPos->type != GL_FLOAT) ||
      (attrPos->normalized != GL_FALSE) ||
      (attrPos->stride != (3 * sizeof(float)))) {
    assert((attrPos->size == 3) && "vtxattr position size should be 3");
    assert((attrPos->type == GL_FLOAT) &&
           "vtxattr position type should be GL_FLOAT");
    assert((attrPos->normalized == GL_FALSE) &&
           "vtxattr position normalized should be GL_FALSE");
    assert((attrPos->stride == (3 * sizeof(float))) &&
           "vtxattr position stride must be equal to 3*sizeof(float)");
    return NULL;
  }

  // first try to locate this mesh in our static cache
  ParticleAccelerator *accel = reinterpret_cast<ParticleAccelerator *>(
      Locate(elembuf, arraybuf, &abinfo, count, offset));
  if (accel != NULL) {
    return accel;
  }

  // allocate a new cache entry and initialize
  CacheEntry *ce = new CacheEntry;
  ce->elembuf = elembuf;
  ce->arraybuf = arraybuf;
  ce->abinfo = abinfo;
  ce->count = count;
  ce->offset = offset;
  ce->hitCount = 1;
  ce->builtTime = clock();

  // build the particle accelerator
  AddParticleData(&ce->primData, elembuf, arraybuf, isDoublePrecisionPos,
                  &abinfo, count, offset, constantWidth, widthBuf, widthBufLen);

  // if both buffers are marked as static, add this mesh to the static cache
  // list, otherwise add it to the dynamic list instead
  if ((elembuf->IsStatic() == true) && (arraybuf->IsStatic() == true)) {
    staticList_.push_back(ce);
    ce->dynamic = false;
  } else {
    dynamicList_.push_back(ce);
    ce->dynamic = true;
  }

  // update used cache size, and add the mesh accelerator to our map
  cacheSizeUsed_ += ce->primData.GetByteSize();
  primDataMap_[ce->primData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] ParticleBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<ParticleAccelerator *>(ce->primData.accel);
}

AccelBuilder::LineAccelerator *AccelBuilder::BuildLineAccel(
    const Buffer *elembuf, const Buffer *arraybuf, bool isDoublePrecisionPos,
    const VertexAttribute *vertexAttributes, GLsizei count, GLuint offset,
    GLfloat constantWidth, const GLfloat *widthBuf, const GLsizei widthBufLen,
    bool cap) {
  // ensure the position vertex attribute is enabled, at minimum
  const VertexAttribute *attrPos = &vertexAttributes[kVtxAttrPosition];
  if (attrPos->enabled == false) {
    return NULL;
  }

  // initialize array buffer info
  ArrayBufInfo abinfo = {0};
  abinfo.offsetPosition = (GLubyte *)attrPos->ptr - (GLubyte *)NULL;

  // FIXME: assume non-interleaved XYZ float for position (for now)
  if ((attrPos->size != 3) || (attrPos->type != GL_FLOAT) ||
      (attrPos->normalized != GL_FALSE) ||
      (attrPos->stride != (3 * sizeof(float)))) {
    assert((attrPos->size == 3) && "vtxattr position size should be 3");
    assert((attrPos->type == GL_FLOAT) &&
           "vtxattr position type should be GL_FLOAT");
    assert((attrPos->normalized == GL_FALSE) &&
           "vtxattr position normalized should be GL_FALSE");
    assert((attrPos->stride == (3 * sizeof(float))) &&
           "vtxattr position stride must be equal to 3*sizeof(float)");
    return NULL;
  }

  // first try to locate this line in our static cache
  LineAccelerator *accel = reinterpret_cast<LineAccelerator *>(
      Locate(elembuf, arraybuf, &abinfo, count, offset));
  if (accel != NULL) {
    return accel;
  }

  // allocate a new cache entry and initialize
  CacheEntry *ce = new CacheEntry;
  ce->elembuf = elembuf;
  ce->arraybuf = arraybuf;
  ce->abinfo = abinfo;
  ce->count = count;
  ce->offset = offset;
  ce->hitCount = 1;
  ce->builtTime = clock();

  // build the particle accelerator
  AddLineData(&ce->primData, elembuf, arraybuf, &abinfo, count, offset,
              constantWidth, widthBuf, widthBufLen, cap);

  // if both buffers are marked as static, add this mesh to the static cache
  // list, otherwise add it to the dynamic list instead
  if ((elembuf->IsStatic() == true) && (arraybuf->IsStatic() == true)) {
    staticList_.push_back(ce);
    ce->dynamic = false;
  } else {
    dynamicList_.push_back(ce);
    ce->dynamic = true;
  }

  // update used cache size, and add the mesh accelerator to our map
  cacheSizeUsed_ += ce->primData.GetByteSize();
  primDataMap_[ce->primData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] LineBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<LineAccelerator *>(ce->primData.accel);
}

AccelBuilder::TetraAccelerator *AccelBuilder::BuildTetraAccel(
    const Buffer *elembuf, const Buffer *arraybuf, bool isDoublePrecisionPos,
    const VertexAttribute *vertexAttributes, GLsizei count, GLuint offset) {
  // ensure the position vertex attribute is enabled, at minimum
  const VertexAttribute *attrPos = &vertexAttributes[kVtxAttrPosition];
  if (attrPos->enabled == false) {
    return NULL;
  }

  // initialize array buffer info
  ArrayBufInfo abinfo = {0};
  abinfo.offsetPosition = (GLubyte *)attrPos->ptr - (GLubyte *)NULL;

  // FIXME: assume non-interleaved XYZ float for position (for now)
  if ((attrPos->size == 3) && (attrPos->type == GL_FLOAT) &&
      (attrPos->normalized == GL_FALSE) &&
      (attrPos->stride == (3 * sizeof(float)))) {
    // OK
  } else if ((attrPos->size == 3) && (attrPos->type == GL_DOUBLE) &&
             (attrPos->normalized == GL_FALSE) &&
             (attrPos->stride == (3 * sizeof(double)))) {
    // OK. Double preicsion
    assert(isDoublePrecisionPos == true);
  } else {
    fprintf(stderr, "[LSGL] Unsupported primitive format.\n");
    return NULL;
  }

  // first try to locate this mesh in our static cache
  TetraAccelerator *accel = reinterpret_cast<TetraAccelerator *>(
      Locate(elembuf, arraybuf, &abinfo, count, offset));
  if (accel != NULL) {
    return accel;
  }

  // allocate a new cache entry and initialize
  CacheEntry *ce = new CacheEntry;
  ce->elembuf = elembuf;
  ce->arraybuf = arraybuf;
  ce->isDoublePrecisionPos = isDoublePrecisionPos;
  ce->abinfo = abinfo;
  ce->count = count;
  ce->offset = offset;
  ce->hitCount = 1;
  ce->builtTime = clock();

  // build the mesh accelerator
  AddTetraData(&ce->primData, elembuf, arraybuf, isDoublePrecisionPos, &abinfo,
               count, offset);

  // if both buffers are marked as static, add this mesh to the static cache
  // list, otherwise add it to the dynamic list instead
  if ((elembuf->IsStatic() == true) && (arraybuf->IsStatic() == true)) {
    staticList_.push_back(ce);
    ce->dynamic = false;
  } else {
    dynamicList_.push_back(ce);
    ce->dynamic = true;
  }

  // update used cache size, and add the mesh accelerator to our map
  cacheSizeUsed_ += ce->primData.GetByteSize();
  primDataMap_[ce->primData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] TetraBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<AccelBuilder::TetraAccelerator *>(ce->primData.accel);
}

AccelBuilder::SolidAccelerator *AccelBuilder::BuildSolidAccel(
    GLenum solidType,
    const Buffer *elembuf, const Buffer *arraybuf, bool isDoublePrecisionPos,
    const VertexAttribute *vertexAttributes, GLsizei count, GLuint offset) {
  // ensure the position vertex attribute is enabled, at minimum
  const VertexAttribute *attrPos = &vertexAttributes[kVtxAttrPosition];
  if (attrPos->enabled == false) {
    return NULL;
  }

  // initialize array buffer info
  ArrayBufInfo abinfo = {0};
  abinfo.offsetPosition = (GLubyte *)attrPos->ptr - (GLubyte *)NULL;

  // FIXME: assume non-interleaved XYZ float for position (for now)
  if ((attrPos->size == 3) && (attrPos->type == GL_FLOAT) &&
      (attrPos->normalized == GL_FALSE) &&
      (attrPos->stride == (3 * sizeof(float)))) {
    // OK
  } else if ((attrPos->size == 3) && (attrPos->type == GL_DOUBLE) &&
             (attrPos->normalized == GL_FALSE) &&
             (attrPos->stride == (3 * sizeof(double)))) {
    // OK. Double preicsion
    assert(isDoublePrecisionPos == true);
  } else {
    fprintf(stderr, "[LSGL] Unsupported primitive format.\n");
    return NULL;
  }

  // first try to locate this mesh in our static cache
  SolidAccelerator *accel = reinterpret_cast<SolidAccelerator *>(
      Locate(elembuf, arraybuf, &abinfo, count, offset));
  if (accel != NULL) {
    return accel;
  }

  // allocate a new cache entry and initialize
  CacheEntry *ce = new CacheEntry;
  ce->elembuf = elembuf;
  ce->arraybuf = arraybuf;
  ce->isDoublePrecisionPos = isDoublePrecisionPos;
  ce->abinfo = abinfo;
  ce->count = count;
  ce->offset = offset;
  ce->hitCount = 1;
  ce->builtTime = clock();

  // build the mesh accelerator
  AddSolidData(solidType, &ce->primData, elembuf, arraybuf, isDoublePrecisionPos, &abinfo,
               count, offset);

  // if both buffers are marked as static, add this mesh to the static cache
  // list, otherwise add it to the dynamic list instead
  if ((elembuf->IsStatic() == true) && (arraybuf->IsStatic() == true)) {
    staticList_.push_back(ce);
    ce->dynamic = false;
  } else {
    dynamicList_.push_back(ce);
    ce->dynamic = true;
  }

  // update used cache size, and add the mesh accelerator to our map
  cacheSizeUsed_ += ce->primData.GetByteSize();
  primDataMap_[ce->primData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] SolidBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<AccelBuilder::SolidAccelerator *>(ce->primData.accel);
}

bool AccelBuilder::AddTriangleData(PrimData *pd, TriangleAccelerator *accel) {
  // locate the cache entry for this accelerator
  PrimDataMap::const_iterator it =
      primDataMap_.find(reinterpret_cast<unsigned char *>(accel));
  if (it == primDataMap_.end()) {
    return true;
  }

  // add the data for this mesh to the specified mesh data
  const CacheEntry *ce = it->second;
  AddTriangleData(pd, ce->elembuf, ce->arraybuf, ce->isDoublePrecisionPos,
                  &ce->abinfo, ce->count, ce->offset);
  return ce->dynamic;
}

void AccelBuilder::Invalidate(const Buffer *buf) {
  // check if this buffer is used by any meshes in the static cache first
  for (int t = 0; t < staticList_.size(); t++) {
    CacheEntry *ce = staticList_[t];
    if ((buf == ce->elembuf) || (buf == ce->arraybuf)) {
      // remove this mesh
      staticList_.erase(staticList_.begin() + t);
      FreePrim(ce);

      // move our index back one so we don't skip over any elements
      t--;
    }
  }

  // next check if this buffer is used by any meshes in the dynamic list next
  for (int t = 0; t < dynamicList_.size(); t++) {
    CacheEntry *ce = dynamicList_[t];
    if ((buf == ce->elembuf) || (buf == ce->arraybuf)) {
      // remove this mesh
      dynamicList_.erase(dynamicList_.begin() + t);
      FreePrim(ce);

      // move our index back one so we don't skip over any elements
      t--;
    }
  }
}

void AccelBuilder::EndFrame(void) {
  // free all meshes in the dynamic list
  for (int t = 0; t < dynamicList_.size(); t++) {
    FreePrim(dynamicList_[t]);
  }
  dynamicList_.clear();

  // if cache used size is larger than the maximum allowed size, release some
  // low-hit static meshes
  while (GetCacheUsed() > GetCacheSize()) {
    // TODO
    break;
  }
}

clock_t AccelBuilder::GetBuiltTime(TriangleAccelerator *accel) {
  PrimDataMap::const_iterator it =
      primDataMap_.find(reinterpret_cast<unsigned char *>(accel));
  return (it != primDataMap_.end()) ? it->second->builtTime : 0;
}

template <typename T>
void AccelBuilder::AddData(DataBuffer<T> &dest, const Buffer *src, GLuint count,
                           GLuint offset, GLuint adjust, T *maxValue) {
  // request memory for new data elements
  const T *sptr = reinterpret_cast<const T *>(src->GetData());
  T *dptr = dest.Add(count);

  // apply source offset
  sptr = reinterpret_cast<const T *>(reinterpret_cast<const GLubyte *>(sptr) +
                                     offset);

  // determine maximum value, if desired
  if (maxValue != NULL) {
    *maxValue = sptr[0];
    for (GLuint t = 1; t < count; t++) {
      if (sptr[t] > *maxValue) {
        *maxValue = sptr[t];
      }
    }
  }

  // if no adjustment is required, just use memcpy
  if (adjust == 0) {
    memcpy(dptr, sptr, count * sizeof(T));
  } else {
    // adjust each element by specified amount
    for (GLuint t = 0; t < count; t++) {
      dptr[t] = sptr[t] + adjust;
    }
  }
}

void AccelBuilder::AddTriangleData(PrimData *pd, const Buffer *elembuf,
                                   const Buffer *arraybuf,
                                   bool isDoublePrecisionPos,
                                   const ArrayBufInfo *abinfo, GLsizei count,
                                   GLuint offset) {
  GLuint maxIndex;

  // first free any existing mesh accelerator
  AccelBuilder::TriangleAccelerator *macc =
      reinterpret_cast<AccelBuilder::TriangleAccelerator *>(pd->accel);
  delete macc;
  pd->accel = NULL;

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->dpositionBuffer.GetCount() / 3, &maxIndex);

    // update the mesh structure with the data from our data buffers
    pd->mesh.nfaces = pd->indexBuffer.GetCount() / 3;
    pd->mesh.faces = pd->indexBuffer.GetBase();

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      // pd->mesh.nvertices = arraybuf->GetSize() / (sizeof(double) * 3);
      pd->mesh.nvertices = (maxIndex + 1);
      pd->mesh.vertices = NULL;
      pd->mesh.dvertices =
          reinterpret_cast<const double *>(arraybuf->GetData());
    } else {
      AddData(pd->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      pd->mesh.nvertices = pd->dpositionBuffer.GetCount() / 3;
      pd->mesh.vertices = NULL;
      pd->mesh.dvertices = pd->dpositionBuffer.GetBase();
    }

    pd->mesh.isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->positionBuffer.GetCount() / 3, &maxIndex);

    // update the mesh structure with the data from our data buffers
    pd->mesh.nfaces = pd->indexBuffer.GetCount() / 3;
    pd->mesh.faces = pd->indexBuffer.GetBase();

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      // pd->mesh.nvertices = arraybuf->GetSize() / (sizeof(float) * 3 * 3) / 3;
      pd->mesh.nvertices = (maxIndex + 1);
      pd->mesh.vertices = reinterpret_cast<const float *>(arraybuf->GetData());
      pd->mesh.dvertices = NULL;

    } else {
      AddData(pd->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      pd->mesh.nvertices = pd->positionBuffer.GetCount() / 3;
      pd->mesh.vertices = pd->positionBuffer.GetBase();
      pd->mesh.dvertices = NULL;
    }

    pd->mesh.isDoublePrecisionPos = false;
  }

  pd->type = PRIMITIVE_TRIANGLES;

  TriangleBuildOptions options;
  printf("[LSGL] Double precision position = %d\n",
         pd->mesh.isDoublePrecisionPos);
  timerutil t;
  t.start();
  TriangleAccel *accel = new TriangleAccel();
  accel->Build32(&pd->mesh, options);
  t.end();
  printf("[LSGL] accel built time = %d msec\n", (int)t.msec());

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL] bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  pd->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddParticleData(PrimData *pd, const Buffer *elembuf,
                                   const Buffer *arraybuf,
                                   bool isDoublePrecisionPos,
                                   const ArrayBufInfo *abinfo, GLsizei count,
                                   GLuint offset, GLfloat constantWidth,
                                   const GLfloat *widthBuf,
                                   GLsizei widthBufLen) {
  GLuint maxIndex;

  // first free any existing mesh accelerator
  AccelBuilder::ParticleAccelerator *pacc =
      reinterpret_cast<AccelBuilder::ParticleAccelerator *>(pd->accel);
  delete pacc;
  pd->accel = NULL;

  if ((widthBufLen >= (count + offset)) && (widthBuf != NULL)) {
    GLfloat *dptr = pd->radiusBuffer.Add(count);

    for (GLuint i = 0; i < count; i++) {
      dptr[i] = widthBuf[i + offset];
    }
  }

  delete pd->points;
  pd->points = new Particles();
  if ((widthBuf == NULL) || (widthBufLen == 0)) {
    pd->points->radius = NULL;
  } else {
    pd->points->radius = pd->radiusBuffer.GetBase();
  }
  pd->points->constantRadius = constantWidth * 0.5; // @fixme {}

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->points->dpositions =
          reinterpret_cast<const double *>(arraybuf->GetData());
      pd->points->positions = NULL;
    } else {
      AddData(pd->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      pd->points->dpositions = pd->dpositionBuffer.GetBase();
      pd->points->positions = NULL;
    }

    pd->points->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->points->positions =
          reinterpret_cast<const float *>(arraybuf->GetData());
      pd->points->dpositions = NULL;

    } else {
      AddData(pd->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      pd->points->positions = pd->positionBuffer.GetBase();
      pd->points->dpositions = NULL;
    }

    pd->points->isDoublePrecisionPos = false;
  }

  pd->points->numParticles = pd->indexBuffer.GetCount();

  pd->mesh.nfaces = 0;
  pd->mesh.faces = NULL;
  pd->mesh.nvertices = 0;
  pd->mesh.vertices = NULL;

  pd->type = PRIMITIVE_POINTS;

  ParticleBuildOptions options;
  timerutil t;
  t.start();
  ParticleAccel *accel = new ParticleAccel();

  // accel->Build(pd->points, options);
  accel->Build32(pd->points, options);

  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Particle accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  pd->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddLineData(PrimData *pd, const Buffer *elembuf,
                               const Buffer *arraybuf,
                               const ArrayBufInfo *abinfo, GLsizei count,
                               GLuint offset, GLfloat constantWidth,
                               const GLfloat *widthBuf,
                               const GLsizei widthBufLen, bool cap) {
  GLuint maxIndex;

  // first free any existing accelerator
  AccelBuilder::LineAccelerator *pacc =
      reinterpret_cast<AccelBuilder::LineAccelerator *>(pd->accel);
  delete pacc;
  pd->accel = NULL;

  if ((widthBufLen >= (count + offset)) && (widthBuf != NULL)) {
    GLfloat *dptr = pd->radiusBuffer.Add(count);

    for (GLuint i = 0; i < count; i++) {
      dptr[i] = widthBuf[i + offset];
    }
  }

  pd->mesh.nfaces = 0;
  pd->mesh.faces = NULL;
  pd->mesh.nvertices = 0;
  pd->mesh.vertices = NULL;
  pd->type = PRIMITIVE_LINES;

  LineBuildOptions options;
  options.cap = cap;

  timerutil t;
  t.start();
  LineAccel *accel = new LineAccel();

  pd->lines = new Lines();
  // lines->positions = pd->positionBuffer.GetBase();
  if ((widthBuf == NULL) || (widthBufLen == 0)) {
    pd->lines->radius = NULL;
  } else {
    pd->lines->radius = pd->radiusBuffer.GetBase();
  }
  pd->lines->constantRadius = constantWidth * 0.5; // @fixme {}

  bool isDoublePrecisionPos = false; // @fixme;

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->lines->dpositions =
          reinterpret_cast<const double *>(arraybuf->GetData());
      pd->lines->positions = NULL;
    } else {
      AddData(pd->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      pd->lines->dpositions = pd->dpositionBuffer.GetBase();
      pd->lines->positions = NULL;
    }

    pd->lines->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->lines->positions =
          reinterpret_cast<const float *>(arraybuf->GetData());
      pd->lines->dpositions = NULL;

    } else {
      AddData(pd->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      pd->lines->positions = pd->positionBuffer.GetBase();
      pd->lines->dpositions = NULL;
    }

    pd->lines->isDoublePrecisionPos = false;
  }

  pd->lines->numLines = pd->indexBuffer.GetCount() / 2; // LINES
  pd->lines->segments = pd->indexBuffer.GetBase();

  accel->Build(
      pd->lines,
      options); // @todo { implement parallel BVH builder for line primitive. }
  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Line accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  pd->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddTetraData(PrimData *pd, const Buffer *elembuf,
                                const Buffer *arraybuf,
                                bool isDoublePrecisionPos,
                                const ArrayBufInfo *abinfo, GLsizei count,
                                GLuint offset) {
  GLuint maxIndex;

  // first free any existing accelerator
  AccelBuilder::TetraAccelerator *pacc =
      reinterpret_cast<AccelBuilder::TetraAccelerator *>(pd->accel);
  delete pacc;
  pd->accel = NULL;

  delete pd->tetras;
  pd->tetras = new Tetrahedron();

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->tetras->dvertices =
          reinterpret_cast<const double *>(arraybuf->GetData());
      pd->tetras->vertices = NULL;
    } else {
      AddData(pd->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      pd->tetras->dvertices = pd->dpositionBuffer.GetBase();
      pd->tetras->vertices = NULL;
    }

    pd->tetras->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->tetras->vertices =
          reinterpret_cast<const float *>(arraybuf->GetData());
      pd->tetras->dvertices = NULL;

    } else {
      AddData(pd->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      pd->tetras->vertices = pd->positionBuffer.GetBase();
      pd->tetras->dvertices = NULL;
    }

    pd->tetras->isDoublePrecisionPos = false;
  }

  pd->tetras->numTetrahedrons = pd->indexBuffer.GetCount() / 4;
  pd->tetras->faces = pd->indexBuffer.GetBase();

  // Take a reference
  pd->mesh.nfaces = 0;
  pd->mesh.faces = NULL;
  pd->mesh.nvertices = 0;
  pd->mesh.vertices = NULL;
  pd->type = PRIMITIVE_TETRAHEDRONS;

  TetraBuildOptions options;
  timerutil t;
  t.start();
  TetraAccel *accel = new TetraAccel();

  // accel->Build(pd->tetras, options);  // slower but rather high quality BVH
  accel->Build32(pd->tetras, options); // faster but rather low quality BVH

  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Tetra accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  pd->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddSolidData(GLenum solidType, PrimData *pd, const Buffer *elembuf,
                                const Buffer *arraybuf,
                                bool isDoublePrecisionPos,
                                const ArrayBufInfo *abinfo, GLsizei count,
                                GLuint offset) {
  GLuint maxIndex;

  // first free any existing accelerator
  AccelBuilder::SolidAccelerator *pacc =
      reinterpret_cast<AccelBuilder::SolidAccelerator *>(pd->accel);
  delete pacc;
  pd->accel = NULL;

  delete pd->solids;
  pd->solids = new Solid();

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->solids->dvertices =
          reinterpret_cast<const double *>(arraybuf->GetData());
      pd->solids->vertices = NULL;
    } else {
      AddData(pd->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      pd->solids->dvertices = pd->dpositionBuffer.GetBase();
      pd->solids->vertices = NULL;
    }

    pd->solids->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(pd->indexBuffer, elembuf, count, offset,
            pd->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      pd->solids->vertices =
          reinterpret_cast<const float *>(arraybuf->GetData());
      pd->solids->dvertices = NULL;

    } else {
      AddData(pd->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      pd->solids->vertices = pd->positionBuffer.GetBase();
      pd->solids->dvertices = NULL;
    }

    pd->solids->isDoublePrecisionPos = false;
  }

  if (solidType == GL_PYRAMIDS_EXT) {
    pd->solids->numSolids = pd->indexBuffer.GetCount() / 5;
    pd->solids->numVertsPerSolid = 5;
    pd->type = PRIMITIVE_PYRAMIDS;
  } else if (solidType == GL_PRISMS_EXT) {
    pd->solids->numSolids = pd->indexBuffer.GetCount() / 6;
    pd->solids->numVertsPerSolid = 6;
    pd->type = PRIMITIVE_PRISMS;
  } else if (solidType == GL_HEXAHEDRONS_EXT) {
    pd->solids->numSolids = pd->indexBuffer.GetCount() / 8;
    pd->solids->numVertsPerSolid = 8;
    pd->type = PRIMITIVE_HEXAHEDRONS;
  } else {
    assert(0);
  }
  pd->solids->indices = pd->indexBuffer.GetBase();

  // Take a reference
  pd->mesh.nfaces = 0;
  pd->mesh.faces = NULL;
  pd->mesh.nvertices = 0;
  pd->mesh.vertices = NULL;

  SolidBuildOptions options;
  timerutil t;
  t.start();
  SolidAccel *accel = new SolidAccel();

  accel->Build(pd->solids, options);  // slower but rather high quality BVH
  //accel->Build32(pd->solids, options); // faster but rather low quality BVH

  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Solid accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  pd->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::DebugDumpMesh(const Mesh *mesh) {
  int t;

  printf("mesh debug dump:\n");
  printf("---------------\n\n");

  printf("number of faces: %lu\n", (unsigned long)mesh->nfaces);
  printf("face buffer:\n");
  for (t = 0; t < mesh->nfaces * 3; t++) {
    printf("%d", mesh->faces[t]);
    if ((t % 3) != 2) {
      printf(" - ");
    } else {
      printf(",\n");
    }
  }
  printf("\n\n");

  printf("number of vertices: %lu\n", (unsigned long)mesh->nvertices);
  printf("vertex buffer:\n");
  for (t = 0; t < mesh->nvertices; t++) {
    printf("%d: pos=(%.3f, %.3f, %.3f),", t, mesh->vertices[t * 3],
           mesh->vertices[(t * 3) + 1], mesh->vertices[(t * 3) + 2]);
    if (mesh->normals != NULL) {
      printf(" nml=(%.3f, %.3f, %.3f),", mesh->normals[t * 3],
             mesh->normals[(t * 3) + 1], mesh->normals[(t * 3) + 2]);
    }
    // if (mesh->texcoords != NULL) {
    //  printf(" tex=(%.3f, %.3f),", mesh->texcoords[t * 2],
    //         mesh->texcoords[(t * 2) + 1]);
    //}
    printf("\n");
  }
  printf("\n\n");
}

unsigned char *AccelBuilder::Locate(const Buffer *elembuf,
                                    const Buffer *arraybuf,
                                    const ArrayBufInfo *abinfo, GLsizei count,
                                    GLuint offset) const {
  // iterate for all cached meshes
  for (int t = 0; t < staticList_.size(); t++) {
    CacheEntry *ce = staticList_[t];

    // check if this cached mesh matches the specified parameters
    if ((ce->elembuf == elembuf) && (ce->arraybuf == arraybuf) &&
        (ce->abinfo == *abinfo) && (ce->count == count) &&
        (ce->offset == offset)) {
      // increment the hit count and return the already built accelerator
      ce->hitCount++;
      return ce->primData.accel;
    }
  }

  // not found...
  return NULL;
}

void AccelBuilder::FreePrim(CacheEntry *ce) {
  // subtract the memory used for the mesh from the cache, and remove the mesh
  // from the map
  cacheSizeUsed_ -= ce->primData.GetByteSize();

  //
  // Remove accel
  //
  if (ce->primData.type == PRIMITIVE_TRIANGLES) {
    TriangleAccel *accel =
        reinterpret_cast<TriangleAccel *>(ce->primData.accel);
    delete accel;
  } else if (ce->primData.type == PRIMITIVE_POINTS) {
    delete ce->primData.points;
    AccelBuilder::ParticleAccelerator *accel =
        reinterpret_cast<AccelBuilder::ParticleAccelerator *>(
            ce->primData.accel);
    delete accel;
  } else if (ce->primData.type == PRIMITIVE_LINES) {
    delete ce->primData.lines;
    AccelBuilder::LineAccelerator *accel =
        reinterpret_cast<AccelBuilder::LineAccelerator *>(ce->primData.accel);
    delete accel;
  } else if (ce->primData.type == PRIMITIVE_TETRAHEDRONS) {
    delete ce->primData.tetras;
    AccelBuilder::TetraAccelerator *accel =
        reinterpret_cast<AccelBuilder::TetraAccelerator *>(ce->primData.accel);
    delete accel;
  } else if ((ce->primData.type == PRIMITIVE_PYRAMIDS) || 
             (ce->primData.type == PRIMITIVE_PRISMS) || 
             (ce->primData.type == PRIMITIVE_HEXAHEDRONS)) {
    delete ce->primData.solids;
    AccelBuilder::SolidAccelerator *accel =
        reinterpret_cast<AccelBuilder::SolidAccelerator *>(ce->primData.accel);
    delete accel;
  } else {
    assert(0 && "Unknown error");
  }

  // free the cache entry
  delete ce;
}
