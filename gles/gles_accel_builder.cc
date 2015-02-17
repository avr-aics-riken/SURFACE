/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cassert>
#include <cstring>

#include "gles_accel_builder.h"
#include "../render/timerutil.h"

using namespace lsgl;

AccelBuilder::MeshDataMap AccelBuilder::meshDataMap_;

AccelBuilder::AccelBuilder(GLuint maxCacheSize)
    : maxCacheSize_(maxCacheSize), cacheSizeUsed_(0) {}

AccelBuilder::~AccelBuilder() {
  // free all dynamic meshes
  EndFrame();

  // free all meshes in our static cache
  for (int t = 0; t < staticList_.size(); t++) {
    FreeMesh(staticList_[t]);
  }
  staticList_.clear();
}

AccelBuilder::MeshAccelerator *AccelBuilder::BuildMeshAccel(
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
  MeshAccelerator *accel = reinterpret_cast<MeshAccelerator *>(
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
  AddMeshData(&ce->meshData, elembuf, arraybuf, isDoublePrecisionPos, &abinfo,
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
  cacheSizeUsed_ += ce->meshData.GetByteSize();
  meshDataMap_[ce->meshData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] MeshBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<AccelBuilder::MeshAccelerator *>(ce->meshData.accel);
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
  AddParticleData(&ce->meshData, elembuf, arraybuf, isDoublePrecisionPos,
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
  cacheSizeUsed_ += ce->meshData.GetByteSize();
  meshDataMap_[ce->meshData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] ParticleBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<ParticleAccelerator *>(ce->meshData.accel);
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
  AddLineData(&ce->meshData, elembuf, arraybuf, &abinfo, count, offset,
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
  cacheSizeUsed_ += ce->meshData.GetByteSize();
  meshDataMap_[ce->meshData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("LineBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<LineAccelerator *>(ce->meshData.accel);
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
  AddTetraData(&ce->meshData, elembuf, arraybuf, isDoublePrecisionPos, &abinfo,
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
  cacheSizeUsed_ += ce->meshData.GetByteSize();
  meshDataMap_[ce->meshData.accel] = ce;

  // print current cache usage info and return a pointer to the accelerator
  printf("[LSGL] TetraBuilder cache: %dKB used / %dKB free (%.2f%%)\n",
         GetCacheUsed() / 1024, GetCacheSize() / 1024, GetCachePercent());
  return reinterpret_cast<AccelBuilder::TetraAccelerator *>(ce->meshData.accel);
}

bool AccelBuilder::AddMeshData(MeshData *md, MeshAccelerator *accel) {
  // locate the cache entry for this accelerator
  MeshDataMap::const_iterator it =
      meshDataMap_.find(reinterpret_cast<unsigned char *>(accel));
  if (it == meshDataMap_.end()) {
    return true;
  }

  // add the data for this mesh to the specified mesh data
  const CacheEntry *ce = it->second;
  AddMeshData(md, ce->elembuf, ce->arraybuf, ce->isDoublePrecisionPos,
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
      FreeMesh(ce);

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
      FreeMesh(ce);

      // move our index back one so we don't skip over any elements
      t--;
    }
  }
}

void AccelBuilder::EndFrame(void) {
  // free all meshes in the dynamic list
  for (int t = 0; t < dynamicList_.size(); t++) {
    FreeMesh(dynamicList_[t]);
  }
  dynamicList_.clear();

  // if cache used size is larger than the maximum allowed size, release some
  // low-hit static meshes
  while (GetCacheUsed() > GetCacheSize()) {
    // TODO
    break;
  }
}

clock_t AccelBuilder::GetBuiltTime(MeshAccelerator *accel) {
  MeshDataMap::const_iterator it =
      meshDataMap_.find(reinterpret_cast<unsigned char *>(accel));
  return (it != meshDataMap_.end()) ? it->second->builtTime : 0;
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

void AccelBuilder::AddMeshData(MeshData *md, const Buffer *elembuf,
                               const Buffer *arraybuf,
                               bool isDoublePrecisionPos,
                               const ArrayBufInfo *abinfo, GLsizei count,
                               GLuint offset) {
  GLuint maxIndex;

  // first free any existing mesh accelerator
  AccelBuilder::MeshAccelerator *macc =
      reinterpret_cast<AccelBuilder::MeshAccelerator *>(md->accel);
  delete macc;
  md->accel = NULL;

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->dpositionBuffer.GetCount() / 3, &maxIndex);

    // update the mesh structure with the data from our data buffers
    md->mesh.nfaces = md->indexBuffer.GetCount() / 3;
    md->mesh.faces = md->indexBuffer.GetBase();

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      // md->mesh.nvertices = arraybuf->GetSize() / (sizeof(double) * 3);
      md->mesh.nvertices = (maxIndex + 1);
      md->mesh.vertices = NULL;
      md->mesh.dvertices =
          reinterpret_cast<const double *>(arraybuf->GetData());
    } else {
      AddData(md->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      md->mesh.nvertices = md->dpositionBuffer.GetCount() / 3;
      md->mesh.vertices = NULL;
      md->mesh.dvertices = md->dpositionBuffer.GetBase();
    }

    md->mesh.isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->positionBuffer.GetCount() / 3, &maxIndex);

    // update the mesh structure with the data from our data buffers
    md->mesh.nfaces = md->indexBuffer.GetCount() / 3;
    md->mesh.faces = md->indexBuffer.GetBase();

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      // md->mesh.nvertices = arraybuf->GetSize() / (sizeof(float) * 3 * 3) / 3;
      md->mesh.nvertices = (maxIndex + 1);
      md->mesh.vertices = reinterpret_cast<const float *>(arraybuf->GetData());
      md->mesh.dvertices = NULL;

    } else {
      AddData(md->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      md->mesh.nvertices = md->positionBuffer.GetCount() / 3;
      md->mesh.vertices = md->positionBuffer.GetBase();
      md->mesh.dvertices = NULL;
    }

    md->mesh.isDoublePrecisionPos = false;
  }

  md->type = PRIMITIVE_TRIANGLES;

  TriangleBuildOptions options;
  printf("[LSGL] Double precision position = %d\n",
         md->mesh.isDoublePrecisionPos);
  timerutil t;
  t.start();
  TriangleAccel *accel = new TriangleAccel();
  accel->Build32(&md->mesh, options);
  t.end();
  printf("[LSGL] accel built time = %d msec\n", (int)t.msec());

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL] bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);
  md->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddParticleData(MeshData *md, const Buffer *elembuf,
                                   const Buffer *arraybuf,
                                   bool isDoublePrecisionPos,
                                   const ArrayBufInfo *abinfo, GLsizei count,
                                   GLuint offset, GLfloat constantWidth,
                                   const GLfloat *widthBuf,
                                   GLsizei widthBufLen) {
  GLuint maxIndex;

  // first free any existing mesh accelerator
  AccelBuilder::ParticleAccelerator *pacc =
      reinterpret_cast<AccelBuilder::ParticleAccelerator *>(md->accel);
  delete pacc;
  md->accel = NULL;

  if ((widthBufLen >= (count + offset)) && (widthBuf != NULL)) {
    GLfloat *dptr = md->radiusBuffer.Add(count);

    for (GLuint i = 0; i < count; i++) {
      dptr[i] = widthBuf[i + offset];
    }
  }

  Particles *part =
      new Particles(); // @fixme { No delete operaton for this object. }
  if ((widthBuf == NULL) || (widthBufLen == 0)) {
    part->radius = NULL;
  } else {
    part->radius = md->radiusBuffer.GetBase();
  }
  part->constantRadius = constantWidth * 0.5; // @fixme {}

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      part->dpositions = reinterpret_cast<const double *>(arraybuf->GetData());
      part->positions = NULL;
    } else {
      AddData(md->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      part->dpositions = md->dpositionBuffer.GetBase();
      part->positions = NULL;
    }

    part->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      part->positions = reinterpret_cast<const float *>(arraybuf->GetData());
      part->dpositions = NULL;

    } else {
      AddData(md->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      part->positions = md->positionBuffer.GetBase();
      part->dpositions = NULL;
    }

    part->isDoublePrecisionPos = false;
  }

  part->numParticles = md->indexBuffer.GetCount();

  md->mesh.nfaces = 0;
  md->mesh.faces = NULL;
  md->mesh.nvertices = 0;
  md->mesh.vertices = NULL;

  md->type = PRIMITIVE_POINTS;

  ParticleBuildOptions options;
  timerutil t;
  t.start();
  ParticleAccel *accel = new ParticleAccel();

  // accel->Build(part, options);
  accel->Build32(part, options);

  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Particle accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  md->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddLineData(MeshData *md, const Buffer *elembuf,
                               const Buffer *arraybuf,
                               const ArrayBufInfo *abinfo, GLsizei count,
                               GLuint offset, GLfloat constantWidth,
                               const GLfloat *widthBuf,
                               const GLsizei widthBufLen, bool cap) {
  GLuint maxIndex;

  // first free any existing accelerator
  AccelBuilder::LineAccelerator *pacc =
      reinterpret_cast<AccelBuilder::LineAccelerator *>(md->accel);
  delete pacc;
  md->accel = NULL;

  if ((widthBufLen >= (count + offset)) && (widthBuf != NULL)) {
    GLfloat *dptr = md->radiusBuffer.Add(count);

    for (GLuint i = 0; i < count; i++) {
      dptr[i] = widthBuf[i + offset];
    }
  }

  md->mesh.nfaces = 0;
  md->mesh.faces = NULL;
  md->mesh.nvertices = 0;
  md->mesh.vertices = NULL;
  md->type = PRIMITIVE_LINES;

  LineBuildOptions options;
  options.cap = cap;

  timerutil t;
  t.start();
  LineAccel *accel = new LineAccel();

  Lines *lines = new Lines(); // @fixme { No delete operaton for this object. }
  // lines->positions = md->positionBuffer.GetBase();
  if ((widthBuf == NULL) || (widthBufLen == 0)) {
    lines->radius = NULL;
  } else {
    lines->radius = md->radiusBuffer.GetBase();
  }
  lines->constantRadius = constantWidth * 0.5; // @fixme {}

  bool isDoublePrecisionPos = false; // @fixme;

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      lines->dpositions = reinterpret_cast<const double *>(arraybuf->GetData());
      lines->positions = NULL;
    } else {
      AddData(md->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      lines->dpositions = md->dpositionBuffer.GetBase();
      lines->positions = NULL;
    }

    lines->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      lines->positions = reinterpret_cast<const float *>(arraybuf->GetData());
      lines->dpositions = NULL;

    } else {
      AddData(md->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      lines->positions = md->positionBuffer.GetBase();
      lines->dpositions = NULL;
    }

    lines->isDoublePrecisionPos = false;
  }

  lines->numLines = md->indexBuffer.GetCount() / 2; // LINES

  accel->Build(lines, options);
  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Line accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  md->accel = reinterpret_cast<unsigned char *>(accel);
}

void AccelBuilder::AddTetraData(MeshData *md, const Buffer *elembuf,
                                const Buffer *arraybuf,
                                bool isDoublePrecisionPos,
                                const ArrayBufInfo *abinfo, GLsizei count,
                                GLuint offset) {
  GLuint maxIndex;

  // first free any existing accelerator
  AccelBuilder::TetraAccelerator *pacc =
      reinterpret_cast<AccelBuilder::TetraAccelerator *>(md->accel);
  delete pacc;
  md->accel = NULL;

  Tetrahedron *tetras =
      new Tetrahedron(); // @fixme { No delete operaton for this object. }

  if (isDoublePrecisionPos) {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->dpositionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      tetras->dvertices = reinterpret_cast<const double *>(arraybuf->GetData());
      tetras->vertices = NULL;
    } else {
      AddData(md->dpositionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);

      tetras->dvertices = md->dpositionBuffer.GetBase();
      tetras->vertices = NULL;
    }

    tetras->isDoublePrecisionPos = true;

  } else {
    // @todo { provide zero-copy version. }
    AddData(md->indexBuffer, elembuf, count, offset,
            md->positionBuffer.GetCount() / 3, &maxIndex);

    if (arraybuf->IsRetained()) {
      // zero-copy vesion.
      assert(abinfo->offsetPosition == 0);
      tetras->vertices = reinterpret_cast<const float *>(arraybuf->GetData());
      tetras->dvertices = NULL;

    } else {
      AddData(md->positionBuffer, arraybuf, (maxIndex + 1) * 3,
              abinfo->offsetPosition);
      tetras->vertices = md->positionBuffer.GetBase();
      tetras->dvertices = NULL;
    }

    tetras->isDoublePrecisionPos = false;
  }

  tetras->numTetrahedrons = md->indexBuffer.GetCount() / 4;
  tetras->faces = md->indexBuffer.GetBase();

  // Take a reference
  md->mesh.nfaces = tetras->numTetrahedrons;
  md->mesh.faces = tetras->faces;

  md->mesh.nvertices = 0;
  md->mesh.vertices = NULL;

  md->type = PRIMITIVE_TETRAHEDRONS;

  TetraBuildOptions options;
  timerutil t;
  t.start();
  TetraAccel *accel = new TetraAccel();

  accel->Build(tetras, options);
  // accel->Build32(tetras, options);

  t.end();

  double bmin[3], bmax[3];
  accel->BoundingBox(bmin, bmax);
  printf("[LSGL] Tetra accel built time = %d msec\n", (int)t.msec());
  printf("[LSGL]   bmin = %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("[LSGL]   bmax = %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  md->accel = reinterpret_cast<unsigned char *>(accel);
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
      return ce->meshData.accel;
    }
  }

  // not found...
  return NULL;
}

void AccelBuilder::FreeMesh(CacheEntry *ce) {
  // subtract the memory used for the mesh from the cache, and remove the mesh
  // from the map
  cacheSizeUsed_ -= ce->meshData.GetByteSize();

  //
  // Remove accel
  //
  if (ce->meshData.type == PRIMITIVE_TRIANGLES) {
    TriangleAccel *accel =
        reinterpret_cast<TriangleAccel *>(ce->meshData.accel);
    delete accel;
  } else if (ce->meshData.type == PRIMITIVE_POINTS) {
    AccelBuilder::ParticleAccelerator *accel =
        reinterpret_cast<AccelBuilder::ParticleAccelerator *>(
            ce->meshData.accel);
    delete accel;
  } else if (ce->meshData.type == PRIMITIVE_LINES) {
    AccelBuilder::LineAccelerator *accel =
        reinterpret_cast<AccelBuilder::LineAccelerator *>(ce->meshData.accel);
    delete accel;
  } else if (ce->meshData.type == PRIMITIVE_TETRAHEDRONS) {
    AccelBuilder::TetraAccelerator *accel =
        reinterpret_cast<AccelBuilder::TetraAccelerator *>(ce->meshData.accel);
    delete accel;
  } else {
    assert(0 && "Unknown error");
  }

  // free the cache entry
  delete ce;
}
