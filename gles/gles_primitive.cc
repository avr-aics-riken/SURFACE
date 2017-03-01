/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2016 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include "GLES2/gl2.h"

#include <cassert>

#include "gles_context.h"
#include "gles_render_graph.h"

using namespace lsgl;

void Context::glDrawArrays(GLenum mode, GLint first, GLsizei count) {
  TRACE_EVENT("(GLenum mode = %d, GLint first = %d, GLsizei count = %d", mode,
              first, count);

  if (count < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (count == 0) {
    return;
  }

  if ((mode != GL_TRIANGLES) && (mode != GL_POINTS) && (mode != GL_LINES) &&
      (mode != GL_TETRAHEDRONS_EXT) &&
      (mode != GL_PYRAMIDS_EXT) &&
      (mode != GL_PRISMS_EXT) &&
      (mode != GL_HEXAHEDRONS_EXT)) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // lookup the currently bound array buffers
  Buffer *arraybuf = HandleToBuffer(GL_ARRAY_BUFFER);
  if (arraybuf == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // Construct index buffer
  std::vector<GLint> indices;

  int n = 0;

  // Truncate # of indices depending on primitive type.
  if (mode == GL_TRIANGLES) {
    n = (count / 3) * 3;
  } else if (mode == GL_POINTS) {
    n = count;
  } else if (mode == GL_LINES) {
    n = (count / 2) * 2;
  } else if (mode == GL_TETRAHEDRONS_EXT) {
    n = (count / 4) * 4;
  } else if (mode == GL_PYRAMIDS_EXT) {
    n = (count / 5) * 5;
  } else if (mode == GL_PRISMS_EXT) {
    n = (count / 6) * 6;
  } else if (mode == GL_HEXAHEDRONS_EXT) {
    n = (count / 8) * 8;
  }

  indices.resize(n);
  for (int i = 0; i < n; i++) {
    indices[i] = i + first;
  }

  // Cerate index buffer.
  GLuint idx;
  glGenBuffers(1, &idx);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, idx);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLint) * indices.size(),
               &indices.at(0), GL_DYNAMIC_DRAW);

  glDrawElements(mode, count, GL_UNSIGNED_INT, 0);

  // Add buffer index to free list to delete this buffer in later phase(at
  // glFinish())
  bufferFreeList_.push_back(idx);
}

void Context::glDrawElements(GLenum mode, GLsizei count, GLenum type,
                             const GLvoid *indices) {
  TRACE_EVENT("(GLenum mode = %d, GLsizei count = %d, GLenum type = %d, const "
              "GLvoid* indices = %p",
              mode, count, type, indices);

  if (count < 0) {
    return SetGLError(GL_INVALID_VALUE);
  } else if (count == 0) {
    return;
  }

  if ((mode != GL_TRIANGLES) && (mode != GL_POINTS) && (mode != GL_LINES) &&
      (mode != GL_TETRAHEDRONS_EXT) &&
      (mode != GL_PYRAMIDS_EXT) &&
      (mode != GL_PRISMS_EXT) &&
      (mode != GL_HEXAHEDRONS_EXT)) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // assume int index
  if ((type != GL_UNSIGNED_INT)) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // lookup the currently bound element buffers
  Buffer *elembuf = HandleToBuffer(GL_ELEMENT_ARRAY_BUFFER);
  if (elembuf == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  int k = state_.currentDrawStackIndex;
  assert(k < kMaxDrawStack);

  // See "position" attribute.
  Buffer *arraybuf = HandleToBuffer(GL_ARRAY_BUFFER);
  if (arraybuf == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }
  const VertexAttribute *attrPos =
      &state_.vertexAttributes[k][kVtxAttrPosition];
  if (attrPos->enabled == false) {
    return SetGLError(GL_INVALID_OPERATION);
  }
  Buffer *posbuf = resourceManager_.GetBuffer(attrPos->handle);
  if (posbuf == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // lookup the current shader program
  const Program *prg = resourceManager_.GetProgram(state_.currentProgram);
  if (prg == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  AccelBuilder::PrimitiveType primTy;
  if (mode == GL_TRIANGLES) {
    primTy = AccelBuilder::PRIMITIVE_TRIANGLES;
  } else if (mode == GL_POINTS) {
    primTy = AccelBuilder::PRIMITIVE_POINTS;
  } else if (mode == GL_LINES) {
    primTy = AccelBuilder::PRIMITIVE_LINES;
  } else if (mode == GL_TETRAHEDRONS_EXT) {
    primTy = AccelBuilder::PRIMITIVE_TETRAHEDRONS;
  } else if (mode == GL_PYRAMIDS_EXT) {
    primTy = AccelBuilder::PRIMITIVE_PYRAMIDS;
  } else if (mode == GL_PRISMS_EXT) {
    primTy = AccelBuilder::PRIMITIVE_PRISMS;
  } else if (mode == GL_HEXAHEDRONS_EXT) {
    primTy = AccelBuilder::PRIMITIVE_HEXAHEDRONS;
  } else {
    assert(0);
  }

  double bmin[3];
  double bmax[3];

  bool isDoublePrecisionPos = false;
  if (attrPos->type == GL_DOUBLE) {
    isDoublePrecisionPos = true;
  }

  unsigned char *accel = NULL;
  if (mode == GL_TRIANGLES) {
    AccelBuilder::TriangleAccelerator *meshAccel =
        accelBuilder_.BuildTriangleAccel(elembuf, posbuf, isDoublePrecisionPos,
                                         &state_.vertexAttributes[k].at(0),
                                         state_.texture2D, count,
                                         (GLubyte *)indices - (GLubyte *)NULL);
    assert(meshAccel);
    meshAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(meshAccel);
  } else if (mode == GL_POINTS) {

    float pointSize = 1.0f;
    GLuint psizePos = prg->GetUniformLocation("lsgl_PointSize");
    if (psizePos != (GLuint)(-1)) {
      const Uniform *psizeUniform = prg->GetUniform(psizePos);
      assert(psizeUniform->data.size() == sizeof(float));
      memcpy(&pointSize, &psizeUniform->data.at(0), sizeof(float));
    }

    // Find a vertex attribute of point size if avaiable.
    const float *pointSizeV = NULL;
    GLsizei pointSizeVLen = 0;

    GLuint pointSizeVPos = prg->GetVaryingLocation("lsgl_PointSize");
    const VertexAttribute *pointSizeAttr =
        &state_.vertexAttributes[k][pointSizeVPos];

    if (pointSizeAttr) {

      if ((pointSizeAttr->size == 1) && (attrPos->type == GL_FLOAT) &&
          (pointSizeAttr->normalized == GL_FALSE) &&
          (pointSizeAttr->stride == (sizeof(float))) &&
          (pointSizeAttr->enabled)) {

        Buffer *buf = resourceManager_.GetBuffer(pointSizeAttr->handle);
        if (buf) {
          pointSizeVLen = buf->GetSize() / sizeof(float);
          pointSizeV = reinterpret_cast<const float *>(buf->GetData());
        }
      }
    }

    AccelBuilder::ParticleAccelerator *particleAccel =
        accelBuilder_.BuildParticleAccel(elembuf, posbuf, isDoublePrecisionPos,
                                         &state_.vertexAttributes[k].at(0),
                                         count,
                                         (GLubyte *)indices - (GLubyte *)NULL,
                                         pointSize, pointSizeV, pointSizeVLen);
    particleAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(particleAccel);
  } else if (mode == GL_LINES) {

    float lineSize =
        state_.lineWidth; // use glLineWith() to set uniform line width.

    bool cap = false;
    GLuint capPos = prg->GetUniformLocation("lsgl_LineCap");
    if (capPos != (GLuint)(-1)) {
      const Uniform *capUniform = prg->GetUniform(capPos);
      assert(capUniform->data.size() == sizeof(int));
      int value = 0;
      memcpy(&value, &capUniform->data.at(0), sizeof(int));
      if (value > 0) {
        cap = true;
      }
    }

    // Find a vertex attribute of line size if avaiable.
    const float *lineSizeV = NULL;
    GLsizei lineSizeVLen = 0;

    GLuint lineSizeVPos = prg->GetVaryingLocation("lsgl_LineWidth");
    const VertexAttribute *lineSizeAttr =
        &state_.vertexAttributes[k][lineSizeVPos];

    if (lineSizeAttr) {

      if ((lineSizeAttr->size == 1) && (attrPos->type == GL_FLOAT) &&
          (lineSizeAttr->normalized == GL_FALSE) &&
          (lineSizeAttr->stride == (sizeof(float))) &&
          (lineSizeAttr->enabled)) {

        Buffer *buf = resourceManager_.GetBuffer(lineSizeAttr->handle);
        if (buf) {
          lineSizeVLen = buf->GetSize() / sizeof(float);
          lineSizeV = reinterpret_cast<const float *>(buf->GetData());
        }
      }
    }

    AccelBuilder::LineAccelerator *lineAccel =
        accelBuilder_.BuildLineAccel(elembuf, posbuf, isDoublePrecisionPos,
                                     &state_.vertexAttributes[k].at(0), count,
                                     (GLubyte *)indices - (GLubyte *)NULL,
                                     lineSize, lineSizeV, lineSizeVLen, cap);
    lineAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(lineAccel);

  } else if (mode == GL_TETRAHEDRONS_EXT) {
    AccelBuilder::TetraAccelerator *tetraAccel =
        accelBuilder_.BuildTetraAccel(elembuf, posbuf, isDoublePrecisionPos,
                                      &state_.vertexAttributes[k].at(0), count,
                                      (GLubyte *)indices - (GLubyte *)NULL);
    assert(tetraAccel);
    tetraAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(tetraAccel);
  } else if ((mode == GL_PRISMS_EXT) || (mode == GL_PYRAMIDS_EXT) || (mode == GL_HEXAHEDRONS_EXT)) {
    AccelBuilder::SolidAccelerator *solidAccel =
        accelBuilder_.BuildSolidAccel(mode, elembuf, posbuf, isDoublePrecisionPos,
                                      &state_.vertexAttributes[k].at(0), count,
                                      (GLubyte *)indices - (GLubyte *)NULL);
    assert(solidAccel);
    solidAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(solidAccel);
  } else {
    assert(0 && "Unsupported primitive type");
  }

  if (accel != NULL) {

    //
    // @todo { check if draw can be cached. }
    // STATIC && same size && same pointer = Can use cache.
    //
    // assert(elembuf->GetUsage() == GL_STATIC_DRAW);

    // Grab local -> world transform from predefined Uniform variable.
    // NOTE: bmin, bmax is defined in local space.
    //
    float Mw[4][4];
    float Mp[4][4];
    float Mv[4][4];

    std::string world = "lsgl_World";
    int locWorld = prg->GetUniformLocation(world);
    assert(locWorld != -1);
    Uniform const *uniformWorld = prg->GetUniform(locWorld);
    assert(uniformWorld->data.size() == 64); // sizeof(float)*4*4
    memcpy(Mw, &uniformWorld->data.at(0), sizeof(float) * 4 * 4);

    std::string proj = "lsgl_Proj";
    int locProj = prg->GetUniformLocation(proj);
    assert(locProj != -1);
    Uniform const *uniformProj = prg->GetUniform(locProj);
    assert(uniformProj->data.size() == 64);
    memcpy(Mp, &uniformProj->data.at(0), sizeof(float) * 4 * 4);

    std::string view = "lsgl_View";
    int locView = prg->GetUniformLocation(view);
    assert(locView != -1);
    Uniform const *uniformView = prg->GetUniform(locView);
    assert(uniformView->data.size() == 64);
    memcpy(Mv, &uniformView->data.at(0), sizeof(float) * 4 * 4);

    double Tw[4][4];
    // double Tp[4][4];
    // double Tv[4][4];

    // float -> double upcast
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        Tw[j][i] = Mw[j][i];
        // Tp[j][i] = Mp[j][i];
        // Tv[j][i] = Mv[j][i];
      }
    }

    //
    // Copy program state
    //
    Program *new_prg = new Program(*prg);

    // Prepare shader state
    FragmentState fragmentState;
    ShadingState shadingState;
    new_prg->PrepareEval(fragmentState, shadingState,
                         state_.vertexAttributes[state_.currentDrawStackIndex],
                         (*this));

    RenderElement renderElement(primTy, mode, accel, bmin, bmax, Tw,
                                state_.currentDrawStackIndex, new_prg,
                                fragmentState, shadingState);
    assert(renderGraph_);
    renderGraph_->AddRenderElement(renderElement);

    state_.currentDrawStackIndex++;

    //
    // Copy vertex attribute state
    //
    int dst = state_.currentDrawStackIndex;
    int src = dst - 1;
    assert(dst < kMaxDrawStack);
    assert(src >= 0);
    for (int indx = 0; indx < kMaxVertexAttribs; indx++) {
      state_.vertexAttributes[dst][indx] = state_.vertexAttributes[src][indx];
    }

    // and mark our frame as dirty
    dirtyFrame_ = true;
  }
}
