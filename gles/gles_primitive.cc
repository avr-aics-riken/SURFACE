/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
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

  if ((mode != GL_TRIANGLES) && (mode != GL_POINTS) && (mode != GL_LINES) && (mode != GL_TETRAHEDRONS_EXT)) {
    return SetGLError(GL_INVALID_ENUM);
  }

  // lookup the currently bound array buffers
  Buffer *arraybuf = HandleToBuffer(GL_ARRAY_BUFFER);
  if (arraybuf == NULL) {
    return SetGLError(GL_INVALID_OPERATION);
  }

  // Construct index buffer
  std::vector<GLint> indices;

  if (mode == GL_TRIANGLES) {
    //   triangle
    for (GLint i = 0; i < count; i++) {
      indices.push_back(3 * (i + first) + 0);
      indices.push_back(3 * (i + first) + 1);
      indices.push_back(3 * (i + first) + 2);
    }
  } else if (mode == GL_POINTS) {
    //   points
    for (GLint i = 0; i < count; i++) {
      indices.push_back(i + first);
    }
  } else if (mode == GL_LINES) {
    //   lines
    for (GLint i = 0; i < count; i++) {
      indices.push_back(2 * (i + first) + 0);
      indices.push_back(2 * (i + first) + 1);
    }
  } else if (mode == GL_TETRAHEDRONS_EXT) {
    //   tetrahedronss
    for (GLint i = 0; i < count; i++) {
      indices.push_back(4 * (i + first) + 0);
      indices.push_back(4 * (i + first) + 1);
      indices.push_back(4 * (i + first) + 2);
      indices.push_back(4 * (i + first) + 3);
    }
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

  // assume triangle or point input.
  if ((mode != GL_TRIANGLES) && (mode != GL_POINTS) && (mode != GL_LINES) && (mode != GL_TETRAHEDRONS_EXT)) {
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

  AccelBuilder::PrimitiveType primTy;
  if (mode == GL_TRIANGLES) {
    primTy = AccelBuilder::PRIMITIVE_TRIANGLES;
  } else if (mode == GL_POINTS) {
    primTy = AccelBuilder::PRIMITIVE_POINTS;
  } else if (mode == GL_LINES) {
    primTy = AccelBuilder::PRIMITIVE_LINES;
  } else if (mode == GL_TETRAHEDRONS_EXT) {
    primTy = AccelBuilder::PRIMITIVE_TETRAHEDRONS;
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
    AccelBuilder::MeshAccelerator *meshAccel =
        meshBuilder_.Build(elembuf, posbuf, isDoublePrecisionPos,
                           &state_.vertexAttributes[k].at(0), state_.texture2D,
                           count, (GLubyte *)indices - (GLubyte *)NULL);
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
        meshBuilder_.BuildParticleAccel(elembuf, posbuf, isDoublePrecisionPos,
                                        &state_.vertexAttributes[k].at(0),
                                        count,
                                        (GLubyte *)indices - (GLubyte *)NULL,
                                        pointSize, pointSizeV, pointSizeVLen);
    particleAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(particleAccel);
  } else if (mode == GL_LINES) {
    AccelBuilder::LineAccelerator *lineAccel = meshBuilder_.BuildLineAccel(
        elembuf, posbuf, isDoublePrecisionPos,
        &state_.vertexAttributes[k].at(0), count,
        (GLubyte *)indices - (GLubyte *)NULL, state_.lineWidth);
    lineAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(lineAccel);
  } else if (mode == GL_TETRAHEDRONS_EXT) {
    AccelBuilder::TetraAccelerator *tetraAccel =
        meshBuilder_.BuildTetraAccel(elembuf, posbuf, isDoublePrecisionPos,
                           &state_.vertexAttributes[k].at(0),
                           count, (GLubyte *)indices - (GLubyte *)NULL);
    assert(tetraAccel);
    tetraAccel->BoundingBox(bmin, bmax);
    accel = reinterpret_cast<unsigned char *>(tetraAccel);
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
    double Tp[4][4];
    double Tv[4][4];

    // float -> double upcast
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        Tw[j][i] = Mw[j][i];
        Tp[j][i] = Mp[j][i];
        Tv[j][i] = Mv[j][i];
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

    RenderElement renderElement(primTy, accel, bmin, bmax, Tw,
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
