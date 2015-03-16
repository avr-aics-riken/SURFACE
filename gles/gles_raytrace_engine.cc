/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "GLES2/gl2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LSGL_ENABLE_K_PROFILE
#include <fj_tool/fipp.h>
#include <fjcoll.h>
#endif

#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cassert>
#include <algorithm>

#include "gles_common.h"
#include "gles_context.h"

#include "../render/render_timerutil.h"

// 1 if measure pure raytracing performance(disable shader call)
#define ENABLE_RAYTRACE_PERF (0)

#ifdef __sparc__
// Disable cancel mechanism on K/FX10 for a while
#define LSGL_ENABLE_RENDER_CANCEL (0)
#else
#define LSGL_ENABLE_RENDER_CANCEL (1)
#endif

using namespace lsgl;

RaytraceEngine *RaytraceEngine::sRaytraceEngine = NULL;

namespace {

static inline unsigned char clamp(float f) {
  // f = pow(f, 1.0/2.2);    // gamma
  int i = (int)(f * 255.0);
  if (i < 0)
    i = 0;
  if (i > 255)
    i = 255;
  return (unsigned char)i;
}

typedef void (*WriteFunction)(char *ptr, float r, float g, float b, float a);

static inline void WriteColorRGBAU8(char *ptr, float r, float g, float b,
                                    float a) {
  unsigned char *out = reinterpret_cast<unsigned char *>(ptr);
  out[0] = clamp(r);
  out[1] = clamp(g);
  out[2] = clamp(b);
  out[3] = clamp(a);
}

static inline void WriteColorRGBAF32(char *ptr, float r, float g, float b,
                                     float a) {
  float *out = reinterpret_cast<float *>(ptr);
  out[0] = r;
  out[1] = g;
  out[2] = b;
  out[3] = a;
}

//
// Screen division utils
//
#if defined(LSGL_ENABLE_MPI) && defined(LSGL_ENABLE_SCREEN_PARALLEL)

// inline unsigned int CalcAreasize(const unsigned long size_u,
//                                 const unsigned long size_v) {
//  return (size_u * size_v);
//}

inline double CalcAspect(const unsigned long size_u,
                         const unsigned long size_v) {
  return (static_cast<double>(size_u) / static_cast<double>(size_v));
}

inline double EvalAspect(const double aspect) { return (fabs(1.0 - aspect)); }

bool DecideDivPattern(const unsigned int numDiv, const unsigned int size_u,
                      const unsigned int size_v, unsigned int *div_u,
                      unsigned int *div_v) {
  if (!div_u || !div_v) {
    return false;
  }

  if (numDiv <= 1) {
    *div_u = *div_v = 1;
    return false;
  }

  double minValue = std::numeric_limits<double>::max();

  unsigned long div = numDiv;
  unsigned long size[2] = {static_cast<unsigned long>(size_u),
                           static_cast<unsigned long>(size_v)};

  unsigned long divPttn[2] = {0, 0};

  for (unsigned long i = 1; i <= div; i++) {
    if (div % i != 0)
      continue;
    unsigned long jmax = div / i;
    for (unsigned long j = 1; j <= jmax; j++) {
      if ((div / i) % j != 0)
        continue;
      if ((size[1] / j) < 1)
        break;
      if (i * j == div) {
        double value = EvalAspect(CalcAspect((size[0] / i), (size[1] / j)));
        // printf("value : %lf [i:%ld] [j:%ld]\n", value, i, j);

        if (value < minValue) {
          minValue = value;
          divPttn[0] = i;
          divPttn[1] = j;
        }
      }
    }
  }

  if (divPttn[0] == 0 || divPttn[1] == 0) {
    divPttn[0] = divPttn[1] = 1;
    if (size_u > size_v) {
      if (size_u < numDiv) {
        *div_u = *div_v = 1;
        return false;
      }
      divPttn[0] = numDiv;
    } else {
      if (size_v < numDiv) {
        *div_u = *div_v = 1;
        return false;
      }
      divPttn[0] = numDiv;
    }
  }

  *div_u = static_cast<unsigned int>(divPttn[0]);
  *div_v = static_cast<unsigned int>(divPttn[1]);

  return true;
}

// For given screen size and N,
// find best 2D division of the screen with N items, and also compute screen
// region of the given rank.
void StaticallyDivideScreen(int region[4], // left, upper, right, bottom
                            int width, int height, int N, int rank) {
  unsigned int divPattern[2];
  if (!DecideDivPattern(N, width, height, &divPattern[0], &divPattern[1])) {
    fprintf(stderr, "[LSGL] Faild to divide\n");
    region[0] = 0;
    region[1] = 0;
    region[2] = 0;
    region[3] = 0;
    assert(0);
    return;
  }

  printf("[LSGL] dbg: div pattern = %d, %d\n", divPattern[0], divPattern[1]);

  unsigned int tileW =
      std::max((unsigned int)1, (unsigned int)(width / divPattern[0]));
  unsigned int tileH =
      std::max((unsigned int)1, (unsigned int)(height / divPattern[1]));

  printf("[LSGL] dbg: tile size = %d, %d\n", tileW, tileH);

  // Find rank's 2D position
  int rankX = rank % divPattern[0];
  int rankY = rank / divPattern[0];

  printf("[LSGL] dbg: rank = %d, xy = (%d, %d)\n", rank, rankX, rankY);

  region[0] = rankX * tileW;
  region[1] = rankY * tileW;
  if (rankX == divPattern[0]) {
    region[2] = width; // include remainder
  } else {
    region[2] = (rankX + 1) * tileW;
  }
  if (rankY == divPattern[1]) {
    region[3] = height; // include remainder
  } else {
    region[3] = (rankY + 1) * tileH;
  }
  printf("[LSGL] dbg: rank = %d, region = (%d, %d), (%d, %d)\n", rank,
         region[0], region[1], region[2], region[3]);

  return;
}
#endif

void CalculateRenderRegion(int region[4], // left, upper, right, bottom
                           int width, int height) {
// Automatic screen-tiling rendering is valid if MPI and SCREEN_PARALLEL flag
// was set.
#if defined(LSGL_ENABLE_MPI) && defined(LSGL_ENABLE_SCREEN_PARALLEL)
  //
  // @todo { Dynamic division }
  //
  int numNodes;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);

  StaticallyDivideScreen(region, width, height, numNodes, rank);
#else
  region[0] = 0;
  region[1] = 0;
  region[2] = width;
  region[3] = height;
#endif
}

#if defined(LSGL_ENABLE_MPI) && defined(LSGL_ENABLE_SCREEN_PARALLEL)
template <typename T>
void MergeScreen(int region[4], int masterRank, int rank,
                 unsigned char *framebuffer, // [inout]
                 int width, int height) {
  int ret;

  //
  // !!! Assume all subregion(rank) has same extent. !!!
  //

  int subWidth = region[2] - region[0];
  int subHeight = region[3] - region[1];

  std::vector<T> buf(subWidth * subHeight * 4);
  T *ptr = reinterpret_cast<T *>(framebuffer);

  int k = 0;
  for (int y = region[1]; y < region[3]; y++) {
    for (int x = region[0]; x < region[2]; x++) {
      int srcIdx = (height - y - 1) * width + x;
      buf[4 * k + 0] = ptr[4 * srcIdx + 0];
      buf[4 * k + 1] = ptr[4 * srcIdx + 1];
      buf[4 * k + 2] = ptr[4 * srcIdx + 2];
      buf[4 * k + 3] = ptr[4 * srcIdx + 3];
      k++;
    }
  }

  //
  // Use MPI gather to collect render subregion from each rank.
  //
  MPI_Barrier(MPI_COMM_WORLD);

  int size = subWidth * subHeight * 4 * sizeof(T);
  ret = MPI_Gather(&buf.at(0), size, MPI_BYTE, framebuffer, size, MPI_BYTE,
                   masterRank, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  MPI_Barrier(MPI_COMM_WORLD);

  // NOTE: Don't reorder framebuffer to save memory.
  // Application must reorder framebuffer.
}
#endif

} // namespace

RaytraceEngine::RaytraceEngine()
    : framebuffer_(NULL), numRays_(0), pixelStep_(1),
      progressCallbackFunc_(NULL), callbackUserData_(NULL) {
  sRaytraceEngine = this;

  // Default = orthographic camera, positioned at (0.0, 0.0, 0.0) and see (0.0,
  // 0.0, -1.0);
  GLfloat eye[3] = {0.0f, 0.0f, 0.0f};
  GLfloat target[3] = {0.0f, 0.0f, -1.0f};
  GLfloat up[3] = {0.0f, 1.0f, -1.0f};
  GLfloat fov = 0.0f; // = orthographic
  SetCamera(eye, target, up, fov);
}

RaytraceEngine::~RaytraceEngine() { sRaytraceEngine = NULL; }

void RaytraceEngine::SetCamera(const GLfloat *eye, const GLfloat *target,
                               const GLfloat *up, GLfloat fov) {
  if (camera_) {
    delete camera_;
  }

  real3 veye(eye[0], eye[1], eye[2]);
  real3 vtarget(target[0], target[1], target[2]);
  real3 vup(up[0], up[1], up[2]);

  // If fov <= 0, assume camera is orthographic.
  bool ortho = false;
  float eps = std::numeric_limits<float>::epsilon();
  if (fov < eps) {
    ortho = true;
  }

  camera_ = new Camera(veye, vtarget, vup, fov, ortho);
}

void RaytraceEngine::SetStereoEnvCamera(const GLfloat *eye,
                                        const GLfloat *target,
                                        const GLfloat *up, GLfloat zeroParallax,
                                        GLfloat eyeSeparation) {
  if (camera_) {
    delete camera_;
  }

  real3 veye(eye[0], eye[1], eye[2]);
  real3 vtarget(target[0], target[1], target[2]);
  real3 vup(up[0], up[1], up[2]);

  bool ortho = false;

  camera_ = new Camera(veye, vtarget, vup, zeroParallax, eyeSeparation);
}

void RaytraceEngine::OnPrepare(const Context *ctx) {
  // simply store pointer to context
  ctx_ = ctx;
}

bool RaytraceEngine::OnStart(Framebuffer *fb, RenderGraph *renderGraph) {
  //
  // @todo { Call camera shader here. }
  //
  assert(camera_);

  // ensure framebuffer is valid
  if ((fb == NULL) || (fb->HasValidSize() == false)) {
    return false;
  }

  // Don't forget to initialize RNGs.
  initialize_random();

  framebuffer_ = fb;
  accumTime_ = 0.0;

  renderGraph_ = renderGraph;

  return true;
}

// C wrapper function to support callback from a shader
bool scene_trace(Intersection &isect, Ray &r) {
  return RaytraceEngine::GetRaytraceEngine()->Trace(isect, r);
}

bool RaytraceEngine::Trace(Intersection &isectRet, Ray &r) {

  Intersection isect;
  isect.clear();
  isect.t = std::numeric_limits<float>::max();
  float tMax = std::numeric_limits<float>::max();

  if ((!renderGraph_) || !renderGraph_->IsBuilt()) {
    assert(0);
    return false;
  }

  // Trace into render graph.
  bool ret = renderGraph_->Trace(isect, r);
  if (ret) {
    isectRet = isect;
  }

  return ret;
}

void RaytraceEngine::ShadeFragment(float shadecol[4], bool &discarded /* out */,
                                   const Intersection &isect, const Ray &ray,
                                   int thread_id) {

  discarded = false;

  if (!isect.renderElement) {
    assert(0);
    return;
  }

  RenderElement *renderElement =
      reinterpret_cast<RenderElement *>(isect.renderElement);

  const Program *prg = renderElement->GetProgram();
  if (!prg) {
    return;
  }

  if (prg->IsLinked() == true) {

    GLfloat fragCol[4] = {0.0f, 0.0f, 0.0f, 0.0f};

    // fragCoord[0] = 0.5f + px;
    // fragCoord[1] = 0.5f + py;
    // fragCoord[2] = isect.t;     // @fixme
    // fragCoord[3] = 1.0f / isect.t;
    //

    FragmentState &fragmentState = renderElement->GetFragmentState();
    ShadingState &shadingState = renderElement->GetShadingState();

    GLfloat fragCoord[4];
    fragCoord[0] = ray.px;
    fragCoord[1] = ray.py;
    fragCoord[2] = isect.t; // @fixme
    fragCoord[3] = 1.0f / isect.t;

    FragmentShader::CameraInfo cameraInfo;
    cameraInfo.frame[0][0] = camera_->getEye()[0];
    cameraInfo.frame[0][1] = camera_->getEye()[1];
    cameraInfo.frame[0][2] = camera_->getEye()[2];
    cameraInfo.frame[1][0] = camera_->getLookat()[0];
    cameraInfo.frame[1][1] = camera_->getLookat()[1];
    cameraInfo.frame[1][2] = camera_->getLookat()[2];
    cameraInfo.frame[2][0] = camera_->getUp()[0];
    cameraInfo.frame[2][1] = camera_->getUp()[1];
    cameraInfo.frame[2][2] = camera_->getUp()[2];
    cameraInfo.fov = camera_->getFov();

    const std::vector<VertexAttribute> &vertexAttributes =
        ctx_->state_
            .vertexAttributes[renderElement->GetDrawStackIndex()]; // @fixme

    assert(prg->GetFragmentShader(0));

    IntersectionState isectState;
    isectState.position = isect.position;
    isectState.normal = isect.normal;
    isectState.geometricNormal = isect.geometric;
    isectState.tangent = isect.tangent;
    isectState.binormal = isect.binormal;
    isectState.raydir = ray.direction();
    isectState.raydepth = ray.depth;
    isectState.px = ray.px;
    isectState.py = ray.py;
    isectState.doubleSided = ray.double_sided;
    isectState.rayattrib = ray.user_attrib;
    isectState.prev_node = isect.renderElement;
    isectState.prev_prim_id = isect.prim_id;
    isectState.u = isect.u;
    isectState.v = isect.v;
    isectState.f0 = isect.f0;
    isectState.f1 = isect.f1;
    isectState.f2 = isect.f2;

    bool ret = prg->GetFragmentShader(0)->Eval(
        fragCol, fragmentState, shadingState, vertexAttributes, fragCoord,
        isectState, cameraInfo, thread_id);

    if (ret) {
      shadecol[0] = fragCol[0];
      shadecol[1] = fragCol[1];
      shadecol[2] = fragCol[2];
      shadecol[3] = fragCol[3];
    } else {
      discarded = true;
    }

  } else {
    // fall back to show_normal shader.
    shadecol[0] = 0.5f * isect.geometric[0] + 0.5f;
    shadecol[1] = 0.5f * isect.geometric[1] + 0.5f;
    shadecol[2] = 0.5f * isect.geometric[2] + 0.5f;
    shadecol[3] = 1.0f;
  }
}

namespace {

inline void ViewportConversion(int &xw, int &yw, int x, int y, int width,
                               int height, float x01, // [0, 1]
                               float y01)             // [0, 1]
{
  // See glViewport
  // http://www.opengl.org/sdk/docs/man/xhtml/glViewport.xml

  // First convert to normalized device coord, [-1, 1]^2
  float xnd = 2.0f * (x01 - 0.5f);
  float ynd = 2.0f * (y01 - 0.5f);

  // flip sign of x since we use raytracing camera.
  xw = (xnd + 1.0f) * (0.5f * width) - x;
  yw = (ynd + 1.0f) * (0.5f * height) + y;
}
}

bool RaytraceEngine::OnData() {
  real3 corner, udir, vdir;
  const int subSamples = (ctx_->state_.sampleCoverage == true)
                             ? ctx_->state_.sampleCoverageValue
                             : 1;
  const int subSamplesSqr = subSamples * subSamples;
  const float invN = 1.0f / (float)subSamplesSqr;

  // HACK: assume framebuffer has a valid render color buffer at index zero
  const int width = framebuffer_->GetWidth();
  const int height = framebuffer_->GetHeight();

  WriteFunction writeFunc = NULL;
  char *colorBufferPtr = NULL;
  int bytesPerPixel = 0;

  if (framebuffer_->GetColorBuffer(0)) {

    colorBufferPtr = framebuffer_->GetColorBuffer(0)->GetBuffer();
    bytesPerPixel = framebuffer_->GetColorBuffer(0)->GetBytesPerPixel();

    // Find suitable pixel write function.
    if ((framebuffer_->GetColorBuffer(0)->GetFormat() == GL_RGBA8_OES) ||
        (framebuffer_->GetColorBuffer(0)->GetFormat() == GL_RGBA)) {
      writeFunc = WriteColorRGBAU8;
      assert(bytesPerPixel == 4);
    } else if ((framebuffer_->GetColorBuffer(0)->GetFormat() ==
                GL_RGBA32F_EXT)) {
      writeFunc = WriteColorRGBAF32;
      assert(bytesPerPixel == 16);
    }
  }

  // Assume float buffer.
  float *depth = NULL;

  if (framebuffer_->GetDepthBuffer()) {
    if ((framebuffer_->GetDepthBuffer()->GetFormat() ==
         GL_DEPTH_COMPONENT32_OES) ||
        (framebuffer_->GetDepthBuffer()->GetFormat() == GL_DEPTH_COMPONENT)) {
      depth = reinterpret_cast<float *>(
          framebuffer_->GetDepthBuffer()->GetBuffer());
    }
  }

  // render region
  int startX = 0, startY = 0;
  int endX = width, endY = height;

  // Override render region if scissor test was set.
  if (ctx_->state_.scissorTest) {
    startX = std::max(0, ctx_->state_.scissorX);
    startY = std::max(0, ctx_->state_.scissorY);
    endX = startX + std::max(0, ctx_->state_.scissorWidth);
    endY = startY + std::max(0, ctx_->state_.scissorHeight);
  }

  printf("[LSGL] DBG: rect (%d, %d) - (%d, %d)\n", startX, startY, endX, endY);

  // Calculate actual render region this (compute) node should consider.
  int region[4];
  CalculateRenderRegion(region, (endX - startX), (endY - startY));
  region[0] += startX;
  region[1] += startY;
  region[2] += startX;
  region[3] += startY;

#ifdef LSGL_ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("[LSGL] DBG: rank[%d] region (%d, %d) - (%d, %d)\n", rank, region[0],
         region[1], region[2], region[3]);
#else
  printf("[LSGL] DBG: region (%d, %d) - (%d, %d)\n", region[0], region[1],
         region[2], region[3]);
#endif

#ifdef LSGL_ENABLE_K_PROFILE
  // fipp_start();
  start_collection("render1");
#endif

  renderTimer_.start();

  {
    camera_->buildFrame(corner, udir, vdir, width, height);
    // printf("corner = %f, %f, %f\n", corner[0], corner[1], corner[2]);
    // printf("udir  = %f, %f, %f\n", udir[0], udir[1], udir[2]);
    // printf("vdir  = %f, %f, %f\n", vdir[0], vdir[1], vdir[2]);
    real3 look = cross(udir, vdir);
    look.normalize();

    const float invWidth = 1.0f / width;
    const float invHeight = 1.0f / height;

    int step = this->pixelStep_;

    int numLines = (region[3] - region[1]) / step;
    if (numLines < 1)
      numLines = 1;
    int showCountTick = numLines / 100;
    if (showCountTick < 1)
      showCountTick = 1;

#if LSGL_ENABLE_RENDER_CANCEL
    int progress = 0;
    bool canceled = false;
    bool hasCanceled = false;
#endif

//#ifdef _OPENMP
// printf("[OMP] maxthreads = %d\n", omp_get_max_threads());
// printf("[OMP] numthreads = %d\n", omp_get_num_threads());
//#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int y = region[1]; y < region[3]; y += step) {

#if LSGL_ENABLE_RENDER_CANCEL
#ifdef _OPENMP
#pragma omp critical
      hasCanceled = canceled;

      if (hasCanceled) {
        continue;
      }
#else
      if (canceled)
        break;
#endif
#endif

#if LSGL_ENABLE_RENDER_CANCEL
#ifdef _OPENMP
#pragma omp atomic
      progress++;
#else
      progress++;
#endif

#ifdef _OPENMP
      if (omp_get_thread_num() == 0) { // master thread
        if (progress % showCountTick == 0) {
          if (progressCallbackFunc_) {
            bool ok = progressCallbackFunc_(
                (int)(100.0 * progress / (double)numLines), y, height, callbackUserData_);

            if (!ok) {
#pragma omp critical
              canceled = true;
            }
          } else {
            printf("\r[LSGL] render %d / %d", y, height);
          }
        }
      }
#else
      if (progress % showCountTick == 0) {
        if (progressCallbackFunc_) {
          bool ok = progressCallbackFunc_(
              (int)(100.0 * progress / (double)numLines), y, height, callbackUserData_);
          if (!ok) {
            canceled = true;
          }
        } else {
          printf("\r[LSGL] render %d / %d", y, height);
        }
      }
#endif
#endif

      int thread_id = 0;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#endif

      for (int x = region[0]; x < region[2]; x += step) {
        float fcol[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float depthval = 0.0f;
        bool fragment = false;

        //
        // First do viewport conversion, and if the pixel doesn't fit
        // into screen, just do nothing for this pixel.
        //
        float x01 = (x + 0.5f) * invWidth;
        float y01 = (y + 0.5f) * invHeight;
        int px = -1, py = -1;
        ViewportConversion(px, py, ctx_->state_.viewportX,
                           ctx_->state_.viewportY, ctx_->state_.viewportWidth,
                           ctx_->state_.viewportHeight, x01, y01);

        if ((px < 0) || (py < 0) || (px >= width) || (py >= height)) {
          continue;
        }

        for (int vv = 0; vv < subSamples; vv++) {
          for (int uu = 0; uu < subSamples; uu++) {

            float ju = 0.0f;
            float jv = 0.0f;
            if (uu != 0)
              ju = randomreal(thread_id);
            if (vv != 0)
              jv = randomreal(thread_id);

            // get eye dir.
            float u = (uu + ju) / (float)subSamples;
            float v = (vv + jv) / (float)subSamples;

            real3 org;
            real3 dir;

#if 0 // Flip Y
            if (camera_->isOrtho()) {
              org = (corner + (x + u + step / 2.0f) * udir +
                     ((height - y - 1) - v - step / 2.0f) * vdir);
              dir = look;
            } else {
              org = camera_->getEye();
              dir = (corner + (x + u + step / 2.0f) * udir +
                     ((height - y - 1) - v - step / 2.0f) * vdir) -
                    camera_->getEye();
            }
#else

            if (camera_->isStereoEnv()) {
              camera_->generateStereoEnvRay(
                  org, dir, (float)px + u + step / 2.0f,
                  (float)py + v - step / 2.0f + 1.0f, width, height);
              // printf("org = %f, %f, %f\n", org[0], org[1], org[2]);
              // printf("dir = %f, %f, %f\n", dir[0], dir[1], dir[2]);

            } else if (camera_->isOrtho()) {
              // Add +1.0f in Y to match OpenGL's rasterizer rule.
              org = (corner + ((float)px + u + step / 2.0f) * udir +
                     ((float)py + v - step / 2.0f + 1.0f) * vdir);
              dir = look;
            } else {
              org = camera_->getEye();
              // Add +1.0f in Y to match OpenGL's rasterizer rule.
              dir = (corner + ((float)px + u + step / 2.0f) * udir +
                     ((float)py + v - step / 2.0f + 1.0f) * vdir) -
                    org;
            }
#endif

            // normalize dir
            dir.normalize();
            Ray ray(org, dir);

            // Store pixel position to the ray for the use in the shader.
            ray.px = px + u;
            ray.py = height - py - 1 + v;
            ray.depth = 0; // eye ray has depth 0.
            ray.user_attrib = 0.0f;
            ray.double_sided = 1; // eye ray is always do intersect testing with
                                  // double-sided enabled.
            ray.prev_prim_id = (unsigned int)(-1);
            ray.prev_node = NULL;

            Intersection isect;

            bool ret = scene_trace(isect, ray);

            if (ret) {
              float shadecol[4] = {0.0f, 0.0f, 0.0f, 0.0f};
#if ENABLE_RAYTRACE_PERF
              shadecol[0] = 1.0f;
              shadecol[1] = 1.0f;
              shadecol[2] = 1.0f;
              shadecol[3] = 1.0f;

              fcol[0] += shadecol[0];
              fcol[1] += shadecol[1];
              fcol[2] += shadecol[2];
              fcol[3] += shadecol[3];

              depthval += isect.t;
              fragment = true;
#else
              bool discarded = false;
              ShadeFragment(shadecol, discarded, isect, ray, thread_id);
              if (!discarded) {
                fcol[0] += shadecol[0];
                fcol[1] += shadecol[1];
                fcol[2] += shadecol[2];
                fcol[3] += shadecol[3];

                // @fixme { Don't take averate of depth. }
                // @fixme { Compute GL-compatible depth value. }
                depthval += isect.t;
                fragment = true;
              }
#endif
            }
          }
        }

        // in GL, we only write pixels to the framebuffer which actually have
        // rasterized fragments (hit rays)
        if (fragment == true) {
          depthval *= invN; // take avarage of depth value.

          // Y-down
          // block
          for (int v = 0; v < step; v++) {
            if ((height - y - 1 + v) > (height - 1))
              continue;
            for (int u = 0; u < step; u++) {
              if (x + u > (width - 1))
                continue;
              if (writeFunc && colorBufferPtr) {
                writeFunc(colorBufferPtr +
                              (bytesPerPixel *
                               ((height - y - 1 + v) * width + (x + u))),
                          fcol[0] * invN, fcol[1] * invN, fcol[2] * invN,
                          fcol[3] * invN);
              }

              // @todo { get depth from a fragment shader output. }
              if (depth) {
                depth[(height - y - 1 + v) * width + (x + u)] = depthval;
              }
            }
          }
        }
      }
    }
  }

  renderTimer_.end();

#if defined(LSGL_ENABLE_MPI) && defined(LSGL_ENABLE_SCREEN_PARALLEL)

  // Merge color buffer
  if (colorBufferPtr) {
    timerutil mergeTimer;
    mergeTimer.start();

    int masterRank = 0;
    assert((bytesPerPixel == 4) || (bytesPerPixel == 16));
    if (bytesPerPixel == 4) { // assume byte RGBA
      MergeScreen<unsigned char>(region, masterRank, rank,
                                 (unsigned char *)colorBufferPtr, width,
                                 height);
    } else if (bytesPerPixel == 16) { // assume float RGBA
      MergeScreen<float>(region, masterRank, rank,
                         (unsigned char *)colorBufferPtr, width, height);
    }
    mergeTimer.end();

    if (rank == masterRank) {
      printf("\n[LSGL] Color buffer merge time: %d ms\n",
             (int)mergeTimer.msec());
      fflush(stdout);
    }
  }

  if (depth) {
    timerutil mergeTimer;
    mergeTimer.start();

    int masterRank = 0;
    // Assume 32bit float depth buffer.
    MergeScreen<float>(region, masterRank, rank,
                       (unsigned char *)colorBufferPtr, width, height);
    mergeTimer.end();

    if (rank == masterRank) {
      printf("\n[LSGL] Depth buffer merge time: %d ms\n",
             (int)mergeTimer.msec());
      fflush(stdout);
    }
  }
#endif

#ifdef LSGL_ENABLE_K_PROFILE
  // fipp_stop();
  stop_collection("render1");
#endif

  accumTime_ += renderTimer_.msec();

  int step = this->pixelStep_;
  numRays_ = subSamplesSqr * width * height / (double)(step * step);

  return true;
}

void RaytraceEngine::OnEnd() {
#ifdef LSGL_ENABLE_MPI
  int ret = MPI_Barrier(MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);
#endif
  // # of eye rays
  double mrays = numRays_ / (double)renderTimer_.msec();
  mrays /= (double)(1000);
  printf("\n[LSGL] render time: %d ms(%.4f Mraycasts)\n", (int)accumTime_,
         mrays);
  fflush(stdout);
}
