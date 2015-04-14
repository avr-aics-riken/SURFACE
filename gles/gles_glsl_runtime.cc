/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science,
 *RIKEN.
 * All rights reserved.
 *
 */

#include "../glsl/glsl_runtime.h"

#include "gles_context.h"
#include "gles_raytrace_engine.h"

using namespace lsgl;

//
// Implementation of GLSL Bulitin functions which interact with the renderer.
//

vec4 texture2D(unsigned long long samplerID, vec2 *coords) {
  const Context &ctx = Context::GetCurrentContext();

  vec4 col = {0.0f, 0.0f, 0.0f, 0.0f};
#ifdef LSGL_OPTIMIZE_GLSL
  printf("samplerID = %lld\n", samplerID);
  if (samplerID != 0) {
    float rgba[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    const Texture *tex = reinterpret_cast<const Texture *>(samplerID);
    if (tex) {
      tex->Fetch(rgba, coords->v[0], coords->v[1]);
      col.v[0] = rgba[0];
      col.v[1] = rgba[1];
      col.v[2] = rgba[2];
      col.v[3] = rgba[3];
    }
  }
#else
  if (ctx.resourceManager_.IsValidTexture(samplerID) == true) {
    float rgba[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    const Texture *tex = ctx.resourceManager_.GetTexture(samplerID);
    if (tex) {
      tex->Fetch(rgba, coords->v[0], coords->v[1]);
      col.v[0] = rgba[0];
      col.v[1] = rgba[1];
      col.v[2] = rgba[2];
      col.v[3] = rgba[3];
    }
  }
#endif

  return col;
}

vec4 texture3D(unsigned long long samplerID, vec3 *coords) {
  Context &ctx = Context::GetCurrentContext();

  vec4 col = {0.0f, 0.0f, 0.0f, 0.0f};
#ifdef LSGL_OPTIMIZE_GLSL
  if (samplerID != 0) {
    float rgba[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    Texture *tex = reinterpret_cast<Texture *>(samplerID);
    if (tex) {
      tex->Fetch(rgba, coords->v[0], coords->v[1], coords->v[2]);
      col.v[0] = rgba[0];
      col.v[1] = rgba[1];
      col.v[2] = rgba[2];
      col.v[3] = rgba[3];
    }
  }
#else
  if (ctx.resourceManager_.IsValidTexture(samplerID) == true) {
    float rgba[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    Texture *tex = ctx.resourceManager_.GetTexture(samplerID);
    if (tex) {
      tex->Fetch(rgba, coords->v[0], coords->v[1], coords->v[2]);
      col.v[0] = rgba[0];
      col.v[1] = rgba[1];
      col.v[2] = rgba[2];
      col.v[3] = rgba[3];
    }
  }
#endif

  return col;
}

float trace(void *fragptr, vec3 *org, vec3 *dir, vec4 *shadecol,
            float *userattrib) {
  Fragment *frag = reinterpret_cast<Fragment *>(fragptr);

  if (shadecol) {
    shadecol->v[0] = 0.0f;
    shadecol->v[1] = 0.0f;
    shadecol->v[2] = 0.0f;
    shadecol->v[3] = 0.0f;
  }

  float userval = 0.0f;
  if (userattrib) {
    userval = (*userattrib);
  }

  if (frag->raydepth > 60) {
    // too many reflection
    fprintf(stderr, "[LSGL] Too many ray reflections: %d\n", frag->raydepth);
    assert(0);
    return -1.0f;
  }

  //
  // Trace ray
  //
  Intersection isect;

  real3 rayorg;
  real3 raydir;

  rayorg[0] = org->v[0];
  rayorg[1] = org->v[1];
  rayorg[2] = org->v[2];

  raydir[0] = dir->v[0];
  raydir[1] = dir->v[1];
  raydir[2] = dir->v[2];

  Ray ray(rayorg, raydir);
  ray.depth = frag->raydepth + 1;
  ray.px = frag->px;
  ray.py = frag->py;
  ray.user_attrib = userval;
  ray.double_sided = frag->doubleSided;
  ray.prev_prim_id = frag->prev_prim_id;
  ray.prev_node = frag->prev_node;
  ray.prev_hit_t = frag->prev_hit_t;
  ray.prev_hit_normal = frag->prev_hit_normal;

  int thread_id = frag->threadID;

  bool hit = RaytraceEngine::GetRaytraceEngine()->Trace(isect, ray);

  //
  // Call fragment shader
  //
  if (hit) {
    float col[4] = {0.0, 0.0, 0.0, 0.0};
    bool discarded = false;

    RaytraceEngine::GetRaytraceEngine()->ShadeFragment(col, discarded, isect,
                                                       ray, thread_id);

    // @fixme { Handle discarded fragment correctly. }
    if (shadecol && !discarded) {
      shadecol->v[0] = col[0];
      shadecol->v[1] = col[1];
      shadecol->v[2] = col[2];
      shadecol->v[3] = col[3];
    }

    return isect.t;

  } else {
    // May hit environment, but fill with zero in this time.
    if (shadecol) {
      shadecol->v[0] = 0.0f;
      shadecol->v[1] = 0.0f;
      shadecol->v[2] = 0.0f;
      shadecol->v[3] = 1.0f;
    }

    return -(std::numeric_limits<float>::max)(); // -max = no hit.
  }
}

float shadow(void *fragptr, vec3 *org, vec3 *dir) {
  Fragment *frag = reinterpret_cast<Fragment *>(fragptr);

  if (frag->raydepth > 31) {
    // too many reflection
    fprintf(stderr, "[LSGL] Too many ray reflections: %d\n", frag->raydepth);
    assert(0);
    return -1.0f;
  }

  Intersection isect;

  real3 rayorg;
  real3 raydir;

  rayorg[0] = org->v[0];
  rayorg[1] = org->v[1];
  rayorg[2] = org->v[2];

  raydir[0] = dir->v[0];
  raydir[1] = dir->v[1];
  raydir[2] = dir->v[2];

  Ray ray(rayorg, raydir);
  ray.depth = frag->raydepth + 1;
  ray.px = frag->px;
  ray.py = frag->py;
  ray.user_attrib = 0.0f;
  ray.double_sided = frag->doubleSided;

  bool hit = RaytraceEngine::GetRaytraceEngine()->Trace(isect, ray);

  if (hit) {
    // no shading. just return distance
    return isect.t;
  } else {
    return -(std::numeric_limits<float>::max)(); // -max = no hit.
  }
}
