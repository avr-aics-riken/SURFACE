/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_GLES_HANDLE_ALLOCATOR_H__
#define __LSGL_GLES_HANDLE_ALLOCATOR_H__

#include <vector>
#include <cassert>

#include "gles_common.h"

#include "GLES2/gl2.h"

namespace lsgl {

/// Simple handle resource management.
class HandleAllocator {
public:
  // id = 0 is reserved.
  HandleAllocator() : counter_(1){};
  ~HandleAllocator(){};

  /// Allocates handle object.
  GLuint Allocate() {

    GLuint handle = 0;

    if (!freeList_.empty()) {
      // Reuse previously issued handle.
      handle = freeList_.back();
      freeList_.pop_back();
      return handle;
    }

    handle = counter_;
    assert(handle >= 1);
    assert(handle < 0xFFFFFFFF);

    counter_++;

    return handle;
  }

  /// Release handle object.
  void Release(GLuint handle) {
    if (handle == counter_ - 1) {
      if (counter_ > 1) {
        counter_--;
      }
    } else {
      assert(handle >= 1);
      freeList_.push_back(handle);
    }
  }

private:
  std::vector<GLuint> freeList_;
  GLuint counter_;
};
}

#endif // __LSGL_GLES_HANDLE_ALLOCATOR_H__
