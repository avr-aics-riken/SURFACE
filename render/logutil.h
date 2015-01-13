/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LOGUTIL_H__
#define __LOGUTIL_H__

#include "timerutil.h"

#include <string>

//
// NOTE: log class is not thread-safe.
//

#define DEBUG_LOG_ENABLE (1)

#if DEBUG_LOG_ENABLE
#include <iostream>
#endif

namespace lsgl {
namespace render {

class Log {
public:
  Log(const std::string &context);
  ~Log();

  void message(const std::string &message);

private: // Members
  void writeIndentation();
  void writeIndentation(const char prefix);

  static int indentation_;
  const std::string context_;
#ifdef DEBUG_LOG_ENABLE
  static std::ostream *stream_;
#endif

  timerutil t_;
};

} // namespace render
} // namespace lsgl

#endif // __LOGUTIL_H__
