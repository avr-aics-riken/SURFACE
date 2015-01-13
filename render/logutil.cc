/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "logutil.h"

#include <exception>

namespace lsgl {
namespace render {

#if DEBUG_LOG_ENABLE

// Static initializer
int Log::indentation_ = 0;
std::ostream *Log::stream_ = &std::cout;

Log::Log(const std::string &context) : context_(context) {
  t_.start();

  writeIndentation();
  *(stream_) << "--> " << context_ << std::endl;
  ++indentation_;
  stream_->flush();
}

Log::~Log() {
  t_.end();
  double msec = t_.msec();
  --indentation_;
  writeIndentation(std::uncaught_exception() ? '*' : ' ');
  *(stream_) << "<-- " << context_;
  *(stream_) << " in " << msec << " [msec(s)]";
  *(stream_) << std::endl;
  stream_->flush();
}

void Log::writeIndentation() { writeIndentation(' '); }

void Log::writeIndentation(const char prefix) {
  *stream_ << prefix;
  for (int i = 0; i < indentation_ * 2; ++i) {
    *stream_ << " ";
  }
}

void Log::message(const std::string &message) {
  writeIndentation();
  *stream_ << message << std::endl;
  stream_->flush();
}

#endif // DEBUG_LOG_ENABLE

} // namespace render
} // namespace lsgl
