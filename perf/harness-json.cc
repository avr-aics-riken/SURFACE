#include "harness.h"

#include <string>
#include <cstring>
#include <cstdio>

char* log_json(const char* title, const char* variant, int64_t cpu_ns)
{
  char buf[1024];
  sprintf(buf, "%lld", cpu_ns);
  std::string s = "{\"title\": \"" + std::string(title) + "\", \"variant\": \"" + std::string(variant) + "\", \"cpu_ns\": " + std::string(buf) + "}";
  return strdup(s.c_str());
}
