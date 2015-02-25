#define _FILE_OFFSET_BITS 64

#include <cstdio>
#include <cassert>

#include "tiny_swc.h"

#include <map>

namespace tinyswc {

TinySWC::TinySWC(std::string filename)
    : parsed_(false),
      filename_(filename) {
}

TinySWC::~TinySWC() {
}

bool TinySWC::Parse() {
  assert(!parsed_ && "do not reuse TinySWC object");

  parsed_ = true;

  FILE* fp;
  if ((fp = fopen(filename_.c_str(), "rb")) == NULL) {
    fprintf(stderr, "cannot open file %s\n", filename_.c_str());
    return true;
  }

  char linebuf[8192];

  std::map<int, SamplePoint> buffer;

  while (fgets(linebuf, 8192, fp) != NULL) {
    if (linebuf[0] == '#') { // comment
      continue;
    }

    int n;
    int type;
    float x, y, z;
    float radius;
    int parent;
    int ret = sscanf(linebuf, "%d %d %f %f %f %f %d\n", &n, &type, &x, &y, &z, &radius, &parent);
    if (ret != 7) {
      fprintf(stderr, ".swc parse failed.\n");
      fclose(fp);
      return false;
    }

    SamplePoint p;
    p.SetX(x);
    p.SetY(y);
    p.SetZ(z);
    p.SetRadius(radius);
    p.SetType(type);
    p.SetParent(parent);
    p.SetNumber(n);

    buffer[n] = p;
  }

  // Construct line segment.
  for (size_t i = 0; i < buffer.size(); i++) {

    SamplePoint& s = buffer[i];
    int parent = s.GetParent();
    if (parent < 0 || parent >= buffer.size()) {
      continue;
    }

    SamplePoint& e = buffer[parent];

    neuronSegments_.push_back(s);
    neuronSegments_.push_back(e);

  }

  return true;
}

}  // namespace tinyswc
