#ifndef TINY_SWC_H__
#define TINY_SWC_H__

#include <string>
#include <vector>

// SWC spec:
// http://research.mssm.edu/cnic/swc.html

namespace tinyswc {

class SamplePoint {
 public:
  int GetNumber() const {
    return number_;
  }
  int GetType() const {
    return type_;
  }
  double GetRadius() const {
    return radius_;
  }
  double GetX() const {
    return x_;
  }
  double GetY() const {
    return y_;
  }
  double GetZ() const {
    return z_;
  }
  int GetParent() const {
    return parent_;
  }
  void SetNumber(int n) {
    number_ = n;
  }
  void SetType(int t) {
    type_ = t;
  }
  void SetRadius(double r) {
    radius_ = r;
  }
  void SetX(double x) {
    x_ = x;
  }
  void SetY(double y) {
    y_ = y;
  }
  void SetZ(double z) {
    z_ = z;
  }
  void SetParent(int p) {
    parent_ = p;
  }
 private:
  int number_;
  int type_;
  double radius_;
  double x_;
  double y_;
  double z_;
  int parent_;
};

class TinySWC {
 public:
  TinySWC(std::string filename);
  ~TinySWC();

  // Parses the PDB file. Returns true if failed.
  bool Parse();

  // Returns the list of segments of neuron sample points
  std::vector<SamplePoint>& GetNeuronSegments() {
    return neuronSegments_;
  }

 private:
  bool parsed_;
  std::string filename_;
  std::vector<SamplePoint> neuronSegments_;
};

}  // namespace tinyswc

#endif // TINY_SWC_H__
