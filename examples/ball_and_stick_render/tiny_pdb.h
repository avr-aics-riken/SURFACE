#ifndef TINY_PDB_H__
#define TINY_PDB_H__

#include <string>
#include <vector>

#include "tiny_pdb_math.hxx"

namespace tinypdb {

class Atom {
 public:
  Atom() : visited_(false) { }
  ~Atom() { }

  int GetSerialNumber() const {
    return serial_number_;
  }
  const std::string& GetName() const {
    return name_;
  }
  const std::string& GetIsomer() const {
    return isomer_;
  }
  const std::string& GetResidueName() const {
    return residue_name_;
  }
  const std::string& GetChainName() const {
    return chain_name_;
  }
  int GetResidueNumber() const {
    return residue_number_;
  }
  const std::string& GetResidueInsertion() const {
    return residue_insertion_;
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
  double GetOccupancy() const {
    return occupancy_;
  }
  double GetTFactor() const {
    return t_factor_;
  }
  const std::string& GetExtra() const {
    return extra_;
  }
  const std::string& GetElementSymbol() const {
    return element_symbol_;
  }

  // For hash grid library from SmallVCM
  Vec3f GetPosition() const {
    return Vec3f(x_, y_, z_);
  }

  double GetDistance(const Atom& target) const {
    return sqrt(
        pow(GetX() - target.GetX(), 2.0) +
        pow(GetY() - target.GetY(), 2.0) +
        pow(GetZ() - target.GetZ(), 2.0));
  }

  const std::vector<Atom*>& GetBonds() {
    return bonds_;
  }

  void SetSerialNumber(int serial_number) {
    serial_number_ = serial_number;
  }
  std::string& MutableName() {
    return name_;
  }
  std::string& MutableIsomer() {
    return isomer_;
  }
  std::string& MutableResidueName() {
    return residue_name_;
  }
  std::string& MutableChainName() {
    return chain_name_;
  }
  void SetResidueNumber(int residue_number) {
    residue_number_ = residue_number;
  }
  std::string& MutableResidueInsertion() {
    return residue_insertion_;
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
  void SetOccupancy(double occupancy) {
    occupancy_ = occupancy;
  }
  void SetTFactor(double t_factor) {
    t_factor_ = t_factor;
  }
  std::string& MutableExtra() {
    return extra_;
  }
  std::string& MutableElementSymbol() {
    return element_symbol_;
  }

  void AddBond(Atom* bonded) {
    bonds_.push_back(bonded);
  }

  bool Visited() const {
    return visited_;
  }

  void SetVisited(bool flag) {
    visited_ = flag;
  }

 private:
  int serial_number_;
  std::string name_;
  std::string isomer_;
  std::string residue_name_;
  std::string chain_name_;
  int residue_number_;
  std::string residue_insertion_;
  double x_;
  double y_;
  double z_;
  double occupancy_;
  double t_factor_;
  std::string extra_;
  std::string element_symbol_;

  bool visited_; // flag used for bond traversal.

  std::vector<Atom*> bonds_;
};

class TinyPDB {
 public:
  TinyPDB(std::string filename);
  ~TinyPDB();

  // Parses the PDB file. Returns true if failed.
  bool Parse(bool isBondGenerationEnabled);

  std::vector<Atom>& GetAtoms() {
    return atoms_;
  }

 private:
  bool ParseLine();
  Atom ParseAtom();

  void Identifier(std::string& identifier);
  void FixedIdentifier(std::string& identifier, int length, bool skipWhiteSpace = true);
  int Integer();
  double Double();
  void Trim();
  void SkipLine();

  bool Match(const char* keyword);
  char Char(int index = 0);
  void Consume(int amount = 1);

  void GenerateBonds();

  bool parsed_;
  std::string filename_;
  FILE* fp_;
  char* content_;
  long long int length_;
  long long int pos_;

  std::vector<Atom> atoms_;
};

}  // namespace tinypdb

#endif
