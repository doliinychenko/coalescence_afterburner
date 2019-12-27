#ifndef COALESCENCE_H
#define COALESCENCE_H

#include <vector>

#include "smash/fourvector.h"

namespace coalescence {

enum class ParticleType : char {
  boring, // hadrons not interesting for coalescence
  p,    // proton
  n,    // neutron
  la,   // lambda
  ap,   // anti-proton
  an,   // anti-neuton
  ala,  // anti-lambda
};

enum class NucleusType : char {
  d,     // deuteron
  t,     // triton
  He3,   // Helium-3
  H3L,   // Hypertriton
  He4_0, // Helium-4 ground state
  // ...
};

struct Particle {
  smash::FourVector momentum;  // 4-momentum
  smash::FourVector origin;    // 4-position of origin
  ParticleType type;
  int32_t pdg_mother1;
  int32_t pdg_mother2;
};

struct Nucleus {
  smash::FourVector momentum;  // 4-momentum
  NucleusType type;
};

class Coalescence {
 public:
  Coalescence(const std::string output_file);
  ~Coalescence();
  void coalesce(const std::vector<Particle> &in,
                std::vector<Nucleus> &out);
  void make_nuclei(const std::string &input_file);
 private:
  ParticleType pdg_to_type(int32_t pdg);
  size_t event_number_ = 0;
  FILE *output_;
  // Coalescence parameters
  const double deuteron_deltap_ = 0.3;  // GeV
  const double deuteron_deltax_ = 3.5;  // fm
};

}  // namespace coalescence
#endif  // COALESCENCE_H
