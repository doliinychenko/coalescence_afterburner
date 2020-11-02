#ifndef COALESCENCE_H
#define COALESCENCE_H

#include <vector>
#include <random>

#include "coalescence/fourvector.h"

namespace coalescence {

enum class ParticleType : char {
  boring, // hadrons not interesting for coalescence
  p,    // proton
  n,    // neutron
  la,   // lambda
  sig0, // Sigma0
  ap,   // anti-proton
  an,   // anti-neuton
  ala,  // anti-lambda
  asig0, // anti-Sigma0
  d,     // deuteron
  t,     // triton
  He3,   // Helium-3
  H3L,   // Hypertriton
  He4_0, // Helium-4 ground state
  // ...
};

struct Particle {
  FourVector momentum;  // 4-momentum
  FourVector origin;    // 4-position of origin
  ParticleType type;
  int32_t pdg_mother1;
  int32_t pdg_mother2;
  bool valid;
};

class Coalescence {
 public:
  Coalescence(const std::string output_file,
              double deuteron_deltap, double deuteron_deltar);
  ~Coalescence();
  static FourVector combined_r(const Particle &h1, const Particle &h2);
  bool check_vicinity(const Particle &h1, const Particle &h2, double deltap, double deltar);
  void coalesce(const std::vector<Particle> &in,
                std::vector<Particle> &out);
  void make_nuclei(const std::string &input_file);
 private:
  static constexpr double hbarc = 0.197327053;

  // random number generation
  std::random_device random_device_;
  std::mt19937 rng_generator_;

  // Particles from how many events will be used for coalescence
  const int n_events_combined_ = 1;

  ParticleType pdg_to_type(int32_t pdg);
  size_t event_number_ = 0;
  FILE *output_;
  // Coalescence parameters
  const double deuteron_deltap_ = 0.44;  // GeV
  const double deuteron_deltar_ = 2.0 * M_PI * hbarc / deuteron_deltap_;  // fm
};

}  // namespace coalescence
#endif  // COALESCENCE_H
