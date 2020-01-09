#include "coalescence/coalescence.h"
#include "coalescence/threevector.h"
#include "coalescence/fourvector.h"

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cassert>

namespace coalescence {

Coalescence::Coalescence(const std::string output_file) {
  output_ = std::fopen(output_file.c_str(), "w");
  if (output_ == NULL) {
    throw std::runtime_error("Can't open file " + output_file);
  }
  // Initialize random number generator
  rng_generator_.seed(this->random_device_());
}

Coalescence::~Coalescence() {
  std::fclose(output_);
}

ParticleType Coalescence::pdg_to_type(int32_t pdg) {
  switch (pdg) {
    case 2212:  return ParticleType::p;
    case 2112:  return ParticleType::n;
    case 3122:  return ParticleType::la;
    case 3212:  return ParticleType::sig0;
    case -2212: return ParticleType::ap;
    case -2112: return ParticleType::an;
    case -3122: return ParticleType::ala;
    case -3212:  return ParticleType::asig0;
    default: return ParticleType::boring;
  };
}

void Coalescence::make_nuclei(const std::string &input_file) {
  /*
   *  1. Read event
   *  2. Perform coalescence over particles from the event
   *  3. Write results to output
   *  4. Repeat until the input file is over
   */
  FILE* input = std::fopen(input_file.c_str(), "rb");
  if (input == NULL) {
    throw std::runtime_error("Can't open file " + input_file);
  }

  // Read header
  char magic_number[5], smash_version[20];
  uint16_t format_version, format_variant;
  uint32_t len;
  std::fread(&magic_number[0], 4, 1, input);
  std::fread(&format_version, sizeof(std::uint16_t), 1, input);
  std::fread(&format_variant, sizeof(std::uint16_t), 1, input);
  std::fread(&len, sizeof(std::uint32_t), 1, input);
  std::fread(&smash_version[0], sizeof(char), len, input);

  if (strcmp(magic_number, "SMSH") != 0) {
    throw std::runtime_error(input_file + " is likely not a SMASH binary:" +
                             " magic number does not match ");
  }

  if (format_variant != 1) {
    throw std::runtime_error(input_file + " is not a file of" +
                             " extended SMASH binary format.");
  }

  // std::cout << magic_number << " " << format_version << " "
  //           << format_variant << " " << smash_version << std::endl;

  while (true) {
    char block_type;
    if (!std::fread(&block_type, sizeof(char), 1, input)) {
      break;
    }
    if (block_type == 'f') {
      uint32_t ev;
      double impact_parameter;
      char empty;
      std::fread(&ev, sizeof(std::uint32_t), 1, input);
      std::fread(&impact_parameter, sizeof(double), 1, input);
      std::fread(&empty, sizeof(char), 1, input);
      event_number_++;
      continue;
    }
    assert(block_type == 'p');
    uint32_t n_part_lines;
    std::fread(&n_part_lines, sizeof(std::uint32_t), 1, input);
    std::vector<Particle> hadrons;
    std::vector<Nucleus> nuclei;
    hadrons.clear();
    nuclei.clear();
    for (size_t i = 0; i < n_part_lines; i++) {
      double t, x, y, z, m, p0, px, py, pz,
             form_time, xsecfac, time_last_coll;
      int32_t pdg, id, charge, ncoll, proc_id_origin,
              proc_type_origin, pdg_mother1, pdg_mother2;
      std::fread(&t, sizeof(double), 1, input);
      std::fread(&x, sizeof(double), 1, input);
      std::fread(&y, sizeof(double), 1, input);
      std::fread(&z, sizeof(double), 1, input);
      std::fread(&m, sizeof(double), 1, input);
      std::fread(&p0, sizeof(double), 1, input);
      std::fread(&px, sizeof(double), 1, input);
      std::fread(&py, sizeof(double), 1, input);
      std::fread(&pz, sizeof(double), 1, input);
      std::fread(&pdg, sizeof(std::int32_t), 1, input);
      std::fread(&id,  sizeof(std::int32_t), 1, input);
      std::fread(&charge, sizeof(std::int32_t), 1, input);
      std::fread(&ncoll,  sizeof(std::int32_t), 1, input);
      std::fread(&form_time, sizeof(double), 1, input);
      std::fread(&xsecfac, sizeof(double), 1, input);
      std::fread(&proc_id_origin, sizeof(std::int32_t), 1, input);
      std::fread(&proc_type_origin, sizeof(std::int32_t), 1, input);
      std::fread(&time_last_coll, sizeof(double), 1, input);
      std::fread(&pdg_mother1, sizeof(std::int32_t), 1, input);
      std::fread(&pdg_mother2, sizeof(std::int32_t), 1, input);

      ParticleType hadron_type = pdg_to_type(pdg);
      if (hadron_type != ParticleType::boring) {
        FourVector r(t, x, y, z), p(p0, px, py, pz);
        FourVector origin(time_last_coll,
            r.threevec() - (t - time_last_coll) * p.velocity());
        hadrons.push_back({p, origin, hadron_type, pdg_mother1, pdg_mother2});
        // std::cout << pdg << " " << static_cast<int>(hadron_type) << " "
        //          << r << " " << p << " " << origin << std::endl;
      }
    }
    // All the physics of coalescence happens inside
    coalesce(hadrons, nuclei);

    // Print out nuclei
    fprintf(output_, "# event %lu %lu\n", event_number_, nuclei.size());
    for (const Nucleus &nucleus : nuclei) {
      const FourVector &p = nucleus.momentum;
      fprintf(output_, "%12.8f %12.8f %12.8f %12.8f %d\n",
          p.x0(), p.x1(), p.x2(), p.x3(), static_cast<int>(nucleus.type));
    }

    // std::cout << "Read in " << n_part_lines << " particles" << std::endl;
  }
  std::cout << event_number_ << " events" <<  std::endl;
  std::fclose(input);
}

bool Coalescence::check_vicinity(const Particle &h1,
                                 const Particle &h2,
                                 const Particle &h3,
                                 double deltap, double deltar) {
  FourVector x1(h1.origin), x2(h2.origin), x3(h3.origin),
             p1(h1.momentum), p2(h2.momentum), p3(h3.momentum);
  // 1. Boost to the center of mass frame
  const ThreeVector vcm = (p1 + p2 + p3).velocity();
  p1 = p1.lorentz_boost(vcm);
  p2 = p2.lorentz_boost(vcm);
  p3 = p3.lorentz_boost(vcm);
  x1 = x1.lorentz_boost(vcm);
  x2 = x2.lorentz_boost(vcm);
  x3 = x3.lorentz_boost(vcm);
  if ((p1.threevec() + p2.threevec() + p3.threevec()).sqr() > 1e-12) {
    std::cout << "Something is wrong with cm frame: "
              << p1 + p2 + p3 << std::endl;
  }

  // 2. Check if momentum difference is too large
  if (p1.abs3() > deltap || p2.abs3() > deltap || p3.abs3() > deltap) {
    return false;
  }

  // 3. Roll to the time, when the last nucleon was born
  const double tmax = std::max({x1.x0(), x2.x0(), x3.x0()});
  ThreeVector r1 = x1.threevec() + (tmax - x1.x0()) * p1.velocity(),
              r2 = x2.threevec() + (tmax - x2.x0()) * p2.velocity(),
              r3 = x3.threevec() + (tmax - x3.x0()) * p3.velocity(),
              rcm = (r1 + r2 + r3) / 3.0;
  r1 -= rcm;
  r2 -= rcm;
  r3 -= rcm;

  // 4. Check if spatial distance is too large
  if (r1.abs() > deltar || r2.abs() > deltar || r3.abs() > deltar) {
    return false;
  }

  return true;
}

void Coalescence::coalesce(const std::vector<Particle> &hadrons,
                           std::vector<Nucleus> &nuclei) {
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  nuclei.clear();
  std::vector<Particle> protons, neutrons, antiprotons, antineutrons;
  for (const Particle &hadron : hadrons) {
    // Avoid spectator nucleons. Even if fragmentation of spectators occurs
    // the corresponding nucleons should collide with something.
    if (hadron.pdg_mother1 == 0 && hadron.pdg_mother2 == 0) {
      continue;
    }

    switch (hadron.type) {
      case ParticleType::p: protons.push_back(hadron); break;
      case ParticleType::n: neutrons.push_back(hadron); break;
      case ParticleType::ap: antiprotons.push_back(hadron); break;
      case ParticleType::an: antineutrons.push_back(hadron); break;
      default: ;
    }
  }
  // std::cout << "Trying to combine " << protons.size() << " protons and "
  //          << neutrons.size() << " neutrons into deuterons." << std::endl;
  for (Particle &proton : protons) {
    for (Particle &neutron : neutrons) {
      FourVector x1(proton.origin), x2(neutron.origin),
                        p1(proton.momentum), p2(neutron.momentum);
      // 1. Boost to the center of mass frame
      const ThreeVector vcm = (p1 + p2).velocity();
      p1 = p1.lorentz_boost(vcm);
      p2 = p2.lorentz_boost(vcm);
      x1 = x1.lorentz_boost(vcm);
      x2 = x2.lorentz_boost(vcm);
      if ((p1.threevec() + p2.threevec()).sqr() > 1e-12) {
        std::cout << "Something is wrong with cm frame: "
                  << p1 + p2 << std::endl;
      }

      // 2. Check if momentum difference is too large
      if ((p1.threevec() - p2.threevec()).abs() > deuteron_deltap_) {
        continue;
      }

      // 3. Roll to the time, when the last nucleon was born
      const double tmax = std::max({x1.x0(), x2.x0()});
      ThreeVector r1 = x1.threevec() + (tmax - x1.x0()) * p1.velocity(),
                         r2 = x2.threevec() + (tmax - x2.x0()) * p2.velocity();

      // 4. Check if spatial distance is too large
      if ((r1 - r2).abs() > deuteron_deltax_) {
        continue;
      }

      // 5. Spin average over initial states (* 1/4),
      //    spin sum over final state (* 3), and
      //    isospin projection (* 1/2), see DOI: 10.1103/PhysRevC.53.367
      //   Therfore accept with probability 3/8.
      if (uniform01(rng_generator_) > 3./8.) {
        continue;
      }

      /*
      std::cout << "Combining " << proton.momentum << " "
                                << proton.origin << " "
                                << proton.pdg_mother1 << " "
                                << proton.pdg_mother2 << " and "
                                << neutron.momentum << " "
                                << neutron.origin << " "
                                << neutron.pdg_mother1 << " "
                                << neutron.pdg_mother2 << " " << std::endl;
      */
      nuclei.push_back({proton.momentum + neutron.momentum, NucleusType::d});
    }
  }

  #pragma omp parallel for
  for (size_t i = 0; i < protons.size(); i++) {
    for (size_t j = 0; j < neutrons.size(); j++) {
      for (size_t k = 0; k < j; k++) {
        if (!check_vicinity(protons[i], neutrons[j], neutrons[k],
                            deuteron_deltap_ / std::sqrt(3.0),
                            deuteron_deltax_ / std::sqrt(3.0))) {
          continue;
        }
        // Triton: Spin average over initial states (* 1/8),
        //         spin sum over final state (* 2),
        if (uniform01(rng_generator_) > 1./4.) {
          continue;
        }
        nuclei.push_back({protons[i].momentum +
                          neutrons[j].momentum +
                          neutrons[k].momentum, NucleusType::t});
      }
    }
  }

  #pragma omp parallel for
  for (size_t i = 0; i < neutrons.size(); i++) {
    for (size_t j = 0; j < protons.size(); j++) {
      for (size_t k = 0; k < j; k++) {
        if (!check_vicinity(neutrons[i], protons[j], protons[k], 
                            deuteron_deltap_ / std::sqrt(3.0),
                            deuteron_deltax_ / std::sqrt(3.0))) {
          continue;
        }
        // He3: Spin average over initial states (* 1/8),
        //           spin sum over final state (* 2),
        if (uniform01(rng_generator_) > 1./4.) {
          continue;
        }
        nuclei.push_back({neutrons[i].momentum + 
                          protons[j].momentum + 
                          protons[k].momentum, NucleusType::He3});
      }
    }
  }
}

}  // namescape coalescence
