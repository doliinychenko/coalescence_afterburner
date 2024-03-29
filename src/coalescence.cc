#include "coalescence/coalescence.h"
#include "coalescence/threevector.h"
#include "coalescence/fourvector.h"

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cassert>

namespace coalescence {

Coalescence::Coalescence(const std::string output_file,
  double deuteron_deltap, double deuteron_deltar,
  bool probabilistic) :
    deuteron_deltap_(deuteron_deltap),
    deuteron_deltar_(deuteron_deltar),
    probabilistic_(probabilistic) {
  output_ = std::fopen(output_file.c_str(), "w");
  if (output_ == NULL) {
    throw std::runtime_error("Can't open file " + output_file);
  }
  // Initialize random number generator
  rng_generator_.seed(this->random_device_());
  for (int i = 0; i < y_nbins_; i++) {
    proton_y_[i] = 0.0;
    deuteron_y_[i] = 0.0;
    triton_y_[i] = 0.0;
  }
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
  magic_number[4] = '\x00';
  std::fread(&format_version, sizeof(std::uint16_t), 1, input);
  // std::cout << "SMASH output version " << format_version << std::endl;
  std::fread(&format_variant, sizeof(std::uint16_t), 1, input);
  std::fread(&len, sizeof(std::uint32_t), 1, input);
  std::fread(&smash_version[0], sizeof(char), len, input);

  if (strcmp(magic_number, "SMSH") != 0) {
    std::cout << "Magic number = " << magic_number << std::endl;
    throw std::runtime_error(input_file + " is likely not a SMASH binary:" +
                             " magic number does not match ");
  }

  if (format_variant != 1) {
    throw std::runtime_error(input_file + " is not a file of" +
                             " extended SMASH binary format.");
  }

  // std::cout << magic_number << " " << format_version << " "
  //           << format_variant << " " << smash_version << std::endl;

  std::vector<Particle> hadrons, nuclei;
  hadrons.clear();
  nuclei.clear();
 
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
      if (format_version > 6) {
        std::fread(&empty, sizeof(char), 1, input);
      }
      event_number_++;
      continue;
    }
    if (block_type != 'p') {
      break;
    }
    uint32_t n_part_lines;
    std::fread(&n_part_lines, sizeof(std::uint32_t), 1, input);
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
        hadrons.push_back({p, origin, hadron_type,
                           pdg_mother1, pdg_mother2, 1.0, true});
        // std::cout << pdg << " " << static_cast<int>(hadron_type) << " "
        //          << r << " " << p << " " << origin << std::endl;
      }
    }
    // std::cout << "Read in " << n_part_lines << " particles" << std::endl;

    // All the physics of coalescence happens inside
    if (event_number_ % n_events_combined_ == 0) {
      if (!probabilistic_) {
        coalesce(hadrons, nuclei);
      } else {
        coalesce_probabilistic(hadrons, nuclei);
      }
      // Print out nuclei
      fprintf(output_, "# event %lu %lu\n", event_number_, nuclei.size());
      for (const Particle &nucleus : nuclei) {
        const FourVector &p = nucleus.momentum;
        add_to_histograms(nucleus);
        fprintf(output_, "%12.8f %12.8f %12.8f %12.8f %d %12.8f\n",
            p.x0(), p.x1(), p.x2(), p.x3(), static_cast<int>(nucleus.type), nucleus.weight);
      }
      for (const Particle &hadron : hadrons) {
        add_to_histograms(hadron);
      }
      hadrons.clear();
      nuclei.clear();
    }
  }
  // std::cout << event_number_ << " events" <<  std::endl;
  std::fclose(input);
}

bool Coalescence::check_vicinity(const Particle &h1,
                                 const Particle &h2,
                                 double deltap,
                                 double deltar) {
  // const double deltar = 2.0 * M_PI * hbarc / deltap;
  FourVector x1(h1.origin), x2(h2.origin),
             p1(h1.momentum), p2(h2.momentum);
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
  if ((p1.threevec() - p2.threevec()).abs() > deltap) {
    return false;
  }

  // 3. Roll to the time, when the last hadron was born
  const double tmax = std::max({x1.x0(), x2.x0()});
  ThreeVector r1 = x1.threevec() + (tmax - x1.x0()) * p1.velocity(),
              r2 = x2.threevec() + (tmax - x2.x0()) * p2.velocity();

  // 4. Check if spatial distance is too large
  if ((r1 - r2).abs() > deltar) {
    return false;
  }

  // 5. Check if any of these particles was already coalesced earlier
  if (!h1.valid || !h2.valid) {
    return false;
  }

  return true;
}

double Coalescence::get_pair_weight(const Particle &h1,
                                  const Particle &h2) {
  FourVector x1(h1.origin), x2(h2.origin),
             p1(h1.momentum), p2(h2.momentum);
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

  // 2. Get momentum difference, 0.25 because q = |p1-p2|/2
  const double dp2 = (p1.threevec() - p2.threevec()).sqr() * 0.25;

  // 3. Roll to the time, when the last hadron was born
  const double tmax = std::max({x1.x0(), x2.x0()});
  ThreeVector r1 = x1.threevec() + (tmax - x1.x0()) * p1.velocity(),
              r2 = x2.threevec() + (tmax - x2.x0()) * p2.velocity();

  // 4. Get spatial distance
  const double dr2 = (r1 - r2).sqr();

  constexpr double d2 = 3.2 * 3.2;  // [fm^2], see 2012.04352
  return 3.0 * std::exp(- dr2 / d2 - dp2 * d2 / (hbarc * hbarc));
}


FourVector Coalescence::combined_r(const Particle &h1, const Particle &h2) {
  FourVector x1(h1.origin), x2(h2.origin),
             p1(h1.momentum), p2(h2.momentum);
  const double tmax = std::max({x1.x0(), x2.x0()});
  ThreeVector r1 = x1.threevec() + (tmax - x1.x0()) * p1.velocity(),
              r2 = x2.threevec() + (tmax - x2.x0()) * p2.velocity();
  return FourVector(tmax, 0.5 * (r1 + r2));
}

void Coalescence::coalesce_probabilistic(const std::vector<Particle> &hadrons,
                           std::vector<Particle> &nuclei) {
  nuclei.clear();
  std::vector<Particle> nucleons;
  nucleons.clear();

  for (const Particle &hadron : hadrons) {
    // Avoid spectator nucleons. Even if fragmentation of spectators occurs
    // the corresponding nucleons should collide with something.
    // Be careful not to reject nucleons born from hydro, that also have
    // pdg_mother == 0.
    if (hadron.pdg_mother1 == 0 && hadron.pdg_mother2 == 0 &&
        hadron.momentum.x1() == 0.0 && hadron.momentum.x2() == 0) {
      continue;
    }

    switch (hadron.type) {
      case ParticleType::p: nucleons.push_back(hadron); break;
      case ParticleType::n: nucleons.push_back(hadron); break;
      default: ;
    }
  }
  // std::cout << "Trying to combine " << protons.size() << " protons and "
  //          << neutrons.size() << " neutrons into deuterons." << std::endl;
  size_t N = nucleons.size();
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < i; j++) {
      const double w = get_pair_weight(nucleons[i], nucleons[j]);
      if (w < 1e-6) {
        continue;
      }
      nuclei.push_back({nucleons[i].momentum + nucleons[j].momentum,
                        combined_r(nucleons[i], nucleons[j]),
                        ParticleType::d, static_cast<int>(nucleons[i].type),
                        static_cast<int>(nucleons[j].type), w, true});
    }
  }


}

void Coalescence::coalesce(const std::vector<Particle> &hadrons,
                           std::vector<Particle> &nuclei) {
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  nuclei.clear();
  std::vector<Particle> protons, neutrons, antiprotons, antineutrons;
  for (const Particle &hadron : hadrons) {
    // Avoid spectator nucleons. Even if fragmentation of spectators occurs
    // the corresponding nucleons should collide with something.
    // Be careful not to reject nucleons born from hydro, that also have
    // pdg_mother == 0.
    if (hadron.pdg_mother1 == 0 && hadron.pdg_mother2 == 0 &&
        hadron.momentum.x1() == 0.0 && hadron.momentum.x2() == 0) {
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
      // Spin average over initial states (* 1/4),
      // spin sum over final state (* 3), and
      // isospin projection (* 1/2), see DOI: 10.1103/PhysRevC.53.367
      // Therfore accept with probability 3/8.
      if (uniform01(rng_generator_) < 3./8. &&
          check_vicinity(proton, neutron, deuteron_deltap_, deuteron_deltar_)) { /*
        std::cout << "Combining " << proton.momentum << " "
                                  << proton.origin << " "
                                  << proton.pdg_mother1 << " "
                                  << proton.pdg_mother2 << " and "
                                  << neutron.momentum << " "
                                  << neutron.origin << " "
                                  << neutron.pdg_mother1 << " "
                                  << neutron.pdg_mother2 << " " << std::endl;
        */
        proton.valid = false;
        neutron.valid = false;
        nuclei.push_back({proton.momentum + neutron.momentum,
                          combined_r(proton, neutron),
                          ParticleType::d, 2212, 2112, 1.0, true});
      }
    }
  }

  std::vector<Particle> deuterons;
  for (const Particle &nucleus : nuclei) {
    if (nucleus.type == ParticleType::d) {
      deuterons.push_back(nucleus);
    }
  }

  for (Particle &deuteron : deuterons) {
    if (!deuteron.valid) {
      continue;
    }
    for (Particle &proton : protons) {
      if (!proton.valid) {
        continue;
      }
      if (uniform01(rng_generator_) < 1./4. &&
        check_vicinity(deuteron, proton, deuteron_deltap_, deuteron_deltar_)) {
        deuteron.valid = false;
        proton.valid = false;
        nuclei.push_back({proton.momentum + deuteron.momentum,
                          combined_r(proton, deuteron),
                          ParticleType::He3, 1000010020, 2212, 1.0, true});
      }
    }
  }

  for (Particle &deuteron : deuterons) {
    if (!deuteron.valid) {
      continue;
    }
    for (Particle &neutron : neutrons) {
      if (!neutron.valid) {
        continue;
      }
      if (uniform01(rng_generator_) < 1./4. &&
        check_vicinity(deuteron, neutron, deuteron_deltap_, deuteron_deltar_)) {
        deuteron.valid = false;
        neutron.valid = false;
        nuclei.push_back({neutron.momentum + deuteron.momentum,
                          combined_r(neutron, deuteron),
                          ParticleType::t, 1000010020, 2112, 1.0, true});
      }
    }
  }

}

void Coalescence::add_to_histograms(const Particle &part) {
  if (!part.valid) {
    return;
  }
  const FourVector p = part.momentum;
  const double y = 0.5 * std::log((p[0] + p[3]) / (p[0] - p[3]));
  const int i = std::floor((y - y_min_) / (y_max_ - y_min_) * y_nbins_);
  if (part.type == ParticleType::p) {
    proton_y_[i] += part.weight;
  } else if (part.type == ParticleType::d) {
    deuteron_y_[i] += part.weight;
  } else if (part.type == ParticleType::t) {
    triton_y_[i] += part.weight;
  }
}

void Coalescence::print_histograms() {
  const double dy = (y_max_ - y_min_) / y_nbins_;
  for (int i = 0; i < y_nbins_; i++) {
    proton_y_[i]   /= (event_number_ * dy);
    deuteron_y_[i] /= (event_number_ * dy);
    triton_y_[i]   /= (event_number_ * dy);
  }
  printf("#y, dN/dy for p,d,t;  p*t/d^2\n");
  for (int i = 0; i < y_nbins_; i++) {
    const double y = y_min_ + (y_max_ - y_min_) / y_nbins_ * (i + 0.5);
    double ptd2 = 0.0;
    if (deuteron_y_[i] > 0.0) {
      ptd2 = proton_y_[i] * triton_y_[i] / deuteron_y_[i] / deuteron_y_[i];
    }
    printf("%8.3f %10.1f %10.1f %10.1f %10.4f\n", y, proton_y_[i], deuteron_y_[i], triton_y_[i], ptd2);
  }
}

}  // namescape coalescence
