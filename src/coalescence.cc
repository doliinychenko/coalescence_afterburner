#include "coalescence/coalescence.h"

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
}

Coalescence::~Coalescence() {
  std::fclose(output_);
}

ParticleType Coalescence::pdg_to_type(int32_t pdg) {
  switch (pdg) {
    case 2212:  return ParticleType::p;
    case 2112:  return ParticleType::n;
    case 3122:  return ParticleType::la;
    case -2212: return ParticleType::ap;
    case -2112: return ParticleType::an;
    case -3122: return ParticleType::ala;
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
        smash::FourVector r(t, x, y, z), p(p0, px, py, pz);
        smash::FourVector origin(time_last_coll,
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
      const smash::FourVector &p = nucleus.momentum;
      fprintf(output_, "%12.8f %12.8f %12.8f %12.8f %d\n",
          p.x0(), p.x1(), p.x2(), p.x3(), static_cast<int>(nucleus.type));
    }

    // std::cout << "Read in " << n_part_lines << " particles" << std::endl;
  }
  std::cout << event_number_ << " events" <<  std::endl;
  std::fclose(input);
}

void Coalescence::coalesce(const std::vector<Particle> &hadrons,
                           std::vector<Nucleus> &nuclei) {
  nuclei.clear();
  std::vector<Particle> nucleons, antinucleons;
  for (const Particle &hadron : hadrons) {
    // Avoid spectator nucleons. Even if fragmentation of spectators occurs
    // the corresponding nucleons should collide with something.
    if (hadron.pdg_mother1 == 0 && hadron.pdg_mother2 == 0) {
      continue;
    }

    if (hadron.type == ParticleType::p ||
        hadron.type == ParticleType::n) {
      nucleons.push_back(hadron);
    } else if (hadron.type == ParticleType::p ||
               hadron.type == ParticleType::n) {
      antinucleons.push_back(hadron);
    }
  }
  const size_t N_nucleons = nucleons.size();
  std::cout << "Trying to combine " << nucleons.size()
            << " nucleons into deuterons." << std::endl;
  for (size_t i = 0; i < N_nucleons; i++) {
    for (size_t j = 0; j < i; j++) {
      smash::FourVector x1(nucleons[i].origin), x2(nucleons[j].origin),
                        p1(nucleons[i].momentum), p2(nucleons[j].momentum);
      // 1. Boost to the center of mass frame
      const smash::ThreeVector vcm = (p1 + p2).velocity();
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
      smash::ThreeVector r1 = x1.threevec() + (tmax - x1.x0()) * p1.velocity(),
                         r2 = x2.threevec() + (tmax - x2.x0()) * p2.velocity();

      // 4. Check if spatial distance is too large
      if ((r1 - r2).abs() > deuteron_deltax_) {
        continue;
      }

      // 5. No isospin magic so far - todo.
      // Simply make a deuteron

      /*
      std::cout << "Combining " << nucleons[i].momentum << " "
                                << nucleons[i].origin << " "
                                << nucleons[i].pdg_mother1 << " "
                                << nucleons[i].pdg_mother2 << " and "
                                << nucleons[j].momentum << " "
                                << nucleons[j].origin << " "
                                << nucleons[j].pdg_mother1 << " "
                                << nucleons[j].pdg_mother2 << " " << std::endl;
      */
      nuclei.push_back({nucleons[i].momentum + nucleons[j].momentum,
                        NucleusType::d});
    }
  }
}

}  // namescape coalescence
