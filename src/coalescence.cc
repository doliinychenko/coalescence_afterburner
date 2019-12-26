#include "coalescence/coalescence.h"

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
      std::fread(&pdg, sizeof(std::uint32_t), 1, input);
      std::fread(&id,  sizeof(std::uint32_t), 1, input);
      std::fread(&charge, sizeof(std::uint32_t), 1, input);
      std::fread(&ncoll,  sizeof(std::uint32_t), 1, input);
      std::fread(&form_time, sizeof(double), 1, input);
      std::fread(&xsecfac, sizeof(double), 1, input);
      std::fread(&proc_id_origin, sizeof(std::uint32_t), 1, input);
      std::fread(&proc_type_origin, sizeof(std::uint32_t), 1, input);
      std::fread(&time_last_coll, sizeof(double), 1, input);
      std::fread(&pdg_mother1, sizeof(std::uint32_t), 1, input);
      std::fread(&pdg_mother2, sizeof(std::uint32_t), 1, input);
    }
    std::cout << "Read in " << n_part_lines << " particles" << std::endl;
  }
  std::cout << event_number_ << " events" <<  std::endl;
  std::fclose(input);
}

void Coalescence::coalesce(const std::vector<Particle> &in,
                           std::vector<Nucleus> &out) {
}

}  // namescape coalescence
