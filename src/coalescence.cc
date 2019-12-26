#include "coalescence/coalescence.h"

#include <stdio.h>

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
}

void Coalescence::coalesce(const std::vector<Particle> &in,
                           std::vector<Nucleus> &out) {
}

}  // namescape coalescence
