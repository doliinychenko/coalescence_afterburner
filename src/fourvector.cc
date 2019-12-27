/*
 *    Copyright (c) 2012-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <iostream>

#include "coalescence/fourvector.h"

namespace coalescence {

FourVector FourVector::lorentz_boost(const ThreeVector& v) const {
  const double velocity_squared = v.sqr();

  const double gamma =
      velocity_squared < 1. ? 1. / std::sqrt(1. - velocity_squared) : 0;

  // this is used four times in the Vector:
  const double xprime_0 = gamma * (this->x0() - this->threevec() * v);
  // this is the part of the space-like components that is always the same:
  const double constantpart = gamma / (gamma + 1) * (xprime_0 + this->x0());
  return FourVector(xprime_0, this->threevec() - v * constantpart);
}

std::ostream& operator<<(std::ostream& out, const FourVector& vec) {
  out << '(';
  for (auto x : vec) {
    out << x << " ";
  }
  return out << ')';
}

}  // namespace coalescence
