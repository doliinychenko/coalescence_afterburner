/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_THREEVECTOR_H_
#define SRC_INCLUDE_THREEVECTOR_H_

#include <array>
#include <cmath>
#include <ostream>

namespace coalescence {

/**
 * \ingroup data
 *
 * The ThreeVector class represents a physical three-vector
 * \f$ \vec{x} = (x_1,x_2,x_3)\f$
 * with the components \f$ x_1,x_2,x_3 \f$.
 * It is related to the classes FourVector and Angles,
 * both of which can be converted into a ThreeVector using
 * the 'threevec()' method.
 */
class ThreeVector {
 public:
  /// default constructor (nulls all components)
  ThreeVector() : x_({0., 0., 0.}) {}

  /**
   * Constructor for ThreeVector that takes 3 doubles to set up a ThreeVector
   * with desired values for the components
   *
   * \param[in] y1 value of the first component
   * \param[in] y2 value of the second component
   * \param[in] y3 value of the third component
   */
  ThreeVector(double y1, double y2, double y3) : x_({y1, y2, y3}) {}

  /// access the component at offset \p i.
  double &operator[](std::size_t i) { return x_[i]; }
  /// const overload of the above.
  double operator[](std::size_t i) const { return x_[i]; }

  /// \return first component
  double inline x1() const;
  /// set first component
  void inline set_x1(double x);
  /// \return second component
  double inline x2() const;
  /// set second component
  void inline set_x2(double y);
  /// \return third component
  double inline x3() const;
  /// set third component
  void inline set_x3(double z);
  /// \return the square of the vector (which is a scalar)
  double inline sqr() const;
  /// \return the absolute value
  double inline abs() const;
  /// negation: Returns \f$-\vec x\f$
  ThreeVector inline operator-() const;
  /// increase this vector by \f$\vec v: \vec x^\prime = \vec x + \vec v\f$
  ThreeVector inline operator+=(const ThreeVector &v);
  /// decrease this vector by \f$\vec v: \vec x^\prime = \vec x - \vec v\f$
  ThreeVector inline operator-=(const ThreeVector &v);
  /// scale this vector by \f$a: \vec x^\prime = a \cdot \vec x\f$
  ThreeVector inline operator*=(const double &a);
  /// divide this vector by \f$a: \vec x^\prime = \frac{1}{a} \cdot \vec x\f$
  ThreeVector inline operator/=(const double &a);

  /// \return whether the vector is identical to another vector
  bool operator==(const ThreeVector &rhs) const { return x_ == rhs.x_; }
  /// \return whether the vector is different from another vector
  bool operator!=(const ThreeVector &rhs) const { return x_ != rhs.x_; }

  /// iterates over the components
  using iterator = std::array<double, 3>::iterator;
  /// iterates over the components
  using const_iterator = std::array<double, 3>::const_iterator;

  /**
   * \return an iterator starting at the 0th component.
   *
   * The iterator implements the randomIterator concept. Thus, you can simply
   * write `begin() + 1` to get an iterator that points to the 1st component.
   */
  iterator begin() { return x_.begin(); }

  /// \return an iterator pointing after the 4th component.
  iterator end() { return x_.end(); }

  /// const overload of the above
  const_iterator begin() const { return x_.begin(); }
  /// const overload of the above
  const_iterator end() const { return x_.end(); }

  /// \see begin
  const_iterator cbegin() const { return x_.cbegin(); }
  /// \see end
  const_iterator cend() const { return x_.cend(); }

 private:
  /// the internal storage of the components.
  std::array<double, 3> x_;
};

/**
 * \ingroup logging
 * Writes the three components of the vector to the output stream.
 */
std::ostream &operator<<(std::ostream &, const ThreeVector &);

double inline ThreeVector::x1() const { return x_[0]; }

void inline ThreeVector::set_x1(const double x) { x_[0] = x; }

double inline ThreeVector::x2() const { return x_[1]; }

void inline ThreeVector::set_x2(const double y) { x_[1] = y; }

double inline ThreeVector::x3() const { return x_[2]; }

void inline ThreeVector::set_x3(const double z) { x_[2] = z; }

ThreeVector inline ThreeVector::operator-() const {
  ThreeVector neg(-x_[0], -x_[1], -x_[2]);
  return neg;
}

ThreeVector inline ThreeVector::operator+=(const ThreeVector &v) {
  x_[0] += v.x_[0];
  x_[1] += v.x_[1];
  x_[2] += v.x_[2];
  return *this;
}

/// \return sum of two three-vectors \f$ \vec{a} + \vec{b} \f$.
ThreeVector inline operator+(ThreeVector a, const ThreeVector &b) {
  a += b;
  return a;
}

ThreeVector inline ThreeVector::operator-=(const ThreeVector &v) {
  x_[0] -= v.x_[0];
  x_[1] -= v.x_[1];
  x_[2] -= v.x_[2];
  return *this;
}

/// \return difference between two three-vectors \f$ \vec{a} - \vec{b} \f$.
ThreeVector inline operator-(ThreeVector a, const ThreeVector &b) {
  a -= b;
  return a;
}

ThreeVector inline ThreeVector::operator*=(const double &a) {
  x_[0] *= a;
  x_[1] *= a;
  x_[2] *= a;
  return *this;
}

/// multiply a three-vector by constant factor \f$ b \vec{a} \f$.
inline ThreeVector operator*(ThreeVector a, const double &b) {
  a *= b;
  return a;
}

/// multiply a three-vector by constant factor \f$ a \vec{b} \f$.
inline ThreeVector operator*(const double &a, ThreeVector b) {
  b *= a;
  return b;
}

/**
 * \return inner product of two three-vectors
 * \f$ \vec{a} \cdot \vec{b} \f$
 */
inline double operator*(ThreeVector a, const ThreeVector &b) {
  return a.x1() * b.x1() + a.x2() * b.x2() + a.x3() * b.x3();
}

ThreeVector inline ThreeVector::operator/=(const double &a) {
  const double a_inv = 1.0 / a;
  x_[0] *= a_inv;
  x_[1] *= a_inv;
  x_[2] *= a_inv;
  return *this;
}

/// divide a three-vector by constant factor \f$ \vec{a} / b \f$.
ThreeVector inline operator/(ThreeVector a, const double &b) {
  a /= b;
  return a;
}

double inline ThreeVector::sqr() const { return (*this) * (*this); }

double inline ThreeVector::abs() const { return std::sqrt((*this) * (*this)); }

}  // namespace coalescence

#endif  // SRC_INCLUDE_THREEVECTOR_H_
