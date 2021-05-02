#pragma once

#include <array>
#include <complex>
#include <vector>

namespace eos {

/// @brief Cubic equation
///
/// \f[ x^3 + a x^2 + b x + c = 0 \f]
class cubic_equation {
 public:
  /// @{
  /// @name Coefficients

  double a;
  double b;
  double c;
  /// @}

  cubic_equation() = default;
  cubic_equation(const cubic_equation&) = default;
  cubic_equation(cubic_equation&&) = default;

  cubic_equation(double a_, double b_, double c_) : a{a_}, b{b_}, c{c_} {}

  cubic_equation& operator=(const cubic_equation&) = default;
  cubic_equation& operator=(cubic_equation&&) = default;

  /// @brief Computes real roots in the ascending order
  std::vector<double> real_roots() const;

  /// @brief Computes complex roots
  std::array<std::complex<double>, 3> complex_roots() const;
};

}  // namespace eos
