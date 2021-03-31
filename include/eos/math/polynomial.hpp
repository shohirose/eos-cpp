#pragma once

#include <array>  // std::array
#include <vector> // std::vector

namespace eos
{
  /// @brief Computes real roots of a cubic equation
  /// @param[in] a Array of coefficients of a cubic equation
  /// @returns An array of real roots in the ascending order
  ///
  /// The cubic equation takes the form of:
  /// \f[
  ///   x^3 + a[0] x^2 + a[1] x + a[2] = 0.
  /// \f]
  std::vector<double> real_roots(const std::array<double, 3> &a) noexcept;

} // namespace eos