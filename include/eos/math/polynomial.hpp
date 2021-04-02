#pragma once

#include <array>  // std::array
#include <vector> // std::vector

namespace eos
{
  /// @brief Computes real roots of a cubic equation
  /// @param[in] a First coefficient
  /// @param[in] b Second coefficient
  /// @param[in] c Third coefficient
  /// @returns An array of real roots in the ascending order
  ///
  /// The cubic equation takes the form of:
  /// \f[
  ///   x^3 + ax^2 + bx + c = 0.
  /// \f]
  std::vector<double> real_roots(double a, double b, double c) noexcept;

} // namespace eos