#pragma once

#include <array>   // std::array
#include <cmath>   // std::sqrt, std::pow, std::fabs
#include <complex> // std::complex
#include <vector>  // std::vector
#include <boost/math/constants/constants.hpp>

namespace eos
{

  /// @brief Computes roots of a cubic equation by using Cardano's formula.
  /// @param[in] a Array of coefficients of a cubic equation
  /// @returns Array of complex roots of a cubic equation
  ///
  /// The cubic equation takes the form of:
  /// \f[
  ///   x^3 + a[0] x^2 + a[1] x + a[2] = 0.
  /// \f]
  template <typename T>
  auto roots(const std::array<T, 3> &a) noexcept -> std::array<std::complex<T>, 3>
  {
    const auto p = (3 * a[1] - a[0] * a[0]) / 9;
    const auto q = (27 * a[2] + a[0] * (2 * a[0] * a[0] - 9 * a[1])) / 54;
    // Discriminant of the cubic equation
    const auto disc = p * p * p + q * q;

    using std::complex;
    using std::pow;
    using std::sqrt;

    const auto s = sqrt(complex<T>(disc, 0));
    const auto u1 = pow(-q + s, 1.0 / 3.0);
    const auto u2 = pow(-q - s, 1.0 / 3.0);

    constexpr auto sqrt3 = boost::math::constants::root_three<T>();
    // The primitive cube root of unity
    const auto w1 = complex<T>(-0.5, sqrt3 / 2);
    const auto w2 = complex<T>(-0.5, -sqrt3 / 2);

    // Roots based on Cardano's formula
    const auto x1 = u1 + u2 - a[0] / 3;
    const auto x2 = w1 * u1 + w2 * u2 - a[0] / 3;
    const auto x3 = w2 * u1 + w1 * u2 - a[0] / 3;

    return {x1, x2, x3};
  }

  /// @brief Computes real roots of a cubic equation
  /// @param[in] a Coefficients of a cubic equation
  /// @return An array of real roots of a cubic equation
  template <typename T>
  std::vector<T> real_roots(const std::array<T, 3> &a) noexcept
  {
    const auto x = roots(a);
    std::vector<T> xreal;
    xreal.reserve(3);
    using std::fabs;
    constexpr auto eps = static_cast<T>(1e-10);
    for (auto &&xi : x)
    {
      if (fabs(xi.imag()) < eps)
      {
        xreal.push_back(xi.real());
      }
    }
    return xreal;
  }

} // namespace eos