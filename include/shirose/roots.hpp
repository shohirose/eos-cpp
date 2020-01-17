// MIT License
//
// Copyright (c) 2019 Sho Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <array>    // std::array
#include <cmath>    // std::sqrt, std::pow, std::fabs
#include <complex>  // std::complex
#include <vector>   // std::vector

namespace shirose {

/// @brief Computes roots of a cubic equation by using Cardano's formula.
/// @param[in] a Array of coefficients of a cubic equation
/// @returns Array of roots of a cubic equation
///
/// The cubic equation takes the form of:
/// \f[
///   x^3 + a[0] x^2 + a[1] x + a[2] = 0.
/// \f]
template <typename T>
std::array<std::complex<T>, 3> roots(const std::array<T, 3>& a) noexcept {
  const auto p = (3 * a[1] - a[0] * a[0]) / 9;
  const auto q = (27 * a[2] + a[0] * (2 * a[0] * a[0] - 9 * a[1])) / 54;
  // Discriminant of the cubic equation
  const auto disc = p * p * p + q * q;

  using std::pow;
  using std::sqrt;

  const auto s = sqrt(std::complex<T>(disc, 0));
  const auto u1 = pow(-q + s, 1.0 / 3.0);
  const auto u2 = pow(-q - s, 1.0 / 3.0);

  constexpr double sqrt3 = 1.7320508075688772935;
  // The primitive cube root of unity
  const auto w1 = std::complex<T>(-0.5, sqrt3 / 2);
  const auto w2 = std::complex<T>(-0.5, -sqrt3 / 2);

  // Roots based on Cardano's formula
  const auto x1 = u1 + u2 - a[0] / 3;
  const auto x2 = w1 * u1 + w2 * u2 - a[0] / 3;
  const auto x3 = w2 * u1 + w1 * u2 - a[0] / 3;

  return {x1, x2, x3};
}

template <typename T>
int num_of_real_roots(const std::array<T, 3>& a) noexcept {
  // Depressed cubic equation:
  // x^3 + 3px + 2q = 0
  const auto p = (3 * a[1] - a[0] * a[0]) / 9;
  const auto q = (27 * a[2] + a[0] * (2 * a[0] * a[0] - 9 * a[1])) / 54;
  const auto det = p * p * p + q * q;

  if (det == 0) {
    if (p == 0)
      return 1;
    else
      return 2;
  } else if (det > 0) {
    return 3;
  } else {
    return 1;
  }
}

template <typename T, std::size_t N>
std::vector<T> real_roots(const std::array<std::complex<T>, N>& x) noexcept {
  std::vector<T> xreal;
  xreal.reserve(N);
  using std::fabs;
  constexpr double eps = 1e-10;
  for (auto&& xi : x)
    if (fabs(xi.imag()) < eps) xreal.push_back(xi.real());
  return xreal;
}

}  // namespace shirose