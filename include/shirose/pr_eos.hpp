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

#ifndef SHIROSE_PR_EOS_HPP
#define SHIROSE_PR_EOS_HPP

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "shirose/alpha_functions.hpp"  // shirose::alpha::peng_robinson_1976
#include "shirose/cubic_eos.hpp"        // shirose::cubic_eos

namespace shirose {

/// @brief Peng-Robinson EoS policy for cubic_eos.
/// @tparam T Value type
template <typename T>
class peng_robinson {
 public:
  /// Constant for attraction parameter
  static constexpr double omega_a = 0.45724;
  /// Constant for repulsion parameter
  static constexpr double omega_b = 0.07780;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / ((v - b) * (v + b) + 2 * b * v);
  }

  /// @brief Computes coeficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static std::array<T, 3> cubic_eq(const T &a, const T &b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static T fugacity_coeff(const T &z, const T &a, const T &b) noexcept {
    using std::exp;
    using std::log;
    /// Square root of 2 (two).
    constexpr double sqrt2 = 1.4142135623730950488;
    constexpr double delta1 = 1 + sqrt2;
    constexpr double delta2 = 1 - sqrt2;
    return exp(z - 1 - log(z - b) -
               a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z - delta2 * b)));
  }
};

/// @brief Peng-Robinson equation of state.
/// @tparam T Value type
template <typename T>
using pr_eos = cubic_eos<T, peng_robinson<T>, alpha::peng_robinson_1976<T>>;

/// @brief Makes Peng-Robinson EoS
/// @tparam T Value type
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename T>
pr_eos<T> make_pr_eos(const T &pc, const T &tc, const T &omega) {
  return {pc, tc, alpha::peng_robinson_1976<T>{omega}};
}

}  // namespace shirose

#endif  // SHIROSE_PR_EOS_HPP