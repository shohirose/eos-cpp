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

#ifndef SHIROSE_VDW_EOS_HPP
#define SHIROSE_VDW_EOS_HPP

#include <array>  // std::array
#include <cmath>  // std::exp, std::log

#include "shirose/cubic_eos.hpp"  // shirose::cubic_eos

namespace shirose {

/// @brief Van der Waals EoS policy for cubic_eos.
/// @tparam T Value type
/// @tparam F Function to compute alpha.
template <typename T>
struct van_der_waals {
  /// Constant for attraction parameter
  static constexpr double omega_a = 0.421875;
  /// Constant for respulsion parameter
  static constexpr double omega_b = 0.125;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / (v * v);
  }

  /// @brief Computes coefficients of cubic equation
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<T, 3> cubic_eq(const T &a, const T &b) noexcept {
    return {-b - 1, a, -a * b};
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static T fugacity_coeff(const T &z, const T &a, const T &b) noexcept {
    using std::exp;
    using std::log;
    return exp(-log(z - b) - a / z + z - 1);
  }

  // Member Functions

  /// @brief Computes temperature correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor
  T alpha(const T &) const noexcept { return 1; }
};

/// @brief Van der Waals equation of state.
/// @tparam T Value type
template <typename T>
using vdw_eos = cubic_eos<T, van_der_waals<T>>;

/// @brief Makes van der Waals EoS
/// @tparam T Value type
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename T>
vdw_eos<T> make_vdw_eos(const T &pc, const T &tc) {
  return {pc, tc, van_der_waals<T>{}};
}

}  // namespace shirose

#endif  // SHIROSE_VDW_EOS_HPP