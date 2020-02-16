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

#include <array>                               // std::array
#include <boost/math/constants/constants.hpp>  // boost::math::constants::root_two
#include <cmath>                               // std::sqrt, std::exp, std::log

#include "eos/cubic_eos.hpp"  // eos::cubic_eos_base

namespace eos {

template <typename T>
class peng_robinson_eos;

/// @brief Policy class for Peng-Robinson EoS.
/// @tparam T Value type
template <typename T>
struct peng_robinson_eos_policy {
  using derived_type = peng_robinson_eos<T>;
  using value_type = T;

  static constexpr double omega_a = 0.45724;
  static constexpr double omega_b = 0.07780;

  static constexpr double sqrt2 = boost::math::constants::root_two<double>();
  static constexpr double delta1 = 1 + sqrt2;
  static constexpr double delta2 = 1 - sqrt2;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / (v * (v + b) + b * (v - b));
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
    return exp(z - 1 - log(z - b) - q(z, a, b));
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static T residual_enthalpy(const T &z, const T &t, const T &a, const T &b,
                             const T &beta) noexcept {
    using std::log;
    return gas_constant * t * (z - 1 - (1 - beta) * q(z, a, b));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static T residual_entropy(const T &z, const T &a, const T &b,
                            const T &beta) noexcept {
    using std::log;
    return gas_constant * (log(z - b) + beta * q(z, a, b));
  }

  /*
  /// @brief Computes residual molar specific heat at constant volume
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] gamma Temperature correction factor
  static T residual_specific_heat_at_const_volume(const T &z, const T &a,
                                                  const T &b,
                                                  const T &gamma) noexcept {
    using std::log;
    return gas_constant * gamma * a / (2 * sqrt2 * b) *
           log((z + delta1 * b) / (z + delta2 * b));
  }
  */

 private:
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static T q(const T &z, const T &a, const T &b) noexcept {
    return a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z + delta2 * b));
  }
};

/// @brief Peng-Robinson EoS.
/// @tparam T Value type
template <typename T>
class peng_robinson_eos : public cubic_eos_base<peng_robinson_eos_policy<T>> {
 public:
  using base_type = cubic_eos_base<peng_robinson_eos_policy<T>>;

  // Constructors

  peng_robinson_eos() = default;

  /// @brief Constructs Peng-Robinson EoS
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  peng_robinson_eos(const T &pc, const T &tc, const T &omega)
      : base_type{pc, tc}, omega_{omega}, m_{m(omega)} {}

  // Member functions

  /// @brief Set parameters
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  void set_params(const T &pc, const T &tc, const T &omega) noexcept {
    this->base_type::set_params(pc, tc);
    omega_ = omega;
    m_ = m(omega);
  }

 private:
  friend base_type;

  // Member functions

  /// @brief Computes the correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  T alpha(const T &tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// @brief Computes \f$ \beta = \frac{\mathrm{d} \ln \alpha}{\mathrm{d} \ln T}
  /// \f$
  /// @param[in] tr Reduced temperature
  T beta(const T &tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -m_ * sqrt_tr / a;
  }

  /*
  /// @brief Computes  \f[ \gamma = \frac{T_r^2}{\alpha} \cdot
  /// \frac{\mathrm{d}^2 \alpha}{\mathrm{d} T_r^2} \f]
  /// @param[in] tr Reduced temperature
  T gamma(const T &tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    const auto alpha = a * a;
    return m_ / (2 * alpha) * (m_ * tr + a * sqrt_tr);
  }
  */

  // Static functions

  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static T m(const T &omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  /// Acentric factor
  T omega_;
  /// \f$ m = 0.3796 + 1.485 \omega - 0.1644 \omega^2 + 0.01667 \omega^3 \f$
  T m_;
};

/// @brief Makes Peng-Robinson EoS
/// @tparam T Value type
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename T>
inline peng_robinson_eos<T> make_pr_eos(const T &pc, const T &tc,
                                        const T &omega) {
  return {pc, tc, omega};
}

}  // namespace eos
