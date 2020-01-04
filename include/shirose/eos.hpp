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

#ifndef SHIROSE_EOS_HPP
#define SHIROSE_EOS_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "shirose/roots.hpp"

namespace shirose {

constexpr double gas_constant = 8.31446261815324;

/// @brief Van der Waals EoS policy for cubic_eos.
/// @tparam T Value type
template <typename T>
struct van_der_waals {
  /// Constant for attraction parameter
  static constexpr double omega_a = 0.421875;
  /// Constant for respulsion parameter
  static constexpr double omega_b = 0.125;

  // Static functions

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

  /// @brief Computes temperature correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  constexpr T alpha(const T &) const noexcept { return 1; }
};

/// @brief Soave-Redlich-Kwong EoS policy for cubic_eos.
/// @tparam T Value type
template <typename T>
class soave_redlich_kwong {
 public:
  /// Constant for attraction parameter
  static constexpr double omega_a = 0.42748;
  /// Constant for repulsion parameter
  static constexpr double omega_b = 0.08664;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / (v * (v + b));
  }

  /// @brief Computes the coefficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static std::array<T, 3> cubic_eq(const T &a, const T &b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  /// @brief Computes the fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static T fugacity_coeff(const T &z, const T &a, const T &b) noexcept {
    using std::exp;
    using std::log;
    return exp(z - 1 - log(z - b) - a * a / b * log(b / z + 1));
  }

  /// @brief Computes correction factor for temperature dependence of attraction
  /// parameter
  static T m(const T &omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  // Constructors

  /// @brief Constructs EoS
  /// @param[in] omega Acentric factor
  soave_redlich_kwong(const T &omega) noexcept : m_{this->m(omega)} {}

  /// @brief Computes correction factor for temperature dependence of attraction
  /// parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor for the attraction parameter
  T alpha(const T &tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

 private:
  /// Correction factor for temperature dependency of attraction parameter
  T m_;
};

/// @brief Peng-Robinson EoS policy for cubic_eos.
/// @tparam T Value type
template <typename T>
class peng_robinson {
 public:
  /// Constant for attraction parameter
  static constexpr double omega_a = 0.45724;
  /// Constant for repulsion parameter
  static constexpr double omega_b = 0.07780;
  /// Square root of 2 (two).
  static constexpr double sqrt2 = 1.4142135623730950488;
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
    return exp(z - 1 - log(z - b) -
               a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z - delta2 * b)));
  }

  /// @brief Correction factor for temperature dependency of attraction
  /// parameter
  /// @param[in] omega Acentric factor
  /// @returns Correction factor for temperature dependence of attraction
  /// parameter
  static T m(const T &omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  // Constructors

  /// @brief Constructs EoS
  /// @param[in] omega Acentric factor
  peng_robinson(const T &omega) : m_{this->m(omega)} {}

  // Member functions

  /// @brief Computes the temperature correction factor for the attraction
  /// parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor for the attraction parameter
  T alpha(const T &tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

 private:
  /// Correction factor for temperature dependency of attraction parameter
  T m_;
};

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam Policy EoS policy
///
/// Policy must have the following static functions:
///    - pressure(t, v, a, b)
///    - cubic_eq(ar, br)
///    - fugacity_coeff(z, ar, br)
/// where t is temperature, v is volume, a is attraction parameter, b is
/// repulsion parameter, ar is reduced attraction parameter, br is reduced
/// repulsion parameter, and z is z-factor.
///
/// In addition, the following member function is required:
///    - alpha(tr)
/// where tr is reduced temperature.
template <typename T, typename Policy>
class cubic_eos {
 public:
  // Static functions

  /// @brief Computes attraction parameter
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @returns Attraction parameter at critical condition
  static T attraction_param(const T &pc, const T &tc) noexcept {
    return (Policy::omega_a * gas_constant * gas_constant) * tc * tc / pc;
  }

  /// @brief Computes repulsion parameter
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @returns Repulsion parameter at critical condition
  static T repulsion_param(const T &pc, const T &tc) noexcept {
    return (Policy::omega_b * gas_constant) * tc / pc;
  }

  /// @brief Computes reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced attraction parameter
  static T reduced_attraction_param(const T &pr, const T &tr) noexcept {
    return Policy::omega_a * pr / (tr * tr);
  }

  /// @brief Computes reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced repulsion parameter
  static T reduced_repulsion_param(const T &pr, const T &tr) noexcept {
    return Policy::omega_b * pr / tr;
  }

  // Constructors

  /// @brief Constructs EoS.
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] policy EoS policy
  cubic_eos(const T &pc, const T &tc, const Policy &policy)
      : pc_{pc},
        tc_{tc},
        ac_{this->attraction_param(pc, tc)},
        bc_{this->repulsion_param(pc, tc)},
        policy_{policy} {}

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    const auto tr = t / tc_;
    const auto a = policy_.alpha(tr) * ac_;
    const auto b = bc_;
    return Policy::pressure(t, v, a, b);
  }

  /// @brief A state at given pressure and temperature expressed by this
  /// equation of state
  class pt_state {
   public:
    /// @brief Constructs a state
    /// @param[in] a Reduced attration parameter
    /// @param[in] b Reduced repulsion parameter
    pt_state(const T &a, const T &b) : a_{a}, b_{b} {}

    /// @brief Computes z-factor
    /// @returns An array of z-factors
    std::vector<T> zfactor() const noexcept {
      const auto p = Policy::cubic_eq(a_, b_);
      const auto x = roots(p);
      return real_roots(x);
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @returns Fugacity coefficient
    ///
    /// Z-factor must be one of solutions of zfactor() function.
    T fugacity_coeff(const T &z) const noexcept {
      return Policy::fugacity_coeff(z, a_, b_);
    }

   private:
    /// Reduced attraction parameter
    T a_;
    /// Reduced repulsion parameter
    T b_;
  };

  /// @brief Creates a state of given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns State
  pt_state state(const T &p, const T &t) noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto a = policy_.alpha(tr) * this->reduced_attraction_param(pr, tr);
    const auto b = this->reduced_repulsion_param(pr, tr);
    return {a, b};
  }

 private:
  T pc_;
  T tc_;
  T ac_;
  T bc_;
  Policy policy_;
};

/// @brief Van der Waals equation of state.
template <typename T>
using vdw_eos = cubic_eos<T, van_der_waals<T>>;

/// @brief Soave-Redlich-Kwong equation of state.
template <typename T>
using srk_eos = cubic_eos<T, soave_redlich_kwong<T>>;

/// @brief Peng-Robinson equation of state.
template <typename T>
using pr_eos = cubic_eos<T, peng_robinson<T>>;

/// @brief Makes van der Waals EoS
template <typename T>
vdw_eos<T> make_vdw_eos(const T &pc, const T &tc) {
  return {pc, tc, van_der_waals<T>{}};
}

/// @brief Makes Soave-Redlich-Kwong EoS
template <typename T>
srk_eos<T> make_srk_eos(const T &pc, const T &tc, const T &omega) {
  return {pc, tc, soave_redlich_kwong<T>{omega}};
}

/// @brief Makes Peng-Robinson EoS
template <typename T>
pr_eos<T> make_pr_eos(const T &pc, const T &tc, const T &omega) {
  return {pc, tc, peng_robinson<T>{omega}};
}

}  // namespace shirose

#endif  // SHIROSE_EOS_HPP