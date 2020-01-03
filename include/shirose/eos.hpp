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

/// @brief Base type of cubic equation of state policy.
/// @tparam T Value type
template <typename T>
class cubic_eos_policy_base {
 public:
  /// Value type
  using value_type = T;

  // Constructors

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] pc Critical attraction parameter
  /// @param[in] tc Critical repulsion parameter
  cubic_eos_policy_base(const T &pc, const T &tc, const T &ac,
                        const T &bc) noexcept
      : pc_{pc}, tc_{tc}, ac_{ac}, bc_{bc} {}

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  T reduced_pressure(const T &p) const noexcept { return p / pc_; }

  /// @brief Computes reduced pressure
  /// @param[in] t Temperature
  T reduced_temperature(const T &t) const noexcept { return t / tc_; }

  /// @brief Computes attraction parameter
  T attraction_param(const T &) const noexcept { return ac_; }

  /// @brief Computes repulsion parameter
  T repulsion_param(const T &) const noexcept { return bc_; }

 protected:
  /// Critical pressure
  T pc_;
  /// Critical temperature
  T tc_;
  /// Critical attraction parameter
  T ac_;
  /// Critical repulsion parameter
  T bc_;
};

/// @brief Van der Waals EoS policy for cubic_eos.
/// @tparam T Value type
template <typename T>
class van_der_waals : public cubic_eos_policy_base<T> {
 public:
  using base_type = cubic_eos_policy_base<T>;
  using value_type = T;

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
  /// @returns An array of Z-factor
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

  // Constructors

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  van_der_waals(const T &pc, const T &tc) noexcept
      : base_type(pc, tc,
                  (omega_a * gas_constant * gas_constant) * tc * tc / pc,
                  (omega_b * gas_constant) * tc / pc) {}

  /// @brief Computes reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  T reduced_attraction_param(const T &pr, const T &tr) const noexcept {
    return omega_a * pr / (tr * tr);
  }

  /// @brief Computes reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  T reduced_repulsion_param(const T &pr, const T &tr) const noexcept {
    return omega_b * pr / tr;
  }

 private:
};

/// @brief Soave-Redlich-Kwong EoS policy for cubic_eos.
/// @tparam T Value type
template <typename T>
class soave_redlich_kwong : public cubic_eos_policy_base<T> {
 public:
  using base_type = cubic_eos_policy_base<T>;
  using value_type = T;

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

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  soave_redlich_kwong(const T &pc, const T &tc, const T &omega) noexcept
      : base_type(pc, tc,
                  (omega_a * gas_constant * gas_constant) * tc * tc / pc,
                  (omega_b * gas_constant) * tc / pc),
        m_{this->m(omega)} {}

  /// @brief Computes attraction parameter
  T attraction_param(const T &tr) const noexcept {
    return this->alpha(tr) * this->ac_;
  }

  /// @brief Computes reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  T reduced_attraction_param(const T &pr, const T &tr) const noexcept {
    return this->alpha(tr) * omega_a * pr / (tr * tr);
  }

  /// @brief Computes reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  T reduced_repulsion_param(const T &pr, const T &tr) const noexcept {
    return omega_b * pr / tr;
  }

 private:
  /// @brief Computes correction factor for temperature dependence of attraction
  /// parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor for the attraction parameter
  T alpha(const T &tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// Correction factor for temperature dependency of attraction parameter
  T m_;
};

/// @brief Peng-Robinson EoS policy for cubic_eos.
/// @tparam T Value type
///
/// This defines function for Peng-Robinson equation of state.
template <typename T>
class peng_robinson : public cubic_eos_policy_base<T> {
 public:
  /// Base type
  using base_type = cubic_eos_policy_base<T>;
  /// Value type
  using value_type = T;

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
  /// @returns An array of coefficients of the cubic equation.
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

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  peng_robinson(const T &pc, const T &tc, const T &omega)
      : base_type(pc, tc,
                  (omega_a * gas_constant * gas_constant) * tc * tc / pc,
                  (omega_b * gas_constant) * tc / pc),
        m_{this->m(omega)} {}

  // Member functions

  /// @brief Computes attraction parameter
  /// @param[in] tr Reduced temperature
  /// @returns Attraction parameter
  T attraction_param(const T &tr) const noexcept {
    return this->alpha(tr) * this->ac_;
  }

  /// @brief Computes reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  T reduced_attraction_param(const T &pr, const T &tr) const noexcept {
    return this->alpha(tr) * omega_a * pr / (tr * tr);
  }

  /// @brief Computes reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  T reduced_repulsion_param(const T &pr, const T &tr) const noexcept {
    return omega_b * pr / tr;
  }

 private:
  /// @brief Computes the temperature correction factor for the attraction
  /// parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor for the attraction parameter
  T alpha(const T &tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// Correction factor for temperature dependency of attraction parameter
  T m_;
};

/// @brief Cubic equation of state
/// @tparam Policy Policy of cubic EoS
///
/// Policy must have the following static functions:
///    - pressure(a, b)
///    - cubic_eq(ar, br)
///    - zfactor(ar, br)
///    - fugacity_coeff(z, ar, br)
/// where a is attraction parameter, b is repulsion parameter, ar is reduced
/// attraction parameter, br is reduced repulsion parameter, and z is z-factor.
///
/// In addition, the following member functions are required:
///    - reduced_pressure(p)
///    - reduced_temperature(t)
///    - attraction_param(tr)
///    - repulsion_param(tr)
///    - reduced_attraction_param(pr, tr)
///    - reduced_repulsion_param(pr, tr)
/// where p is pressure, t is temperature, pr is reduced pressure, tr is reduced
/// temperature.
template <typename Policy>
class cubic_eos {
 public:
  using value_type = typename Policy::value_type;

  /// @brief Constructs EoS.
  cubic_eos(const Policy &policy) : policy_{policy} {}

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  value_type pressure(const value_type &t, const value_type &v) noexcept {
    const auto tr = policy_.reduced_temperature(t);
    const auto a = policy_.attraction_param(tr);
    const auto b = policy_.repulsion_param(tr);
    return Policy::pressure(t, v, a, b);
  }

  /// @brief A state at given pressure and temperature expressed by this
  /// equation of state
  class pt_state {
   public:
    /// @brief Constructs a state
    /// @param[in] a Reduced attration parameter
    /// @param[in] b Reduced repulsion parameter
    pt_state(const value_type &a, const value_type &b) : a_{a}, b_{b} {}

    /// @brief Computes z-factor
    /// @returns An array of z-factors
    std::vector<value_type> zfactor() const noexcept {
      const auto p = Policy::cubic_eq(a_, b_);
      const auto x = roots(p);
      return real_roots(x);
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @returns Fugacity coefficient
    ///
    /// Z-factor must be one of solutions of zfactor() function.
    value_type fugacity_coeff(const value_type &z) const noexcept {
      return Policy::fugacity_coeff(z, a_, b_);
    }

   private:
    /// Reduced attraction parameter
    value_type a_;
    /// Reduced repulsion parameter
    value_type b_;
  };

  /// @brief Creates a state of given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns State
  pt_state state(const value_type &p, const value_type &t) noexcept {
    const auto pr = policy_.reduced_pressure(p);
    const auto tr = policy_.reduced_temperature(t);
    const auto a = policy_.reduced_attraction_param(pr, tr);
    const auto b = policy_.reduced_repulsion_param(pr, tr);
    return {a, b};
  }

 private:
  Policy policy_;
};

/// @brief Van der Waals equation of state.
template <typename T>
using vdw_eos = cubic_eos<van_der_waals<T>>;

/// @brief Soave-Redlich-Kwong equation of state.
template <typename T>
using srk_eos = cubic_eos<soave_redlich_kwong<T>>;

/// @brief Peng-Robinson equation of state.
template <typename T>
using pr_eos = cubic_eos<peng_robinson<T>>;

template <typename Policy>
cubic_eos<Policy> make_cubic_eos(Policy &&policy) {
  return cubic_eos<Policy>(std::forward<Policy>(policy));
}

}  // namespace shirose

#endif  // SHIROSE_EOS_HPP