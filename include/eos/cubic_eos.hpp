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

#include <vector>  // std::vector

#include "eos/constants.hpp"  // eos::gas_constant
#include "eos/roots.hpp"      // eos::real_roots

namespace eos {

template <typename T>
struct isothermal_state {
  T t;   /// Temperature
  T tr;  /// Reduced temperature
  T a;   /// Attraction parameter
  T b;   /// Repulsion parameter
};

template <typename T>
struct isobaric_isothermal_state {
  T p;   /// Pressure
  T t;   /// Temperature
  T pr;  /// Reduced pressure
  T tr;  /// Reduced temperature
  T ar;  /// Reduced attraction parameter
  T br;  /// Reduced repulsion parameter
};

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam Policy EoS policy
///
/// Policy must have the following types:
///    - derived_type: Derived EoS type
///    - value_type: Scalar type
///
/// Policy must have the following constants:
///    - omega_a: Constant for attraction parameter
///    - omega_b: Constant for repulsion parameter
///
/// Policy must have the following static functions:
///    - pressure(t, v, a, b)
///    - cubic_eq(ar, br)
///    - fugacity_coeff(z, ar, br)
///    - residual_enthalpy(z, t, ar, br, beta)
///    - residual_entropy(z, ar, br, beta)
/// where t is temperature, v is volume, a is attraction parameter, b is
/// repulsion parameter, ar is reduced attraction parameter, br is reduced
/// repulsion parameter, z is z-factor, and beta is a temperature correction
/// factor.
template <typename Policy>
class cubic_eos_base {
 public:
  using derived_type = typename Policy::derived_type;
  using value_type = typename Policy::value_type;

  static constexpr auto omega_a = Policy::omega_a;
  static constexpr auto omega_b = Policy::omega_b;

  // Constructors

  cubic_eos_base() = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  cubic_eos_base(const value_type &pc, const value_type &tc) noexcept
      : pc_{pc},
        tc_{tc},
        ac_{critical_attraction_param(pc, tc)},
        bc_{critical_repulsion_param(pc, tc)} {}

  // Member functions

  /// @brief Set parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void set_params(const value_type &pc, const value_type &tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    ac_ = critical_attraction_param(pc, tc);
    bc_ = critical_repulsion_param(pc, tc);
  }

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  isothermal_state<value_type> state(const value_type &t) const noexcept {
    const auto tr = t / tc_;
    const auto a = this->derived().alpha(tr) * ac_;
    const auto b = bc_;
    return {t, tr, a, b};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  isobaric_isothermal_state<value_type> state(const value_type &p,
                                              const value_type &t) const
      noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto ar =
        this->derived().alpha(tr) * reduced_attraction_param(pr, tr);
    const auto br = reduced_repulsion_param(pr, tr);
    return {p, t, pr, tr, ar, br};
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  value_type pressure(const value_type &t, const value_type &v) const noexcept {
    const auto s = this->state(t);
    return this->pressure(v, s);
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] v Volume
  /// @param[in] s Isothermal state
  value_type pressure(const value_type &v,
                      const isothermal_state<value_type> &s) const noexcept {
    return Policy::pressure(s.t, v, s.a, s.b);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  std::vector<value_type> zfactor(const value_type &p,
                                  const value_type &t) const noexcept {
    const auto s = this->state(p, t);
    return this->zfactor(s);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] s Isobaric-isothermal state
  /// @return A list of Z-factors
  std::vector<value_type> zfactor(
      const isobaric_isothermal_state<value_type> &s) const noexcept {
    return real_roots(Policy::cubic_eq(s.ar, s.br));
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] s Isobaric-isothermal state
  value_type fugacity_coeff(
      const value_type &z, const isobaric_isothermal_state<value_type> &s) const
      noexcept {
    return Policy::fugacity_coeff(z, s.ar, s.br);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] s Isobaric-isothermal state
  value_type residual_enthalpy(
      const value_type &z, const isobaric_isothermal_state<value_type> &s) const
      noexcept {
    const auto beta = this->derived().beta(s.tr);
    return Policy::residual_enthalpy(z, s.ar, s.br, beta);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] s Isobaric-isothermal state
  value_type residual_entropy(
      const value_type &z, const isobaric_isothermal_state<value_type> &s) const
      noexcept {
    const auto beta = this->derived().beta(s.tr);
    return Policy::residual_entropy(z, s.ar, s.br, beta);
  }

 private:
  // Static functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static value_type critical_attraction_param(const value_type &pc,
                                              const value_type &tc) noexcept {
    return (omega_a * gas_constant * gas_constant) * tc * tc / pc;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static value_type critical_repulsion_param(const value_type &pc,
                                             const value_type &tc) noexcept {
    return (omega_b * gas_constant) * tc / pc;
  }

  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static value_type reduced_attraction_param(const value_type &pr,
                                             const value_type &tr) noexcept {
    return omega_a * pr / (tr * tr);
  }

  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static value_type reduced_repulsion_param(const value_type &pr,
                                            const value_type &tr) noexcept {
    return omega_b * pr / tr;
  }

  // Member functions

  /// @brief Get reference to derived class object
  derived_type &derived() noexcept {
    return static_cast<derived_type &>(*this);
  }

  /// @brief Get const reference to derived class object
  const derived_type &derived() const noexcept {
    return static_cast<const derived_type &>(*this);
  }

  value_type pc_;  /// Critical pressure
  value_type tc_;  /// Critical temperature
  value_type ac_;  /// Critical attraction parameter
  value_type bc_;  /// Critical repulsion parameter
};

}  // namespace eos