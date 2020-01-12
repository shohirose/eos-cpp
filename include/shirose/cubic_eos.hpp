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

#ifndef SHIROSE_CUBIC_EOS_HPP
#define SHIROSE_CUBIC_EOS_HPP

#include <vector>  // std::vector

#include "shirose/constants.hpp"  // shirose::gas_constant
#include "shirose/roots.hpp"      // shirose::roots

namespace shirose {

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam Eos EoS policy
/// @tparam CorrectionPolicy Temperature correction policy for attraction
/// parameter
///
/// Eos must have the following static functions:
///    - pressure(t, v, a, b)
///    - cubic_eq(ar, br)
///    - fugacity_coeff(z, ar, br)
///    - residual_enthalpy(z, t, ar, br, beta)
///    - residual_entropy(z, ar, br, beta)
/// where t is temperature, v is volume, a is attraction parameter, b is
/// repulsion parameter, ar is reduced attraction parameter, br is reduced
/// repulsion parameter, z is z-factor, and beta is a temperature correction
/// factor.
///
/// CorrectionPolicy must have the following member functions:
///    - alpha(tr): temperature correction factor for attraction parameter
///    - beta(tr): \f$ \beta = \frac{T_r}{\alpha} \frac{\mathrm{d}
///    \alpha}{\mathrm{d} T_r} \f$
/// where tr is reduced temperature.
template <typename T, typename Eos, typename CorrectionPolicy>
class cubic_eos {
 public:
  // Static functions

  /// @brief Computes attraction parameter
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @returns Attraction parameter at critical condition
  static T attraction_param(const T &pc, const T &tc) noexcept {
    return (Eos::omega_a * gas_constant * gas_constant) * tc * tc / pc;
  }

  /// @brief Computes repulsion parameter
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @returns Repulsion parameter at critical condition
  static T repulsion_param(const T &pc, const T &tc) noexcept {
    return (Eos::omega_b * gas_constant) * tc / pc;
  }

  /// @brief Computes reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced attraction parameter
  static T reduced_attraction_param(const T &pr, const T &tr) noexcept {
    return Eos::omega_a * pr / (tr * tr);
  }

  /// @brief Computes reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced repulsion parameter
  static T reduced_repulsion_param(const T &pr, const T &tr) noexcept {
    return Eos::omega_b * pr / tr;
  }

  // Constructors

  /// @brief Constructs EoS.
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] policy Temperature correction policy for attraction parameter
  cubic_eos(const T &pc, const T &tc, const CorrectionPolicy &policy)
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
    return Eos::pressure(t, v, a, b);
  }

  /// @brief A state at given pressure and temperature expressed by this
  /// equation of state
  class pt_state {
   public:
    /// @brief Constructs a state
    /// @param[in] t Temperature
    /// @param[in] ar Reduced attration parameter
    /// @param[in] br Reduced repulsion parameter
    /// @param[in] beta Temperature correction factor
    /// @param[in] gamma Temperature correction factor
    pt_state(const T &t, const T &ar, const T &br, const T &beta,
             const T &gamma)
        : t_{t}, ar_{ar}, br_{br}, beta_{beta}, gamma_{gamma} {}

    /// @brief Computes z-factor
    /// @return An array of z-factors
    std::vector<T> zfactor() const noexcept {
      const auto p = Eos::cubic_eq(ar_, br_);
      const auto x = roots(p);
      return real_roots(x);
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @return Fugacity coefficient
    T fugacity_coeff(const T &z) const noexcept {
      return Eos::fugacity_coeff(z, ar_, br_);
    }

    /// @brief Computes residual enthalpy
    /// @param[in] z Z-factor
    /// @return Residual enthalpy
    T residual_enthalpy(const T &z) const noexcept {
      return Eos::residual_enthalpy(z, t_, ar_, br_, beta_);
    }

    /// @brief Computes residual entropy
    /// @param[in] z Z-factor
    /// @return Residual entropy
    T residual_entropy(const T &z) const noexcept {
      return Eos::residual_entropy(z, ar_, br_, beta_);
    }

    /// @brief Computes residual molar specifiec heat at constant volume
    /// @param[in] z Z-factor
    T residual_specific_heat_v(const T &z) const noexcept {
      return Eos::residual_specific_heat_v(z, ar_, br_, gamma_);
    }

   private:
    /// Temperature
    T t_;
    /// Reduced attraction parameter
    T ar_;
    /// Reduced repulsion parameter
    T br_;
    /// Temperature correction factor
    T beta_;
    /// Temperature correction factor
    T gamma_;
  };

  /// @brief Creates a state of given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns State
  pt_state state(const T &p, const T &t) const noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto ar = policy_.alpha(tr) * this->reduced_attraction_param(pr, tr);
    const auto br = this->reduced_repulsion_param(pr, tr);
    const auto beta = policy_.beta(tr);
    const auto gamma = policy_.gamma(tr);
    return {t, ar, br, beta, gamma};
  }

 private:
  /// Critical pressure
  T pc_;
  /// Critical temperature
  T tc_;
  /// Attraction parameter at critical condition
  T ac_;
  /// Repulsion parameter at critical condition
  T bc_;
  /// Temperature correction policy for attraction parameter
  CorrectionPolicy policy_;
};

}  // namespace shirose

#endif  // SHIROSE_CUBIC_EOS_HPP