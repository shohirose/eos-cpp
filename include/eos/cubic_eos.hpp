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

/// @brief Base class of two-parameter cubic EoS
/// @tparam T Value type
/// @tparam Eos EoS policy
///
/// Eos must define constants, `omega_a` and `omega_b`.
template <typename T, typename Eos>
class CubicEosBase {
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
};

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam T Value type
/// @tparam Eos EoS policy
/// @tparam Corrector Temperature corrector for attraction parameter
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
/// Corrector except DefaultCorrector must have the following member functions:
///    - alpha(tr): temperature correction factor for attraction parameter
///    - beta(tr): \f$ \beta = \frac{T_r}{\alpha} \frac{\mathrm{d}
///    \alpha}{\mathrm{d} T_r} \f$
///    - gamma(tr): \f$ \gamma = \frac{T_r^2}{\alpha} \frac{\mathrm{d}^2
///    \alpha}{\mathrm{d} T_r^2} \f$
/// where tr is reduced temperature.
template <typename T, typename Eos, typename Corrector>
class CubicEos : public CubicEosBase<T, Eos> {
 public:
  // Constructors

  /// @brief Constructs EoS.
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] corrector Temperature corrector for attraction parameter
  CubicEos(const T &pc, const T &tc, const Corrector &corrector)
      : pc_{pc},
        tc_{tc},
        ac_{this->attraction_param(pc, tc)},
        bc_{this->repulsion_param(pc, tc)},
        corrector_{corrector} {}

  /// @brief Constructs EoS.
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] corrector Temperature corrector for attraction parameter
  CubicEos(const T &pc, const T &tc, Corrector &&corrector)
      : pc_{pc},
        tc_{tc},
        ac_{this->attraction_param(pc, tc)},
        bc_{this->repulsion_param(pc, tc)},
        corrector_{std::move(corrector)} {}

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    const auto tr = t / tc_;
    const auto a = corrector_.alpha(tr) * ac_;
    const auto b = bc_;
    return Eos::pressure(t, v, a, b);
  }

  /// @brief Isothermal state of a cubic EoS
  class IsothermalState {
   public:
    /// @brief Constructs isothermal state
    /// @param[in] t Temperature
    /// @param[in] a Attraction parameter
    /// @param[in] b Repulsion parameter
    IsothermalState(const T &t, const T &a, const T &b) : t_{t}, a_{a}, b_{b} {}

    /// @brief Computes pressure at given volume along this isothermal line
    /// @param[in] v Volume
    /// @return Pressure
    T pressure(const T &v) const noexcept {
      return Eos::pressure(t_, v, a_, b_);
    }

   private:
    T t_;  /// Temperature
    T a_;  /// Attraction parameter
    T b_;  /// Repulsion parameter
  };

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  /// @return Isothermal state
  IsothermalState state(const T &t) const noexcept {
    const auto tr = t / tc_;
    const auto a = corrector_.alpha(tr) * ac_;
    const auto b = bc_;
    return {t, a, b};
  }

  /// @brief Computes pressure at given temperature and pressure
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @return Pressure
  T pressure(const T &t, const T &v) const noexcept {
    const auto tr = t / tc_;
    const auto a = corrector_.alpha(tr) * ac_;
    const auto b = bc_;
    return Eos::pressure(t, v, a, b);
  }

  /// @brief A state at given pressure and temperature expressed by this
  /// equation of state
  class IsobaricIsothermalState {
   public:
    /// @brief Constructs a state
    /// @param[in] p Pressure
    /// @param[in] t Temperature
    /// @param[in] ar Reduced attration parameter
    /// @param[in] br Reduced repulsion parameter
    /// @param[in] beta Temperature correction factor
    /// @param[in] gamma Temperature correction factor
    IsobaricIsothermalState(const T &p, const T &t, const T &ar, const T &br,
                            const T &beta, const T &gamma)
        : p_{p}, t_{t}, ar_{ar}, br_{br}, beta_{beta}, gamma_{gamma} {}

    /// @brief Computes z-factor
    /// @return An array of z-factors
    std::vector<T> zfactor() const noexcept {
      const auto p = Eos::cubic_eq(ar_, br_);
      return real_roots(p);
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
    /// Pressure
    T p_;
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
  IsobaricIsothermalState state(const T &p, const T &t) const noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto ar =
        corrector_.alpha(tr) * this->reduced_attraction_param(pr, tr);
    const auto br = this->reduced_repulsion_param(pr, tr);
    const auto beta = corrector_.beta(tr);
    const auto gamma = corrector_.gamma(tr);
    return {p, t, ar, br, beta, gamma};
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
  Corrector corrector_;
};

/// @brief Partial specialization of CubicEos for DefaultCorrector.
template <typename T, typename Eos>
class CubicEos<T, Eos, DefaultCorrector<T>> : public CubicEosBase<T, Eos> {
 public:
  // Constructors

  /// @brief Constructs EoS.
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEos(const T &pc, const T &tc)
      : pc_{pc},
        tc_{tc},
        ac_{this->attraction_param(pc, tc)},
        bc_{this->repulsion_param(pc, tc)} {}

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    return Eos::pressure(t, v, ac_, bc_);
  }

  /// @brief Isothermal state of a cubic EoS
  class IsothermalState {
   public:
    /// @brief Constructs isothermal state
    /// @param[in] t Temperature
    /// @param[in] a Attraction parameter
    /// @param[in] b Repulsion parameter
    IsothermalState(const T &t, const T &a, const T &b) : t_{t}, a_{a}, b_{b} {}

    /// @brief Computes pressure at given volume along this isothermal line
    /// @param[in] v Volume
    /// @return Pressure
    T pressure(const T &v) const noexcept {
      return Eos::pressure(t_, v, a_, b_);
    }

   private:
    T t_;  /// Temperature
    T a_;  /// Attraction parameter
    T b_;  /// Repulsion parameter
  };

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  /// @return Isothermal state
  IsothermalState state(const T &t) const noexcept { return {t, ac_, bc_}; }

  /// @brief Computes pressure at given temperature and pressure
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @return Pressure
  T pressure(const T &t, const T &v) const noexcept {
    return Eos::pressure(t, v, ac_, bc_);
  }

  /// @brief A state at given pressure and temperature expressed by this
  /// equation of state
  class IsobaricIsothermalState {
   public:
    /// @brief Constructs a state
    /// @param[in] p Pressure
    /// @param[in] t Temperature
    /// @param[in] ar Reduced attration parameter
    /// @param[in] br Reduced repulsion parameter
    IsobaricIsothermalState(const T &p, const T &t, const T &ar, const T &br)
        : p_{p}, t_{t}, ar_{ar}, br_{br} {}

    /// @brief Computes z-factor
    /// @return An array of z-factors
    std::vector<T> zfactor() const noexcept {
      const auto p = Eos::cubic_eq(ar_, br_);
      return real_roots(p);
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
      return Eos::residual_enthalpy(z, t_, ar_, br_, 0.0);
    }

    /// @brief Computes residual entropy
    /// @param[in] z Z-factor
    /// @return Residual entropy
    T residual_entropy(const T &z) const noexcept {
      return Eos::residual_entropy(z, ar_, br_, 0.0);
    }

    /// @brief Computes residual molar specifiec heat at constant volume
    /// @param[in] z Z-factor
    T residual_specific_heat_v(const T &z) const noexcept {
      return Eos::residual_specific_heat_v(z, ar_, br_, 0.0);
    }

   private:
    /// Pressure
    T p_;
    /// Temperature
    T t_;
    /// Reduced attraction parameter
    T ar_;
    /// Reduced repulsion parameter
    T br_;
  };

  /// @brief Creates a state of given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns State
  IsobaricIsothermalState state(const T &p, const T &t) const noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto ar = this->reduced_attraction_param(pr, tr);
    const auto br = this->reduced_repulsion_param(pr, tr);
    return {p, t, ar, br};
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
};

}  // namespace eos