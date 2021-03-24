#pragma once

#include <vector>  // std::vector

#include "eos/constants.hpp"  // eos::gas_constant
#include "eos/roots.hpp"      // eos::real_roots

namespace eos {

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

  class isothermal_state {
   public:
    /// @param[in] t Temperature
    /// @param[in] tr Reduced temperature
    /// @param[in] a Attraction parameter
    /// @param[in] b Repulsion parameter
    isothermal_state(const value_type &t, const value_type &tr,
                     const value_type &a, const value_type &b) noexcept
        : t_{t}, tr_{tr}, a_{a}, b_{b} {}

    /// @brief Computes pressure at given temperature and volume
    /// @param[in] v Volume
    value_type pressure(const value_type &v) const noexcept {
      return Policy::pressure(t_, v, a_, b_);
    }

   private:
    value_type t_;   /// Temperature
    value_type tr_;  /// Reduced temperature
    value_type a_;   /// Attraction parameter
    value_type b_;   /// Repulsion parameter
  };

  class isobaric_isothermal_state {
   public:
    isobaric_isothermal_state(const value_type &p, const value_type &t,
                              const value_type &pr, const value_type &tr,
                              const value_type &ar, const value_type &br,
                              const value_type &beta) noexcept
        : p_{p}, t_{t}, pr_{pr}, tr_{tr}, ar_{ar}, br_{br}, beta_{beta} {}

    /// @brief Computes Z-factor at given pressure and temperature
    /// @param[in] s Isobaric-isothermal state
    /// @return A list of Z-factors
    std::vector<value_type> zfactor() const noexcept {
      return real_roots(Policy::cubic_eq(ar_, br_));
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @param[in] s Isobaric-isothermal state
    value_type fugacity_coeff(const value_type &z) const noexcept {
      return Policy::fugacity_coeff(z, ar_, br_);
    }

    /// @brief Computes residual enthalpy
    /// @param[in] z Z-factor
    /// @param[in] s Isobaric-isothermal state
    value_type residual_enthalpy(const value_type &z) const noexcept {
      return Policy::residual_enthalpy(z, ar_, br_, beta_);
    }

    /// @brief Computes residual entropy
    /// @param[in] z Z-factor
    /// @param[in] s Isobaric-isothermal state
    value_type residual_entropy(const value_type &z) const noexcept {
      return Policy::residual_entropy(z, ar_, br_, beta_);
    }

   private:
    value_type p_;   /// Pressure
    value_type t_;   /// Temperature
    value_type pr_;  /// Reduced pressure
    value_type tr_;  /// Reduced temperature
    value_type ar_;  /// Reduced attraction parameter
    value_type br_;  /// Reduced repulsion parameter
    value_type beta_;
  };

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  isothermal_state state(const value_type &t) const noexcept {
    const auto tr = t / tc_;
    const auto a = this->attraction_param(tr);
    const auto b = this->repulsion_param();
    return {t, tr, a, b};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  isobaric_isothermal_state state(const value_type &p,
                                  const value_type &t) const noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto ar = this->reduced_attraction_param(pr, tr);
    const auto br = this->reduced_repulsion_param(pr, tr);
    const auto beta = this->derived().beta(tr);
    return {p, t, pr, tr, ar, br, beta};
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  value_type pressure(const value_type &t, const value_type &v) const noexcept {
    return this->state(t).pressure(v);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  std::vector<value_type> zfactor(const value_type &p,
                                  const value_type &t) const noexcept {
    return this->state(p, t).zfactor();
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

  // Member functions

  /// @brief Returns attraction parameter at a given temperature.
  /// @param[in] tr Reduced temperature
  value_type attraction_param(const value_type &tr) const noexcept {
    return this->derived().alpha(tr) * ac_;
  }

  /// @brief Returns repulsion parameter.
  value_type repulsion_param() const noexcept { return bc_; }

  /// @brief Returns reduced attraction parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  value_type reduced_attraction_param(const value_type &pr,
                                      const value_type &tr) const noexcept {
    return this->derived().alpha(tr) * omega_a * pr / (tr * tr);
  }

  /// @brief Returns reduced repulsion parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  value_type reduced_repulsion_param(const value_type &pr,
                                     const value_type &tr) const noexcept {
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