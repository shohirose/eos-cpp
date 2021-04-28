#pragma once

#include <vector>  // std::vector

#include "eos/constants.hpp"        // eos::gas_constant
#include "eos/math/polynomial.hpp"  // eos::real_roots

namespace eos {

namespace detail {

template <typename Eos>
struct cubic_eos_traits {};

}  // namespace detail

template <typename Eos>
class isothermal_line {
 public:
  /// @param[in] t Temperature
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  isothermal_line(double t, double a, double b) noexcept
      : t_{t}, a_{a}, b_{b} {}

  isothermal_line() = default;
  isothermal_line(const isothermal_line &) = default;
  isothermal_line(isothermal_line &&) = default;

  isothermal_line &operator=(const isothermal_line &) = default;
  isothermal_line &operator=(isothermal_line &&) = default;

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] v Volume
  double pressure(double v) const noexcept {
    return Eos::pressure_impl(t_, v, a_, b_);
  }

 private:
  double t_;  /// Temperature
  double a_;  /// Attraction parameter
  double b_;  /// Repulsion parameter
};

template <typename Eos>
class isobaric_isothermal_state {
 public:
  isobaric_isothermal_state(double ar, double br, double beta) noexcept
      : ar_{ar}, br_{br}, beta_{beta} {}

  isobaric_isothermal_state() = default;
  isobaric_isothermal_state(const isobaric_isothermal_state &) = default;
  isobaric_isothermal_state(isobaric_isothermal_state &&) = default;

  isobaric_isothermal_state &operator=(const isobaric_isothermal_state &) =
      default;
  isobaric_isothermal_state &operator=(isobaric_isothermal_state &&) = default;

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] s Isobaric-isothermal state
  /// @return A list of Z-factors
  std::vector<double> zfactor() const noexcept {
    const auto a = Eos::zfactor_cubic_eq_impl(ar_, br_);
    return real_roots(a[0], a[1], a[2]);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  double fugacity_coeff(double z) const noexcept {
    return Eos::fugacity_coeff_impl(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  double residual_enthalpy(double z) const noexcept {
    return Eos::residual_enthalpy_impl(z, ar_, br_, beta_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  double residual_entropy(double z) const noexcept {
    return Eos::residual_entropy_impl(z, ar_, br_, beta_);
  }

 private:
  double ar_;    /// Reduced attraction parameter
  double br_;    /// Reduced repulsion parameter
  double beta_;  /// The derivative of temperature correction factor for
                 /// attraction parameter
};

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam Derived Concrete EoS class
///
/// Derived EoS classes must have the following static functions:
///    - pressure_impl(t, v, a, b)
///    - zfactor_cubic_eq_impl(ar, br)
///    - fugacity_coeff_impl(z, ar, br)
///    - residual_enthalpy_impl(z, t, ar, br, beta)
///    - residual_entropy_impl(z, ar, br, beta)
///    - alpha(tr)
///    - beta()
/// where t is temperature, v is volume, a is attraction parameter, b is
/// repulsion parameter, ar is reduced attraction parameter, br is reduced
/// repulsion parameter, z is z-factor, and beta is a temperature correction
/// factor.
///
/// cubic_eos_traits class specialized for each concrete EoS class must be
/// defined in the detail namespace. cubic_eos_traits class must define the
/// following types and constants:
///    - double: double type
///    - omega_a: Constant for attraction parameter
///    - omega_b: Constant for repulsion parameter
///
template <typename Derived>
class cubic_eos_base {
 public:
  static constexpr auto omega_a = detail::cubic_eos_traits<Derived>::omega_a;
  static constexpr auto omega_b = detail::cubic_eos_traits<Derived>::omega_b;

  // Constructors

  cubic_eos_base() = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  cubic_eos_base(double pc, double tc) noexcept
      : pc_{pc},
        tc_{tc},
        ac_{this->critical_attraction_param(pc, tc)},
        bc_{this->critical_repulsion_param(pc, tc)} {}

  cubic_eos_base(const cubic_eos_base &) = default;
  cubic_eos_base(cubic_eos_base &&) = default;

  cubic_eos_base &operator=(const cubic_eos_base &) = default;
  cubic_eos_base &operator=(cubic_eos_base &&) = default;

  // Member functions

  /// @brief Set parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void set_params(double pc, double tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    ac_ = this->critical_attraction_param(pc, tc);
    bc_ = this->critical_repulsion_param(pc, tc);
  }

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  isothermal_line<Derived> create_isothermal_state(double t) const noexcept {
    const auto tr = this->reduced_temperature(t);
    const auto a = this->attraction_param(tr);
    const auto b = this->repulsion_param();
    return {t, a, b};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  isobaric_isothermal_state<Derived> create_isobaric_isothermal_state(
      double p, double t) const noexcept {
    const auto pr = this->reduced_pressure(p);
    const auto tr = this->reduced_temperature(t);
    const auto ar = this->reduced_attraction_param(pr, tr);
    const auto br = this->reduced_repulsion_param(pr, tr);
    const auto beta = this->derived().beta(tr);
    return {ar, br, beta};
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  double pressure(double t, double v) const noexcept {
    const auto tr = this->reduced_temperature(t);
    const auto a = this->attraction_param(tr);
    const auto b = this->repulsion_param();
    return Derived::pressure_impl(t, v, a, b);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  std::vector<double> zfactor(double p, double t) const noexcept {
    return this->create_isobaric_isothermal_state(p, t).zfactor();
  }

 private:
  // Static functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static double critical_attraction_param(double pc, double tc) noexcept {
    constexpr auto R = gas_constant<double>();
    return (omega_a * R * R) * tc * tc / pc;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static double critical_repulsion_param(double pc, double tc) noexcept {
    constexpr auto R = gas_constant<double>();
    return (omega_b * R) * tc / pc;
  }

  // Member functions

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  double reduced_pressure(double p) const noexcept { return p / pc_; }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature
  double reduced_temperature(double t) const noexcept { return t / tc_; }

  /// @brief Returns attraction parameter at a given temperature.
  /// @param[in] tr Reduced temperature
  double attraction_param(double tr) const noexcept {
    return this->derived().alpha(tr) * ac_;
  }

  /// @brief Returns repulsion parameter.
  double repulsion_param() const noexcept { return bc_; }

  /// @brief Returns reduced attraction parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  double reduced_attraction_param(double pr, double tr) const noexcept {
    return this->derived().alpha(tr) * omega_a * pr / (tr * tr);
  }

  /// @brief Returns reduced repulsion parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  double reduced_repulsion_param(double pr, double tr) const noexcept {
    return omega_b * pr / tr;
  }

  /// @brief Get reference to derived class object
  Derived &derived() noexcept { return static_cast<Derived &>(*this); }

  /// @brief Get const reference to derived class object
  const Derived &derived() const noexcept {
    return static_cast<const Derived &>(*this);
  }

  double pc_;  /// Critical pressure
  double tc_;  /// Critical temperature
  double ac_;  /// Critical attraction parameter
  double bc_;  /// Critical repulsion parameter
};

}  // namespace eos