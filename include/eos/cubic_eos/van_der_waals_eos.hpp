#pragma once

#include <array>  // std::array
#include <cmath>  // std::exp, std::log

#include "eos/cubic_eos/cubic_eos_base.hpp"  // eos::cubic_eos_base

namespace eos {

class van_der_waals_eos;

namespace detail {

template <>
struct cubic_eos_traits<van_der_waals_eos> {
  static constexpr auto omega_a = 0.421875;
  static constexpr auto omega_b = 0.125;
};

}  // namespace detail

/// @brief Van der Waals Equations of State
class van_der_waals_eos : public cubic_eos_base<van_der_waals_eos> {
 public:
  using base_type = cubic_eos_base<van_der_waals_eos>;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static double pressure(double t, double v, double a, double b) noexcept {
    return gas_constant<double>() * t / (v - b) - a / (v * v);
  }

  /// @brief Computes coefficients of the cubic equation of Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<double, 3> zfactor_cubic_eq(double a,
                                                     double b) noexcept {
    return {-b - 1, a, -a * b};
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static double fugacity_coeff(double z, double a, double b) noexcept {
    return std::exp(-std::log(z - b) - a / z + z - 1);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static double residual_enthalpy(double z, double t, double a,
                                       [[maybe_unused]] double b,
                                       double beta) noexcept {
    return gas_constant<double>() * t * (z - 1 - a * (1 - beta) / z);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static double residual_entropy(double z, double a, double b,
                                      double beta) noexcept {
    return gas_constant<double>() * (std::log(z - b) + a * beta / z);
  }

  /*
  /// @brief Computes residual molar specific heat at constant volume
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] gamma Temperature correction factor
  static double residual_specific_heat_at_const_volume(double z, const
  double &a,
                                                  [[maybe_unused]] const double
  &b, double gamma) noexcept { return gas_constant * gamma * a / z;
  }
  */

  van_der_waals_eos() = default;

  van_der_waals_eos(double pc, double tc) noexcept : base_type{pc, tc} {}

  van_der_waals_eos(const van_der_waals_eos &) = default;
  van_der_waals_eos(van_der_waals_eos &&) = default;

  van_der_waals_eos &operator=(const van_der_waals_eos &) = default;
  van_der_waals_eos &operator=(van_der_waals_eos &&) = default;

  void set_params(double pc, double tc) noexcept {
    this->base_type::set_params(pc, tc);
  }

  constexpr double alpha(double) const noexcept { return 1.0; }

  constexpr double beta(double) const noexcept { return 0.0; }

  // double gamma(double ) const noexcept { return 0.0; }
};

/// @brief Makes van der Waals EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
inline van_der_waals_eos make_van_der_waals_eos(double pc, double tc) {
  return {pc, tc};
}

}  // namespace eos