#pragma once

#include <array>  // std::array
#include <cmath>  // std::exp, std::log

#include "eos/cubic_eos.hpp"  // eos::cubic_eos_base

namespace eos {

template <typename T>
class van_der_waals_eos;

/// @brief Van der Waals EoS policy for CubicEos.
/// @tparam T Value type
template <typename T>
struct van_der_waals_eos_policy {
  using derived_type = van_der_waals_eos<T>;
  using value_type = T;

  static constexpr double omega_a = 0.421875;
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

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static T residual_enthalpy(const T &z, const T &t, const T &a, const T &b,
                             const T &beta) noexcept {
    return gas_constant * t * (z - 1 - a * (1 - beta) / z);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static T residual_entropy(const T &z, const T &a, const T &b,
                            const T &beta) noexcept {
    using std::log;
    return gas_constant * (log(z - b) + a * beta / z);
  }

  /*
  /// @brief Computes residual molar specific heat at constant volume
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] gamma Temperature correction factor
  static T residual_specific_heat_at_const_volume(const T &z, const T &a,
                                                  [[maybe_unused]] const T &b,
                                                  const T &gamma) noexcept {
    return gas_constant * gamma * a / z;
  }
  */
};

template <typename T>
class van_der_waals_eos : public cubic_eos_base<van_der_waals_eos_policy<T>> {
 public:
  using base_type = cubic_eos_base<van_der_waals_eos_policy<T>>;

  van_der_waals_eos() = default;

  van_der_waals_eos(const T &pc, const T &tc) noexcept : base_type{pc, tc} {}

  void set_params(const T &pc, const T &tc) noexcept {
    this->base_type::set_params(pc, tc);
  }

  T alpha(const T &) const noexcept { return 1.0; }

  T beta(const T &) const noexcept { return 0.0; }

  // T gamma(const T &) const noexcept { return 0.0; }
};

/// @brief Makes van der Waals EoS
/// @tparam T Value type
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
template <typename T>
inline van_der_waals_eos<T> make_vdw_eos(const T &pc, const T &tc) {
  return {pc, tc};
}

}  // namespace eos