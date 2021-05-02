#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos/cubic_eos_base.hpp"  // eos::cubic_eos_base

namespace eos {

class peng_robinson_eos;

namespace detail {

template <>
struct cubic_eos_traits<peng_robinson_eos> {
  static constexpr auto omega_a = 0.45724;
  static constexpr auto omega_b = 0.07780;
};

}  // namespace detail

/// @brief Peng-Robinson EoS.
class peng_robinson_eos : public cubic_eos_base<peng_robinson_eos> {
 public:
  using base_type = cubic_eos_base<peng_robinson_eos>;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static double pressure(double t, double v, double a, double b) noexcept {
    constexpr auto R = gas_constant<double>();
    return R * t / (v - b) - a / (v * (v + b) + b * (v - b));
  }

  /// @brief Computes coeficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static std::array<double, 3> zfactor_cubic_eq(double a, double b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static double ln_fugacity_coeff(double z, double a, double b) noexcept {
    return z - 1 - std::log(z - b) - q(z, a, b);
  }

  /// @brief Computes a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static double fugacity_coeff(double z, double a, double b) noexcept {
    return std::exp(ln_fugacity_coeff(z, a, b));
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static double residual_enthalpy(double z, double t, double a, double b,
                                  double beta) noexcept {
    constexpr auto R = gas_constant<double>();
    return R * t * (z - 1 - (1 - beta) * q(z, a, b));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static double residual_entropy(double z, double a, double b,
                                 double beta) noexcept {
    constexpr auto R = gas_constant<double>();
    return R * (std::log(z - b) + beta * q(z, a, b));
  }

  // Constructors

  peng_robinson_eos() = default;

  /// @brief Constructs Peng-Robinson EoS
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  peng_robinson_eos(double pc, double tc, double omega)
      : base_type{pc, tc}, omega_{omega}, m_{m(omega)} {}

  peng_robinson_eos(const peng_robinson_eos&) = default;
  peng_robinson_eos(peng_robinson_eos&&) = default;

  peng_robinson_eos& operator=(const peng_robinson_eos&) = default;
  peng_robinson_eos& operator=(peng_robinson_eos&&) = default;

  // Member functions

  /// @brief Set parameters
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  void set_params(double pc, double tc, double omega) noexcept {
    this->base_type::set_params(pc, tc);
    omega_ = omega;
    m_ = m(omega);
  }

  // Member functions

  /// @brief Computes the correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  double alpha(double tr) const noexcept {
    const auto a = 1 + m_ * (1 - std::sqrt(tr));
    return a * a;
  }

  /// @brief Computes \f$ \beta = \frac{\mathrm{d} \ln \alpha}{\mathrm{d} \ln
  /// T_r } \f$
  /// @param[in] tr Reduced temperature
  double beta(double tr) const noexcept {
    const auto sqrt_tr = std::sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -m_ * sqrt_tr / a;
  }

 private:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static double m(double omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  /// @brief Computes a coefficient appearing in the calculation of fugacity
  /// coefficient, residual enthalpy and entropy.
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double q(double z, double a, double b) noexcept {
#ifndef M_SQRT2
#error M_SQRT2 is not defined!
#endif
    constexpr double sqrt2 = M_SQRT2;
    constexpr auto delta1 = 1 + sqrt2;
    constexpr auto delta2 = 1 - sqrt2;
    return a / (2 * sqrt2 * b) * std::log((z + delta1 * b) / (z + delta2 * b));
  }

  /*
  /// @brief Computes  \f[ \gamma = \frac{T_r^2}{\alpha} \cdot
  /// \frac{\mathrm{d}^2 \alpha}{\mathrm{d} T_r^2} \f]
  /// @param[in] tr Reduced temperature
  double  gamma(double tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    const auto alpha = a * a;
    return m_ / (2 * alpha) * (m_ * tr + a * sqrt_tr);
  }
  */

  /// Acentric factor
  double omega_;
  /// \f$ m = 0.3796 + 1.485 \omega - 0.1644 \omega^2 + 0.01667 \omega^3 \f$
  double m_;
};

/// @brief Makes Peng-Robinson EoS
/// @tparam T Scalar type
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
inline peng_robinson_eos make_peng_robinson_eos(double pc, double tc,
                                                double omega) {
  return {pc, tc, omega};
}

}  // namespace eos
