#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos/cubic_eos_base.hpp"  // eos::cubic_eos_base
#include "eos/math/cubic_equation.hpp"       // eos::cubic_equation

namespace eos {

class soave_redlich_kwong_eos;

template <>
struct cubic_eos_traits<soave_redlich_kwong_eos> {
  static constexpr double omega_a = 0.42748;
  static constexpr double omega_b = 0.08664;
};

/// @brief Soave-Redlich-Kwong EoS.
class soave_redlich_kwong_eos : public cubic_eos_base<soave_redlich_kwong_eos> {
 public:
  using base_type = cubic_eos_base<soave_redlich_kwong_eos>;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static double pressure(double t, double v, double a, double b) noexcept {
    constexpr auto R = gas_constant<double>();
    return R * t / (v - b) - a / (v * (v + b));
  }

  /// @brief Computes the coefficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static cubic_equation zfactor_cubic_eq(double a, double b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static double ln_fugacity_coeff(double z, double a, double b) noexcept {
    return z - 1 - std::log(z - b) - a / b * std::log((z + b) / z);
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
    return R * t * (z - 1 - a / b * (1 - beta) * std::log((z + b) / z));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double residual_entropy(double z, double a, double b,
                                 double beta) noexcept {
    constexpr auto R = gas_constant<double>();
    return R * (std::log(z - b) + a / b * beta * std::log((z + b) / z));
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double residual_helmholtz_energy(double z, double t, double a,
                                          double b) noexcept {
    constexpr auto R = gas_constant<double>();
    return R * t * (std::log(z - b) + a / b * std::log((z + b) / z));
  }

  /*
  /// @brief Computes residual molar specific heat at constant volume
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] gamma Temperature correction factor
  static double residual_specific_heat_at_const_volume(double z, const
  double &a, double b, double gamma) noexcept { using std::log;
    return gas_constant * gamma * a / b * log((z + b) / z);
  }
  */

  // Constructors

  soave_redlich_kwong_eos() = default;

  /// @brief Constructs Soave-Redlich-Kwong EoS
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  soave_redlich_kwong_eos(double pc, double tc, double omega)
      : base_type{pc, tc}, omega_{omega}, m_{m(omega)} {}

  soave_redlich_kwong_eos(const soave_redlich_kwong_eos&) = default;
  soave_redlich_kwong_eos(soave_redlich_kwong_eos&&) = default;

  soave_redlich_kwong_eos& operator=(const soave_redlich_kwong_eos&) = default;
  soave_redlich_kwong_eos& operator=(soave_redlich_kwong_eos&&) = default;

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

  /// @brief Computes the correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  double alpha(double tr) const noexcept {
    const auto a = 1 + m_ * (1 - std::sqrt(tr));
    return a * a;
  }

  /// @brief Computes \f$ \beta = \frac{\mathrm{d} \ln \alpha}{\mathrm{d} \ln
  /// double} \f$
  /// @param[in] tr Reduced temperature
  double beta(double tr) const noexcept {
    const auto sqrt_tr = std::sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -m_ * sqrt_tr / a;
  }

  /*
  /// @brief Computes  \f[ \gamma = \frac{T_r^2}{\alpha} \cdot
  /// \frac{\mathrm{d}^2 \alpha}{\mathrm{d} T_r^2} \f]
  /// @param[in] tr Reduced temperature
  double gamma(double tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    const auto alpha = a * a;
    return m_ / (2 * alpha) * (m_ * tr + a * sqrt_tr);
  }
  */

 private:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static double m(double omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  /// Acentric factor
  double omega_;
  /// \f$ m = 0.3796 + 1.485 \omega - 0.1644 \omega^2 + 0.01667 \omega^3 \f$
  double m_;
};

/// @brief Makes Soave-Redlich-Kwong EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
inline soave_redlich_kwong_eos make_soave_redlich_kwong_eos(double pc,
                                                            double tc,
                                                            double omega) {
  return {pc, tc, omega};
}

}  // namespace eos