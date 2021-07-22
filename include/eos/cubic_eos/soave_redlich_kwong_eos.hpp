#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos/cubic_eos_base.hpp"  // eos::CubicEosBase
#include "eos/math/cubic_equation.hpp"       // eos::cubic_equation

namespace eos {

class SoaveRedlichKwongEos;

template <>
struct CubicEosTraits<SoaveRedlichKwongEos> {
  static constexpr double omegaA = 0.42748;
  static constexpr double omegaB = 0.08664;
};

/// @brief Soave-Redlich-Kwong EoS.
class SoaveRedlichKwongEos
    : public CubicEosBase<SoaveRedlichKwongEos, true> {
 public:
  using Base = CubicEosBase<SoaveRedlichKwongEos, true>;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static double pressure(double t, double v, double a, double b) noexcept {
    constexpr auto R = gasConstant<double>();
    return R * t / (v - b) - a / (v * (v + b));
  }

  /// @brief Computes the coefficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static cubic_equation zfactorCubicEq(double a, double b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static double lnFugacityCoeff(double z, double a, double b) noexcept {
    return z - 1 - std::log(z - b) - a / b * std::log((z + b) / z);
  }

  /// @brief Computes a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static double fugacityCoeff(double z, double a, double b) noexcept {
    return std::exp(lnFugacityCoeff(z, a, b));
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static double residualEnthalpy(double z, double t, double a, double b,
                                  double beta) noexcept {
    constexpr auto R = gasConstant<double>();
    return R * t * (z - 1 - a / b * (1 - beta) * std::log((z + b) / z));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double residualEntropy(double z, double a, double b,
                                 double beta) noexcept {
    constexpr auto R = gasConstant<double>();
    return R * (std::log(z - b) + a / b * beta * std::log((z + b) / z));
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double residualHelmholtzEnergy(double z, double t, double a,
                                          double b) noexcept {
    constexpr auto R = gasConstant<double>();
    return R * t * (std::log(z - b) + a / b * std::log((z + b) / z));
  }

  // Constructors

  SoaveRedlichKwongEos() = default;

  /// @brief Constructs Soave-Redlich-Kwong EoS
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  SoaveRedlichKwongEos(double pc, double tc, double omega)
      : Base{pc, tc}, omega_{omega}, m_{m(omega)} {}

  SoaveRedlichKwongEos(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos(SoaveRedlichKwongEos&&) = default;

  SoaveRedlichKwongEos& operator=(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos& operator=(SoaveRedlichKwongEos&&) = default;

  // Member functions

  /// @brief Set parameters
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  void setParams(double pc, double tc, double omega) noexcept {
    this->Base::setParams(pc, tc);
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
inline SoaveRedlichKwongEos makeSoaveRedlichKwongEos(double pc,
                                                            double tc,
                                                            double omega) {
  return {pc, tc, omega};
}

}  // namespace eos