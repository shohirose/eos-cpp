#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos/cubic_eos_base.hpp"  // eos::CubicEosBase

namespace eos {

template <typename Scalar>
class SoaveRedlichKwongEos;

template <typename Scalar_>
struct CubicEosTraits<SoaveRedlichKwongEos<Scalar_>> {
  using Scalar = Scalar_;
  static constexpr Scalar omegaA = 0.42748;
  static constexpr Scalar omegaB = 0.08664;
};

/// @brief Soave-Redlich-Kwong EoS.
template <typename Scalar>
class SoaveRedlichKwongEos
    : public CubicEosBase<SoaveRedlichKwongEos<Scalar>, true> {
 public:
  using Base = CubicEosBase<SoaveRedlichKwongEos<Scalar>, true>;

  // Static functions

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t / (v - b) - a / (v * (v + b));
  }

  /// @brief Computes the coefficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static std::array<Scalar, 3> zfactorCubicEq(const Scalar& a,
                                              const Scalar& b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static Scalar lnFugacityCoeff(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    return z - 1 - std::log(z - b) - a / b * std::log((z + b) / z);
  }

  /// @brief Computes a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static Scalar fugacityCoeff(const Scalar& z, const Scalar& a,
                              const Scalar& b) noexcept {
    return std::exp(lnFugacityCoeff(z, a, b));
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static Scalar residualEnthalpy(const Scalar& z, const Scalar& t,
                                 const Scalar& a, const Scalar& b,
                                 const Scalar& beta) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (z - 1 - a / b * (1 - beta) * std::log((z + b) / z));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * (std::log(z - b) + a / b * beta * std::log((z + b) / z));
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (std::log(z - b) + a / b * std::log((z + b) / z));
  }

  // Constructors

  SoaveRedlichKwongEos() = default;

  /// @brief Constructs Soave-Redlich-Kwong EoS
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  SoaveRedlichKwongEos(const Scalar& pc, const Scalar& tc, const Scalar& omega)
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
  void setParams(const Scalar& pc, const Scalar& tc,
                 const Scalar& omega) noexcept {
    this->Base::setParams(pc, tc);
    omega_ = omega;
    m_ = m(omega);
  }

  /// @brief Computes the correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  Scalar alpha(const Scalar& tr) const noexcept {
    const auto a = 1 + m_ * (1 - std::sqrt(tr));
    return a * a;
  }

  /// @brief Computes \f$ \beta = \frac{\mathrm{d} \ln \alpha}{\mathrm{d} \ln
  /// const Scalar&} \f$
  /// @param[in] tr Reduced temperature
  Scalar beta(const Scalar& tr) const noexcept {
    const auto sqrt_tr = std::sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -m_ * sqrt_tr / a;
  }

 private:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static Scalar m(const Scalar& omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  /// Acentric factor
  Scalar omega_;
  /// \f$ m = 0.3796 + 1.485 \omega - 0.1644 \omega^2 + 0.01667 \omega^3 \f$
  Scalar m_;
};

/// @brief Makes Soave-Redlich-Kwong EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename Scalar>
inline SoaveRedlichKwongEos<Scalar> makeSoaveRedlichKwongEos(
    const Scalar& pc, const Scalar& tc, const Scalar& omega) {
  return {pc, tc, omega};
}

}  // namespace eos