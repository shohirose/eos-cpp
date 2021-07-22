#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos_base.hpp"           // eos::CubicEosBase
#include "eos/soave_correction_policy.hpp"  // eos::SoaveCorrectionPolicy

namespace eos {

/// @brief Soave-Redlich-Kwong EoS.
template <typename Scalar>
struct SoaveRedlichKwongEosPolicy {
  static constexpr Scalar omegaA = 0.42748;
  static constexpr Scalar omegaB = 0.08664;

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
    using std::log;
    return z - 1 - log(z - b) - a / b * log((z + b) / z);
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
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (z - 1 - a / b * (1 - beta) * log((z + b) / z));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * (log(z - b) + a / b * beta * log((z + b) / z));
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (log(z - b) + a / b * log((z + b) / z));
  }
};

template <typename Scalar>
class SoaveRedlichKwongEos
    : public CubicEosBase<Scalar, SoaveRedlichKwongEosPolicy<Scalar>,
                          SoaveCorrectionPolicy<Scalar>> {
 public:
  using Base = CubicEosBase<Scalar, SoaveRedlichKwongEosPolicy<Scalar>,
                            SoaveCorrectionPolicy<Scalar>>;

  SoaveRedlichKwongEos() = default;

  SoaveRedlichKwongEos(const Scalar& pc, const Scalar& tc, const Scalar& omega)
      : Base{pc, tc, SoaveCorrectionPolicy{calcM(omega)}}, omega_{omega} {}

  SoaveRedlichKwongEos(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos(SoaveRedlichKwongEos&&) = default;

  SoaveRedlichKwongEos& operator=(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos& operator=(SoaveRedlichKwongEos&&) = default;

  void setAcentricFactor(const Scalar& omega) {
    omega_ = omega;
    this->correctionPolicy().m() = calcM(omega);
  }

 private:
  static Scalar calcM(const Scalar& omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  Scalar omega_;
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