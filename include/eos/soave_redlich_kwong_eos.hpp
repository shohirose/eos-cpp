#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos_base.hpp"           // eos::CubicEosBase
#include "eos/soave_correction_policy.hpp"  // eos::SoaveCorrectionPolicy

namespace eos {

/**
 * @brief EoS policy for Soave-Redlich-Kwong EoS.
 *
 * @tparam Scalar scalar
 */
template <typename Scalar>
struct SoaveRedlichKwongEosPolicy {
  /// Coefficient for attraction parameter
  static constexpr Scalar omegaA = 0.42748;
  /// Coefficient for repulsion parameter
  static constexpr Scalar omegaB = 0.08664;

  /// @name Public static functions
  //@{

  /**
   * @brief Compute pressure at given temperature and volume
   *
   * @param t temperature
   * @param v volume
   * @param a attraction parameter
   * @param b repulsion parameter
   * @return Scalar pressure
   */
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t / (v - b) - a / (v * (v + b));
  }

  /**
   * @brief Compute coefficients of the cubic equation of Z-factor.
   *
   * @param a reduced attraction parameter
   * @param b reduced repulsion parameter
   * @return std::array<Scalar, 3> coefficients of the cubic equation of
   * Z-factor
   */
  static std::array<Scalar, 3> zfactorCubicEq(const Scalar& a,
                                              const Scalar& b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  /**
   * @brief Computes natural log of fugacity coefficients
   *
   * @param z Z-factor
   * @param a reduced attraction parameter
   * @param b reduced repulsion parameter
   * @return Scalar natural log of fugacity coefficients
   */
  static Scalar lnFugacityCoeff(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    using std::log;
    return z - 1 - log(z - b) - a / b * log((z + b) / z);
  }

  /**
   * @brief Compute residual enthalpy
   *
   * @param z Z-factor
   * @param t temperature
   * @param a reduced attraction parameter
   * @param b reduced repulsion parameter
   * @param beta derivative of correction factor
   * @return Scalar residual enthalpy
   */
  static Scalar residualEnthalpy(const Scalar& z, const Scalar& t,
                                 const Scalar& a, const Scalar& b,
                                 const Scalar& beta) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (z - 1 - a / b * (1 - beta) * log((z + b) / z));
  }

  /**
   * @brief Compute residual entropy
   *
   * @param z Z-factor
   * @param a reduced attraction parameter
   * @param b reduced repulsion parameter
   * @param beta derivative of correction factor
   * @return Scalar residual entropy
   */
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * (log(z - b) + a / b * beta * log((z + b) / z));
  }

  /**
   * @brief Compute residual Helmholtz energy
   *
   * @param z Z-factor
   * @param t temperature
   * @param a reduced attraction parameter
   * @param b reduced repulsion parameter
   * @return Scalar residual Helmholtz energy
   */
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (log(z - b) + a / b * log((z + b) / z));
  }
  //@}
};

/**
 * @brief Soave-Redlich-Kwong EoS
 *
 * EoS proposed by Soave (1972).
 *
 * @tparam Scalar scalar
 */
template <typename Scalar>
class SoaveRedlichKwongEos
    : public CubicEosBase<Scalar, SoaveRedlichKwongEosPolicy<Scalar>,
                          SoaveCorrectionPolicy<Scalar>> {
 public:
  /// Base class
  using Base = CubicEosBase<Scalar, SoaveRedlichKwongEosPolicy<Scalar>,
                            SoaveCorrectionPolicy<Scalar>>;
  /// @name Constructors
  //@{
  SoaveRedlichKwongEos() = default;

  /**
   * @brief Construct a new Soave Redlich Kwong Eos object
   *
   * @param pc critical pressure
   * @param tc critical temperature
   * @param omega acentric factor
   */
  SoaveRedlichKwongEos(const Scalar& pc, const Scalar& tc, const Scalar& omega)
      : Base{pc, tc, SoaveCorrectionPolicy{calcM(omega)}}, omega_{omega} {}

  SoaveRedlichKwongEos(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos(SoaveRedlichKwongEos&&) = default;
  //@}

  SoaveRedlichKwongEos& operator=(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos& operator=(SoaveRedlichKwongEos&&) = default;

  /**
   * @brief Set acentric factor
   *
   * @param omega acentric factor
   */
  void setAcentricFactor(const Scalar& omega) {
    omega_ = omega;
    this->correctionFactor().m() = calcM(omega);
  }

 private:
  /**
   * @brief Compute coefficient m for SoaveCorrectionPolicy
   *
   * @param omega acentric factor
   * @return Scalar coefficient m
   */
  static Scalar calcM(const Scalar& omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  Scalar omega_;  ///< acentric factor
};

/**
 * @brief Make a new Soave-Redlich-Kwong EoS object
 *
 * @tparam Scalar scalar
 * @param pc critical pressure
 * @param tc critical temperature
 * @param omega acentric factor
 * @return SoaveRedlichKwongEos<Scalar>
 */
template <typename Scalar>
inline SoaveRedlichKwongEos<Scalar> makeSoaveRedlichKwongEos(
    const Scalar& pc, const Scalar& tc, const Scalar& omega) {
  return {pc, tc, omega};
}

}  // namespace eos