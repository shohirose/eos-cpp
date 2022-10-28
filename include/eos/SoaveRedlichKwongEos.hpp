#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/CubicEosBase.hpp"           // eos::CubicEosBase
#include "eos/SoaveCorrectionFactor.hpp"  // eos::SoaveCorrectionFactor

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

  /**
   * @brief Compute pressure at given temperature and volume
   *
   * @param t temperature
   * @param v volume
   * @param a attraction parameter
   * @param b repulsion parameter
   * @param[in] alpha temperature correction factor
   * @return Scalar pressure
   */
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b, const Scalar& alpha) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t / (v - b) - alpha * a / (v * (v + b));
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
    return z - 1 - log(z - b) + a / b * log(z / (z + b));
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
    return R * t * (z - 1 + a / b * (1 - beta) * log(z / (z + b)));
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
    return R * (-log(z / (z - b)) - a / b * beta * log(z / (z + b)));
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
    return R * t * (log(z / (z - b)) + a / b * log(z / (z + b)));
  }
};

/**
 * @brief Soave-Redlich-Kwong EoS
 *
 * Soave-Redlich-Kwong EoS proposed by Soave (1972) is defined by
 * \f[
 *    P = \frac{RT}{V - b} - \frac{\alpha(T) a}{V(V + b)}
 * \f]
 * where \f$P\f$ is pressure, \f$T\f$ is temperature, \f$V\f$ is volume, \f$a\f$
 * is attraction parameter, \f$b\f$ is repulsion parameter, and \f$\alpha\f$ is
 * correction factor for attraction parameter.
 *
 * The correction parameter for SRK EoS is defined by
 * \f[
 *    \alpha(T) := \left[ 1 + m \left( 1 - \sqrt{T_r} \right) \right]^2
 * \f]
 * where \f$m\f$ is a coefficient, and \f$T_r\f$ is reduced temperature.
 *
 * @tparam Scalar scalar
 */
template <typename Scalar>
class SoaveRedlichKwongEos
    : public CubicEosBase<Scalar, SoaveRedlichKwongEosPolicy<Scalar>,
                          SoaveCorrectionFactor<Scalar>> {
 public:
  /// Base class
  using Base = CubicEosBase<Scalar, SoaveRedlichKwongEosPolicy<Scalar>,
                            SoaveCorrectionFactor<Scalar>>;

  SoaveRedlichKwongEos() = default;

  /**
   * @brief Default calculator of coefficient m for SoaveCorrectionFactor for
   * SoaveRedlichKwongEos
   *
   * \f[
   *    m = 0.48 + 1.574 \omega - 0.176 \omega^2
   * \f]
   * where \f$\omega\f$ is acentric factor.
   *
   * @tparam Scalar scalar
   */
  struct DefaultCalculator {
    /**
     * @brief Compute coeffcient m for SoaveCorrecionFactor
     *
     * @param[in] omega acentric factor
     * @return Scalar coefficient m
     */
    Scalar operator()(const Scalar& omega) const noexcept {
      return 0.48 + (1.574 - 0.176 * omega) * omega;
    }
  };

  /**
   * @brief Construct a new SoaveRedlichKwongEos object
   *
   * @tparam F function to calculate m in SoaveCorrectionFactor
   * @param[in] pc critical pressure
   * @param[in] tc critical temperature
   * @param[in] omega acentric factor
   * @param[in] f function to calculate m in SoaveCorrectionFactor (optional)
   */
  template <typename F = DefaultCalculator>
  SoaveRedlichKwongEos(const Scalar& pc, const Scalar& tc, const Scalar& omega,
                       F&& f = F{})
      : Base{pc, tc, SoaveCorrectionFactor{f(omega)}}, omega_{omega} {}

  SoaveRedlichKwongEos(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos(SoaveRedlichKwongEos&&) = default;

  SoaveRedlichKwongEos& operator=(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos& operator=(SoaveRedlichKwongEos&&) = default;

  /**
   * @brief Set acentric factor
   *
   * @tparam F function to calculate m in SoaveCorrectionFactor
   * @param[in] omega acentric factor
   * @param[in] f function to calculate m in SoaveCorrectionFactor (optional)
   */
  template <typename F = DefaultCalculator>
  void setAcentricFactor(const Scalar& omega, F&& f = F{}) {
    omega_ = omega;
    this->correctionFactor().m() = f(omega);
  }

 private:
  Scalar omega_;  ///< acentric factor
};

/**
 * @brief Make a new SoaveRedlichKwongEos object
 *
 * @tparam Scalar scalar
 * @tparam F function to calculate m in SoaveCorrectionFactor
 * @param[in] pc critical pressure
 * @param[in] tc critical temperature
 * @param[in] omega acentric factor
 * @param[in] f function to calculate m in SoaveCorrectionFactor (optional)
 * @return SoaveRedlichKwongEos<Scalar>
 */
template <typename Scalar,
          typename F = typename SoaveRedlichKwongEos<Scalar>::DefaultCalculator>
inline SoaveRedlichKwongEos<Scalar> makeSoaveRedlichKwongEos(
    const Scalar& pc, const Scalar& tc, const Scalar& omega, F&& f = F{}) {
  return {pc, tc, omega, std::forward<F>(f)};
}

}  // namespace eos