#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/CubicEosBase.hpp"           // eos::CubicEosBase
#include "eos/MathematicalConstants.hpp"  // eos::sqrtTwo
#include "eos/SoaveCorrectionFactor.hpp"  // eos::SoaveCorrectionFactor

namespace eos {

/**
 * @brief EoS policy for PengRobinsonEos
 *
 * @tparam Scalar
 */
template <typename Scalar>
struct PengRobinsonEosPolicy {
  /// Coefficient for attraction parameter
  static constexpr Scalar omegaA = 0.45724;
  /// Coefficient for repulsion parameter
  static constexpr Scalar omegaB = 0.07780;

  /**
   * @brief Compute pressure at given temperature and volume
   *
   * @param[in] t temperature
   * @param[in] v volume
   * @param[in] a attraction parameter
   * @param[in] b repulsion parameter
   * @returns Scalar pressure
   */
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t / (v - b) - a / (v * (v + b) + b * (v - b));
  }

  /**
   * @brief Compute coefficients of the cubic equation of z-factor
   *
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @returns std::array<Scalar, 3> coefficients of the cubic equation of
   * z-factor.
   */
  static std::array<Scalar, 3> zfactorCubicEq(const Scalar& a,
                                              const Scalar& b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  /**
   * @brief Compute the natural log of a fugacity coefficient
   *
   * @param[in] z Z-factor
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @returns Scalar the natural log of a fugacity coefficient
   */
  static Scalar lnFugacityCoeff(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    using std::log;
    return z - 1 - log(z - b) - calcQ(z, a, b);
  }

  /**
   * @brief Compute residual enthalpy
   *
   * @param[in] z Z-factor
   * @param[in] t temperature
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @param[in] beta the derivative of correction factor
   * @returns Scalar residual enthalpy
   */
  static Scalar residualEnthalpy(const Scalar& z, const Scalar& t,
                                 const Scalar& a, const Scalar& b,
                                 const Scalar& beta) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (z - 1 - (1 - beta) * calcQ(z, a, b));
  }

  /**
   * @brief Compute residual entropy
   *
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @param[in] beta the derivative of correction factor
   * @returns Scalar residual entropy
   */
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * (log(z - b) + beta * calcQ(z, a, b));
  }

  /**
   * @brief Compute redisual Helmholtz energy
   *
   * @param[in] z Z-factor
   * @param[in] t Temperature
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns Scalar residual Helmholtz energy
   */
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (log(z - b) + calcQ(z, a, b));
  }

 private:
  /**
   * @brief Compute the coefficient Q
   *
   * Q is defined by
   * \f[
   *    Q := \frac{A}{2\sqrt{2}B} \log \frac{Z + (1 + \sqrt{2})B}{Z +
   *         (1 - \sqrt{2})B}
   * \f]
   *
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns Scalar the coefficient Q
   */
  static Scalar calcQ(const Scalar& z, const Scalar& a,
                      const Scalar& b) noexcept {
    using std::log;
    constexpr auto sqrt2 = sqrtTwo<Scalar>();
    constexpr auto delta1 = 1 + sqrt2;
    constexpr auto delta2 = 1 - sqrt2;
    return a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z + delta2 * b));
  }
};

/**
 * @brief Peng-Robinson EoS
 *
 * EoS proposed by Peng and Robinson (1976).
 *
 * @tparam Scalar scalar
 */
template <typename Scalar>
class PengRobinsonEos
    : public CubicEosBase<Scalar, PengRobinsonEosPolicy<Scalar>,
                          SoaveCorrectionFactor<Scalar>> {
 public:
  /// Base class
  using Base = CubicEosBase<Scalar, PengRobinsonEosPolicy<Scalar>,
                            SoaveCorrectionFactor<Scalar>>;

  PengRobinsonEos() = default;

  /**
   * @brief Default calculator of coefficient m for SoaveCorrectionFactor
   *
   * \f[
   *    m = 0.37464 + 1.54226 \omega - 0.2699 \omega^2
   * \f]
   * where \f$\omega\f$ is acentric factor.
   *
   * @tparam Scalar scalar
   */
  struct DefaultCalculator {
    /**
     * @brief Compute coefficient m for SoaveCorrectionFactor
     *
     * @param[in] omega acentric factor
     * @return Scalar coefficient m
     */
    Scalar operator()(const Scalar& omega) const noexcept {
      return 0.37464 + omega * (1.54226 - omega * 0.2699);
    }
  };

  /**
   * @brief Construct a new Peng Robinson Eos object
   *
   * @tparam F function to calculate m in SoaveCorrectionFactor
   * @param[in] pc critical pressure
   * @param[in] tc critical temperature
   * @param[in] omega acentric factor
   * @param[in] f functor to calculate m in SoaveCorrectionFactot (optional)
   */
  template <typename F = DefaultCalculator>
  PengRobinsonEos(const Scalar& pc, const Scalar& tc, const Scalar& omega,
                  F&& f = F{})
      : Base{pc, tc, SoaveCorrectionFactor{f(omega)}}, omega_{omega} {}

  PengRobinsonEos(const PengRobinsonEos&) = default;
  PengRobinsonEos(PengRobinsonEos&&) = default;

  PengRobinsonEos& operator=(const PengRobinsonEos&) = default;
  PengRobinsonEos& operator=(PengRobinsonEos&&) = default;

  /**
   * @brief Set acentric factor
   *
   * @tparam F function to calculate m in SoaveCorrectionFactor
   * @param[in] omega acentric factor
   * @param[in] f functor to calculate m in SoaveCorrectionFactor (optional)
   */
  template <typename F = DefaultCalculator>
  void setAcentricFactor(const Scalar& omega, F&& f = F{}) {
    omega_ = omega;
    this->setCorrectionFactor(SoaveCorrectionFactor{f(omega)});
  }

 private:
  Scalar omega_;  ///< Acentric factor
};

/**
 * @brief Create a new PengRobinsonEos object
 *
 * @tparam Scalar scalar
 * @tparam F function to calculate m in SoaveCorrectionFactor
 * @param[in] pc critical pressure
 * @param[in] tc critical temperature
 * @param[in] omega acentric factor
 * @returns PengRobinsonEos<Scalar> PengRobinsonEos object
 */
template <typename Scalar,
          typename F = PengRobinsonEos<Scalar>::DefaultCalculator>
inline PengRobinsonEos<Scalar> makePengRobinsonEos(const Scalar& pc,
                                                   const Scalar& tc,
                                                   const Scalar& omega,
                                                   F&& func = F{}) {
  return {pc, tc, omega, std::move(func)};
}

}  // namespace eos
