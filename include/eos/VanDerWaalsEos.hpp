#pragma once

#include <array>  // std::array
#include <cmath>  // std::log

#include "eos/CubicEosBase.hpp"              // eos::CubicEosBase
#include "eos/IdentityCorrectionFactor.hpp"  // eos::IdentityCorrectionFactor

namespace eos {

/**
 * @brief EoS policy for Van der Waals EoS
 *
 * @tparam Scalar scalar
 */
template <typename Scalar>
struct VanDerWaalsEosPolicy {
  /// Coefficient for attraction parameter
  static constexpr Scalar omegaA = 0.421875;
  /// Coefficient for repulsion parameter
  static constexpr Scalar omegaB = 0.125;

  /**
   * @brief Compute pressure at given temperature and volume.
   *
   * @param t temperature
   * @param v volume
   * @param a attraction parameter
   * @param b respulsion parameter
   * @param alpha temperature correction factor
   * @return Scalar pressure
   */
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b, const Scalar& alpha) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t / (v - b) - alpha * a / (v * v);
  }

  /**
   * @brief Compute coefficients of the cubic equation of Z-factor
   *
   * @param a reduced attraction parameter
   * @param b reduced repulsion parameter
   * @return std::array<Scalar, 3> coefficients of the cubic equation of
   * Z-factor
   */
  static std::array<Scalar, 3> zfactorCubicEq(const Scalar& a,
                                              const Scalar& b) noexcept {
    return {-b - 1, a, -a * b};
  }

  /**
   * @brief Compute the natural logarithm of a fugacity coefficient
   * @param[in] z Z-factor
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @return The natural log of a fugacity coefficient
   */
  static Scalar lnFugacityCoeff(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    using std::log;
    return z - 1 - log(z - b) - a / z;
  }

  /**
   * @brief Compute residual enthalpy
   * @param[in] z Z-factor
   * @param[in] t temperature
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @param[in] beta derivative of correction factor
   * @return Scalar residual enthalpy
   */
  static Scalar residualEnthalpy(const Scalar& z, const Scalar& t,
                                 const Scalar& a, const Scalar& b,
                                 const Scalar& beta) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (z - 1 - a / z * (1 - beta));
  }

  /**
   * @brief Compute residual entropy
   * @param[in] z Z-factor
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @param[in] beta derivative of correction factor
   * @return Scalar residual entropy
   */
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * (-log(z / (z - b)) + a / z * beta);
  }

  /**
   * @brief Compute residual Helmholtz energy
   * @param[in] z Z-factor
   * @param[in] t temperature
   * @param[in] a reduced attraction parameter
   * @param[in] b reduced repulsion parameter
   * @return Scalar residual Helmholtz energy
   */
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (log(z / (z - b)) - a / z);
  }
};

/**
 * @brief Van der Waals EoS
 *
 * Van der Waals EoS is expressed by
 * \f[
 *    P = \frac{RT}{V - b} - \frac{a}{V^2}
 * \f]
 * where \f$P\f$ is pressure, \f$T\f$ is temperature, \f$V\f$ is volume, \f$a\f$
 * is attraction parameter, \f$b\f$ is repulsion parameter, and \f$R\f$ is gas
 * constant.
 *
 * @tparam Scalar
 */
template <typename Scalar>
class VanDerWaalsEos : public CubicEosBase<Scalar, VanDerWaalsEosPolicy<Scalar>,
                                           IdentityCorrectionFactor<Scalar>> {
 public:
  /// Base class
  using Base = CubicEosBase<Scalar, VanDerWaalsEosPolicy<Scalar>,
                            IdentityCorrectionFactor<Scalar>>;

  VanDerWaalsEos() = default;

  /**
   * @brief Construct a new VanDerWaalsEos object
   *
   * @param pc critical pressure
   * @param tc critical temperature
   */
  VanDerWaalsEos(const Scalar& pc, const Scalar& tc)
      : Base{pc, tc, IdentityCorrectionFactor<Scalar>{}} {}

  VanDerWaalsEos(const VanDerWaalsEos&) = default;
  VanDerWaalsEos(VanDerWaalsEos&&) = default;

  VanDerWaalsEos& operator=(const VanDerWaalsEos&) = default;
  VanDerWaalsEos& operator=(VanDerWaalsEos&&) = default;
};

/**
 * @brief Create a new VanDerWaalsEos object
 *
 * @tparam Scalar scalar
 * @param[in] pc critical pressure
 * @param[in] tc critical temperature
 * @return VanDerWaalsEos<Scalar>
 */
template <typename Scalar>
inline VanDerWaalsEos<Scalar> makeVanDerWaalsEos(const Scalar& pc,
                                                 const Scalar& tc) {
  return {pc, tc};
}

}  // namespace eos