#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/cubic_eos_base.hpp"           // eos::CubicEosBase
#include "eos/mathematical_constants.hpp"   // eos::sqrtTwo
#include "eos/soave_correction_policy.hpp"  // eos::SoaveCorrectionPolicy

namespace eos {

/// @brief Peng-Robinson EoS.
template <typename Scalar>
struct PengRobinsonEosPolicy {
  static constexpr Scalar omegaA = 0.45724;
  static constexpr Scalar omegaB = 0.07780;

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
    return R * t / (v - b) - a / (v * (v + b) + b * (v - b));
  }

  /// @brief Computes coefficients of the cubic equation of z-factor.
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor.
  static std::array<Scalar, 3> zfactorCubicEq(const Scalar& a,
                                              const Scalar& b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static Scalar lnFugacityCoeff(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    using std::log;
    return z - 1 - log(z - b) - calcQ(z, a, b);
  }

  /// @brief Computes a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static Scalar fugacityCoeff(const Scalar& z, const Scalar& a,
                              const Scalar& b) noexcept {
    using std::exp;
    return exp(lnFugacityCoeff(z, a, b));
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
    return R * t * (z - 1 - (1 - beta) * calcQ(z, a, b));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * (log(z - b) + beta * calcQ(z, a, b));
  }

  /// @brief Computes redisual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (log(z - b) + calcQ(z, a, b));
  }

 private:
  /// @brief Computes a coefficient appearing in the calculation of fugacity
  /// coefficient, residual enthalpy and entropy.
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar calcQ(const Scalar& z, const Scalar& a,
                      const Scalar& b) noexcept {
    using std::log;
    constexpr auto sqrt2 = sqrtTwo<Scalar>();
    constexpr auto delta1 = 1 + sqrt2;
    constexpr auto delta2 = 1 - sqrt2;
    return a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z + delta2 * b));
  }
};

template <typename Scalar>
class PengRobinsonEos
    : public CubicEosBase<Scalar, PengRobinsonEosPolicy<Scalar>,
                          SoaveCorrectionPolicy<Scalar>> {
 public:
  using Base = CubicEosBase<Scalar, PengRobinsonEosPolicy<Scalar>,
                            SoaveCorrectionPolicy<Scalar>>;

  PengRobinsonEos() = default;

  PengRobinsonEos(const Scalar& pc, const Scalar& tc, const Scalar& omega)
      : Base{pc, tc, SoaveCorrectionPolicy{calcM(omega)}}, omega_{omega} {}

  PengRobinsonEos(const PengRobinsonEos&) = default;
  PengRobinsonEos(PengRobinsonEos&&) = default;

  PengRobinsonEos& operator=(const PengRobinsonEos&) = default;
  PengRobinsonEos& operator=(PengRobinsonEos&&) = default;

  void setAcentricFactor(const Scalar& omega) {
    omega_ = omega;
    this->correctionPolicy().m() = calcM(omega);
  }

 private:
  static Scalar calcM(const Scalar& omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  Scalar omega_;
};

/// @brief Makes Peng-Robinson EoS
/// @tparam T Scalar type
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename Scalar>
inline PengRobinsonEos<Scalar> makePengRobinsonEos(const Scalar& pc,
                                                   const Scalar& tc,
                                                   const Scalar& omega) {
  return {pc, tc, omega};
}

}  // namespace eos
