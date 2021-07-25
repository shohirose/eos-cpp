#pragma once

#include <array>  // std::array
#include <cmath>  // std::exp, std::log

#include "eos/cubic_eos_base.hpp"              // eos::CubicEosBase
#include "eos/identity_correction_factor.hpp"  // eos::IdentityCorrectionFactor

namespace eos {

/// @brief Policy for Van der Waals EoS
template <typename Scalar>
struct VanDerWaalsEosPolicy {
  static constexpr Scalar omegaA = 0.421875;
  static constexpr Scalar omegaB = 0.125;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t / (v - b) - a / (v * v);
  }

  /// @brief Computes coefficients of the cubic equation of Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<Scalar, 3> zfactorCubicEq(const Scalar& a,
                                              const Scalar& b) noexcept {
    return {-b - 1, a, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static Scalar lnFugacityCoeff(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    using std::log;
    return z - 1 - log(z - b) - a / z;
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualEnthalpy(const Scalar& z, const Scalar& t,
                                 const Scalar& a, const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (z - 1 - a / z);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    using std::log;
    constexpr auto R = gasConstant<Scalar>();
    return R * (log(z - b));
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
    return R * t * (log(z - b) + a / z);
  }
};

template <typename Scalar>
class VanDerWaalsEos : public CubicEosBase<Scalar, VanDerWaalsEosPolicy<Scalar>,
                                           IdentityCorrectionFactor> {
 public:
  using Base = CubicEosBase<Scalar, VanDerWaalsEosPolicy<Scalar>,
                            IdentityCorrectionFactor>;

  VanDerWaalsEos() = default;

  VanDerWaalsEos(const Scalar& pc, const Scalar& tc)
      : Base{pc, tc, IdentityCorrectionFactor{}} {}

  VanDerWaalsEos(const VanDerWaalsEos&) = default;
  VanDerWaalsEos(VanDerWaalsEos&&) = default;

  VanDerWaalsEos& operator=(const VanDerWaalsEos&) = default;
  VanDerWaalsEos& operator=(VanDerWaalsEos&&) = default;
};

/// @brief Makes van der Waals EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
template <typename Scalar>
inline VanDerWaalsEos<Scalar> makeVanDerWaalsEos(const Scalar& pc,
                                                 const Scalar& tc) {
  return {pc, tc};
}

}  // namespace eos