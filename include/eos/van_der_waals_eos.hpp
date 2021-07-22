#pragma once

#include <array>  // std::array
#include <cmath>  // std::exp, std::log

#include "eos/cubic_eos_base.hpp"  // eos::CubicEosBase

namespace eos {

template <typename Scalar>
class VanDerWaalsEos;

template <typename Scalar_>
struct CubicEosTraits<VanDerWaalsEos<Scalar_>> {
  using Scalar = Scalar_;
  static constexpr Scalar omegaA = 0.421875;
  static constexpr Scalar omegaB = 0.125;
};

/// @brief Van der Waals Equations of State
template <typename Scalar>
class VanDerWaalsEos : public CubicEosBase<VanDerWaalsEos<Scalar>, false> {
 public:
  using Base = CubicEosBase<VanDerWaalsEos, false>;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static Scalar pressure(const Scalar& t, const Scalar& v, const Scalar& a,
                         const Scalar& b) noexcept {
    return gasConstant<Scalar>() * t / (v - b) - a / (v * v);
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
    return -std::log(z - b) - a / z + z - 1;
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
  static Scalar residualEnthalpy(const Scalar& z, const Scalar& t,
                                 const Scalar& a, const Scalar& b) noexcept {
    return gasConstant<Scalar>() * t * (z - 1 - a / z);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b) noexcept {
    return gasConstant<Scalar>() * (std::log(z - b));
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
    return R * t * (std::log(z - b) + a / z);
  }

  VanDerWaalsEos() = default;

  VanDerWaalsEos(const Scalar& pc, const Scalar& tc) noexcept : Base{pc, tc} {}

  VanDerWaalsEos(const VanDerWaalsEos&) = default;
  VanDerWaalsEos(VanDerWaalsEos&&) = default;

  VanDerWaalsEos& operator=(const VanDerWaalsEos&) = default;
  VanDerWaalsEos& operator=(VanDerWaalsEos&&) = default;

  void setParams(const Scalar& pc, const Scalar& tc) noexcept {
    this->Base::setParams(pc, tc);
  }
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