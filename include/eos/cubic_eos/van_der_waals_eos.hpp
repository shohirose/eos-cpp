#pragma once

#include <array>  // std::array
#include <cmath>  // std::exp, std::log

#include "eos/cubic_eos/cubic_eos_base.hpp"  // eos::CubicEosBase

namespace eos {

class VanDerWaalsEos;

template <>
struct CubicEosTraits<VanDerWaalsEos> {
  static constexpr double omegaA = 0.421875;
  static constexpr double omegaB = 0.125;
};

/// @brief Van der Waals Equations of State
class VanDerWaalsEos : public CubicEosBase<VanDerWaalsEos, false> {
 public:
  using Base = CubicEosBase<VanDerWaalsEos, false>;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static double pressure(double t, double v, double a, double b) noexcept {
    return gasConstant<double>() * t / (v - b) - a / (v * v);
  }

  /// @brief Computes coefficients of the cubic equation of Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<double, 3> zfactorCubicEq(double a, double b) noexcept {
    return {-b - 1, a, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static double lnFugacityCoeff(double z, double a, double b) noexcept {
    return -std::log(z - b) - a / z + z - 1;
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
  static double residualEnthalpy(double z, double t, double a,
                                  double b) noexcept {
    return gasConstant<double>() * t * (z - 1 - a / z);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double residualEntropy(double z, double a, double b) noexcept {
    return gasConstant<double>() * (std::log(z - b));
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static double residualHelmholtzEnergy(double z, double t, double a,
                                          double b) noexcept {
    constexpr auto R = gasConstant<double>();
    return R * t * (std::log(z - b) + a / z);
  }

  VanDerWaalsEos() = default;

  VanDerWaalsEos(double pc, double tc) noexcept : Base{pc, tc} {}

  VanDerWaalsEos(const VanDerWaalsEos &) = default;
  VanDerWaalsEos(VanDerWaalsEos &&) = default;

  VanDerWaalsEos &operator=(const VanDerWaalsEos &) = default;
  VanDerWaalsEos &operator=(VanDerWaalsEos &&) = default;

  void setParams(double pc, double tc) noexcept {
    this->Base::setParams(pc, tc);
  }
};

/// @brief Makes van der Waals EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
inline VanDerWaalsEos makeVanDerWaalsEos(double pc, double tc) {
  return {pc, tc};
}

}  // namespace eos