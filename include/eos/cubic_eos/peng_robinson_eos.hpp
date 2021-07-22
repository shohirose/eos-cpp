#pragma once

#include <array>  // std::array
#include <cmath>  // std::sqrt, std::exp, std::log

#include "eos/common/mathematical_constants.hpp"  // eos::sqrtTwo
#include "eos/cubic_eos/cubic_eos_base.hpp"       // eos::CubicEosBase

namespace eos {

template <typename Scalar>
class PengRobinsonEos;

template <typename Scalar_>
struct CubicEosTraits<PengRobinsonEos<Scalar_>> {
  using Scalar = Scalar_;
  static constexpr Scalar omegaA = 0.45724;
  static constexpr Scalar omegaB = 0.07780;
};

/// @brief Peng-Robinson EoS.
template <typename Scalar>
class PengRobinsonEos : public CubicEosBase<PengRobinsonEos<Scalar>, true> {
 public:
  using Base = CubicEosBase<PengRobinsonEos, true>;

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
    return z - 1 - std::log(z - b) - q(z, a, b);
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
    return R * t * (z - 1 - (1 - beta) * q(z, a, b));
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @param[in] beta Temperature correction factor
  static Scalar residualEntropy(const Scalar& z, const Scalar& a,
                                const Scalar& b, const Scalar& beta) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * (std::log(z - b) + beta * q(z, a, b));
  }

  /// @brief Computes redisual Helmholtz energy
  /// @param[in] z Z-factor
  /// @param[in] t Temperature
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar residualHelmholtzEnergy(const Scalar& z, const Scalar& t,
                                        const Scalar& a,
                                        const Scalar& b) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return R * t * (std::log(z - b) + q(z, a, b));
  }

  // Constructors

  PengRobinsonEos() = default;

  /// @brief Constructs Peng-Robinson EoS
  /// @param[in] pc Critical pressrue
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  PengRobinsonEos(const Scalar& pc, const Scalar& tc, const Scalar& omega)
      : Base{pc, tc}, omega_{omega}, m_{m(omega)} {}

  PengRobinsonEos(const PengRobinsonEos&) = default;
  PengRobinsonEos(PengRobinsonEos&&) = default;

  PengRobinsonEos& operator=(const PengRobinsonEos&) = default;
  PengRobinsonEos& operator=(PengRobinsonEos&&) = default;

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

  // Member functions

  /// @brief Computes the correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  Scalar alpha(const Scalar& tr) const noexcept {
    const auto a = 1 + m_ * (1 - std::sqrt(tr));
    return a * a;
  }

  /// @brief Computes \f$ \beta = \frac{d \ln \alpha}{d \ln T } \f$
  /// @param[in] tr Reduced temperature
  Scalar beta(const Scalar& tr) const noexcept {
    const auto sqrt_tr = std::sqrt(tr);
    return -m_ * sqrt_tr * (1 + m_ * (1 - sqrt_tr));
  }

 private:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static Scalar m(const Scalar& omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  /// @brief Computes a coefficient appearing in the calculation of fugacity
  /// coefficient, residual enthalpy and entropy.
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  static Scalar q(const Scalar& z, const Scalar& a, const Scalar& b) noexcept {
    constexpr auto sqrt2 = sqrtTwo<Scalar>();
    constexpr auto delta1 = 1 + sqrt2;
    constexpr auto delta2 = 1 - sqrt2;
    return a / (2 * sqrt2 * b) * std::log((z + delta1 * b) / (z + delta2 * b));
  }

  /*
  /// @brief Computes  \f[ \gamma = \frac{T_r^2}{\alpha} \cdot
  /// \frac{\mathrm{d}^2 \alpha}{\mathrm{d} T_r^2} \f]
  /// @param[in] tr Reduced temperature
  double  gamma(double tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    const auto alpha = a * a;
    return m_ / (2 * alpha) * (m_ * tr + a * sqrt_tr);
  }
  */

  /// Acentric factor
  Scalar omega_;
  /// \f$ m = 0.3796 + 1.485 \omega - 0.1644 \omega^2 + 0.01667 \omega^3 \f$
  Scalar m_;
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
