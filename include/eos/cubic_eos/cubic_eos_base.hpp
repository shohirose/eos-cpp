#pragma once

#include <vector>  // std::vector

#include "eos/common/thermodynamic_constants.hpp"  // eos::gasConstant
#include "eos/cubic_eos/isobaric_isothermal_state.hpp"
#include "eos/cubic_eos/isothermal_line.hpp"

namespace eos {

template <typename Eos>
struct CubicEosTraits {};

template <typename Derived>
class CubicEosCrtpBase {
 public:
  static constexpr auto omegaA = CubicEosTraits<Derived>::omegaA;
  static constexpr auto omegaB = CubicEosTraits<Derived>::omegaB;

  CubicEosCrtpBase() = default;
  CubicEosCrtpBase(const CubicEosCrtpBase &) = default;
  CubicEosCrtpBase(CubicEosCrtpBase &&) = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosCrtpBase(double pc, double tc) noexcept
      : pc_{pc},
        tc_{tc},
        ac_{this->criticalAttractionParam(pc, tc)},
        bc_{this->criticalRepulsionParam(pc, tc)} {}

  CubicEosCrtpBase &operator=(const CubicEosCrtpBase &) = default;
  CubicEosCrtpBase &operator=(CubicEosCrtpBase &&) = default;

  // Static functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static double criticalAttractionParam(double pc, double tc) noexcept {
    constexpr auto R = gasConstant<double>();
    return (omegaA * R * R) * tc * tc / pc;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static double criticalRepulsionParam(double pc, double tc) noexcept {
    constexpr auto R = gasConstant<double>();
    return (omegaB * R) * tc / pc;
  }

  /// @brief Returns reduced attraction parameter at a given pressure and
  /// temperature without temperature correction.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static double reducedAttractionParam(double pr, double tr) noexcept {
    return omegaA * pr / (tr * tr);
  }

  /// @brief Returns reduced repulsion parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static double reducedRepulsionParam(double pr, double tr) noexcept {
    return omegaB * pr / tr;
  }

  // Member functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void setParams(double pc, double tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    ac_ = criticalAttractionParam(pc, tc);
    bc_ = criticalRepulsionParam(pc, tc);
  }

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  double reducedPressure(double p) const noexcept { return p / pc_; }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature
  double reducedTemperature(double t) const noexcept { return t / tc_; }

 protected:
  /// @brief Get reference to derived class object
  Derived &derived() noexcept { return static_cast<Derived &>(*this); }

  /// @brief Get const reference to derived class object
  const Derived &derived() const noexcept {
    return static_cast<const Derived &>(*this);
  }

  double pc_;  /// Critical pressure
  double tc_;  /// Critical temperature
  double ac_;  /// Critical attraction parameter
  double bc_;  /// Critical repulsion parameter
};

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam Derived Concrete EoS class
///
/// Derived EoS classes must have the following static functions:
///    - pressure(t, v, a, b)
///    - zfactorCubicEq(ar, br)
///    - fugacityCoeff(z, ar, br)
///    - residualEnthalpy(z, t, ar, br, beta)
///    - residualEntropy(z, ar, br, beta)
///    - alpha(tr)
///    - beta()
/// where t is temperature, v is volume, a is attraction parameter, b is
/// repulsion parameter, ar is reduced attraction parameter, br is reduced
/// repulsion parameter, z is z-factor, and beta is a temperature correction
/// factor.
///
/// CubicEosTraits class specialized for each concrete EoS class must be
/// defined in the detail namespace. CubicEosTraits class must define the
/// following types and constants:
///    - double: double type
///    - omegaA: Constant for attraction parameter
///    - omegaB: Constant for repulsion parameter
///
template <typename Derived, bool UseTemperatureCorrectionFactor>
class CubicEosBase : public CubicEosCrtpBase<Derived> {
 public:
  using Base = CubicEosCrtpBase<Derived>;
  static constexpr auto omegaA = Base::omegaA;
  static constexpr auto omegaB = Base::omegaB;

  // Constructors

  CubicEosBase() = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosBase(double pc, double tc) noexcept : Base{pc, tc} {}

  CubicEosBase(const CubicEosBase &) = default;
  CubicEosBase(CubicEosBase &&) = default;

  CubicEosBase &operator=(const CubicEosBase &) = default;
  CubicEosBase &operator=(CubicEosBase &&) = default;

  // Member functions

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  IsothermalLine<Derived> createIsothermalLine(double t) const noexcept {
    const auto tr = this->reducedTemperature(t);
    const auto alpha = this->derived().alpha(tr);
    return {t, alpha * ac_, bc_};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  IsobaricIsothermalState<Derived, UseTemperatureCorrectionFactor>
  createIsobaricIsothermalState(double p, double t) const noexcept {
    const auto pr = this->reducedPressure(p);
    const auto tr = this->reducedTemperature(t);
    const auto ar =
        this->derived().alpha(tr) * this->reducedAttractionParam(pr, tr);
    const auto br = this->reducedRepulsionParam(pr, tr);
    const auto beta = this->derived().beta(tr);
    return {t, ar, br, beta};
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  double pressure(double t, double v) const noexcept {
    const auto tr = this->reducedTemperature(t);
    const auto a = this->derived().alpha(tr) * this->attractionParam();
    const auto b = this->repulsionParam();
    return Derived::pressure(t, v, a, b);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  std::vector<double> zfactor(double p, double t) const noexcept {
    return this->createIsobaricIsothermalState(p, t).zfactor();
  }
};

template <typename Derived>
class CubicEosBase<Derived, false> : public CubicEosCrtpBase<Derived> {
 public:
  using Base = CubicEosCrtpBase<Derived>;
  static constexpr auto omegaA = Base::omegaA;
  static constexpr auto omegaB = Base::omegaB;

  // Constructors

  CubicEosBase() = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosBase(double pc, double tc) noexcept : Base{pc, tc} {}

  CubicEosBase(const CubicEosBase &) = default;
  CubicEosBase(CubicEosBase &&) = default;

  CubicEosBase &operator=(const CubicEosBase &) = default;
  CubicEosBase &operator=(CubicEosBase &&) = default;

  // Member functions

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  IsothermalLine<Derived> createIsothermalLine(double t) const noexcept {
    return {t, ac_, bc_};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  IsobaricIsothermalState<Derived, false> createIsobaricIsothermalState(
      double p, double t) const noexcept {
    const auto pr = this->reducedPressure(p);
    const auto tr = this->reducedTemperature(t);
    const auto ar = this->reducedAttractionParam(pr, tr);
    const auto br = this->reducedRepulsionParam(pr, tr);
    return {t, ar, br};
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  double pressure(double t, double v) const noexcept {
    return Derived::pressure(t, v, ac_, bc_);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  std::vector<double> zfactor(double p, double t) const noexcept {
    return this->createIsobaricIsothermalState(p, t).zfactor();
  }
};

}  // namespace eos