#pragma once

#include <vector>  // std::vector

#include "eos/thermodynamic_constants.hpp"  // eos::gasConstant
#include "eos/isobaric_isothermal_state.hpp"
#include "eos/isothermal_line.hpp"

namespace eos {

template <typename Eos>
struct CubicEosTraits {};

template <typename Derived>
class CubicEosCrtpBase {
 public:
  using Scalar = typename CubicEosTraits<Derived>::Scalar;
  static constexpr auto omegaA = CubicEosTraits<Derived>::omegaA;
  static constexpr auto omegaB = CubicEosTraits<Derived>::omegaB;

  CubicEosCrtpBase() = default;
  CubicEosCrtpBase(const CubicEosCrtpBase&) = default;
  CubicEosCrtpBase(CubicEosCrtpBase&&) = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosCrtpBase(const Scalar& pc, const Scalar& tc) noexcept
      : pc_{pc},
        tc_{tc},
        ac_{this->criticalAttractionParam(pc, tc)},
        bc_{this->criticalRepulsionParam(pc, tc)} {}

  CubicEosCrtpBase& operator=(const CubicEosCrtpBase&) = default;
  CubicEosCrtpBase& operator=(CubicEosCrtpBase&&) = default;

  // Static functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static Scalar criticalAttractionParam(const Scalar& pc,
                                        const Scalar& tc) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return (omegaA * R * R) * tc * tc / pc;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static Scalar criticalRepulsionParam(const Scalar& pc,
                                       const Scalar& tc) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return (omegaB * R) * tc / pc;
  }

  /// @brief Returns reduced attraction parameter at a given pressure and
  /// temperature without temperature correction.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static Scalar reducedAttractionParam(const Scalar& pr,
                                       const Scalar& tr) noexcept {
    return omegaA * pr / (tr * tr);
  }

  /// @brief Returns reduced repulsion parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static Scalar reducedRepulsionParam(const Scalar& pr,
                                      const Scalar& tr) noexcept {
    return omegaB * pr / tr;
  }

  // Member functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void setParams(const Scalar& pc, const Scalar& tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    ac_ = criticalAttractionParam(pc, tc);
    bc_ = criticalRepulsionParam(pc, tc);
  }

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  Scalar reducedPressure(const Scalar& p) const noexcept { return p / pc_; }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature
  Scalar reducedTemperature(const Scalar& t) const noexcept { return t / tc_; }

 protected:
  /// @brief Get reference to derived class object
  Derived& derived() noexcept { return static_cast<Derived&>(*this); }

  /// @brief Get const reference to derived class object
  const Derived& derived() const noexcept {
    return static_cast<const Derived&>(*this);
  }

  Scalar pc_;  /// Critical pressure
  Scalar tc_;  /// Critical temperature
  Scalar ac_;  /// Critical attraction parameter
  Scalar bc_;  /// Critical repulsion parameter
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
///    - omegaA: Constant for attraction parameter
///    - omegaB: Constant for repulsion parameter
///
template <typename Derived, bool UseTemperatureCorrectionFactor>
class CubicEosBase : public CubicEosCrtpBase<Derived> {
 public:
  using Base = CubicEosCrtpBase<Derived>;
  using Scalar = typename CubicEosTraits<Derived>::Scalar;
  static constexpr auto omegaA = Base::omegaA;
  static constexpr auto omegaB = Base::omegaB;

  // Constructors

  CubicEosBase() = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosBase(const Scalar& pc, const Scalar& tc) noexcept : Base{pc, tc} {}

  CubicEosBase(const CubicEosBase&) = default;
  CubicEosBase(CubicEosBase&&) = default;

  CubicEosBase& operator=(const CubicEosBase&) = default;
  CubicEosBase& operator=(CubicEosBase&&) = default;

  // Member functions

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  IsothermalLine<Derived> createIsothermalLine(const Scalar& t) const noexcept {
    const auto tr = this->reducedTemperature(t);
    const auto alpha = this->derived().alpha(tr);
    return {t, alpha * ac_, bc_};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  IsobaricIsothermalState<Derived, UseTemperatureCorrectionFactor>
  createIsobaricIsothermalState(const Scalar& p,
                                const Scalar& t) const noexcept {
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
  Scalar pressure(const Scalar& t, const Scalar& v) const noexcept {
    const auto tr = this->reducedTemperature(t);
    const auto a = this->derived().alpha(tr) * this->attractionParam();
    const auto b = this->repulsionParam();
    return Derived::pressure(t, v, a, b);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  template <typename CubicEquationSolver>
  std::vector<Scalar> zfactor(
      const Scalar& p, const Scalar& t,
      const CubicEquationSolver& solver) const noexcept {
    return this->createIsobaricIsothermalState(p, t)
        .zfactor<CubicEquationSolver>(solver);
  }
};

template <typename Derived>
class CubicEosBase<Derived, false> : public CubicEosCrtpBase<Derived> {
 public:
  using Base = CubicEosCrtpBase<Derived>;
  using Scalar = typename CubicEosTraits<Derived>::Scalar;
  static constexpr auto omegaA = Base::omegaA;
  static constexpr auto omegaB = Base::omegaB;

  // Constructors

  CubicEosBase() = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosBase(const Scalar& pc, const Scalar& tc) noexcept : Base{pc, tc} {}

  CubicEosBase(const CubicEosBase&) = default;
  CubicEosBase(CubicEosBase&&) = default;

  CubicEosBase& operator=(const CubicEosBase&) = default;
  CubicEosBase& operator=(CubicEosBase&&) = default;

  // Member functions

  /// @brief Creates isothermal state
  /// @param[in] t Temperature
  IsothermalLine<Derived> createIsothermalLine(const Scalar& t) const noexcept {
    return {t, ac_, bc_};
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  IsobaricIsothermalState<Derived, false> createIsobaricIsothermalState(
      const Scalar& p, const Scalar& t) const noexcept {
    const auto pr = this->reducedPressure(p);
    const auto tr = this->reducedTemperature(t);
    const auto ar = this->reducedAttractionParam(pr, tr);
    const auto br = this->reducedRepulsionParam(pr, tr);
    return {t, ar, br};
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  Scalar pressure(const Scalar& t, const Scalar& v) const noexcept {
    return Derived::pressure(t, v, ac_, bc_);
  }

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  template <typename CubicEquationSolver>
  std::vector<Scalar> zfactor(
      const Scalar& p, const Scalar& t,
      const CubicEquationSolver& solver) const noexcept {
    return this->createIsobaricIsothermalState(p, t)
        .zfactor<CubicEquationSolver>(solver);
  }
};

}  // namespace eos