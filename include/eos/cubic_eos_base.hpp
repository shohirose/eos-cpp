#pragma once

#include <utility>  // std::pair
#include <vector>   // std::vector

#include "eos/thermodynamic_constants.hpp"  // eos::gasConstant

namespace eos {

/// @brief Two-parameter cubic equation of state (EoS)
template <typename Scalar, typename EosPolicy, typename CorrectionPolicy>
class CubicEosBase {
 public:
  static constexpr auto omegaA = EosPolicy::omegaA;
  static constexpr auto omegaB = EosPolicy::omegaB;

  CubicEosBase() = default;
  CubicEosBase(const CubicEosBase&) = default;
  CubicEosBase(CubicEosBase&&) = default;

  /// @brief Constructs cubic EoS
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosBase(const Scalar& pc, const Scalar& tc,
               const CorrectionPolicy& corrector) noexcept
      : pc_{pc},
        tc_{tc},
        a_{attractionParam(pc, tc)},
        b_{repulsionParam(pc, tc)},
        corrector_{corrector} {}

  CubicEosBase& operator=(const CubicEosBase&) = default;
  CubicEosBase& operator=(CubicEosBase&&) = default;

  // Member functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void setCriticalParams(const Scalar& pc, const Scalar& tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    a_ = attractionParam(pc, tc);
    b_ = repulsionParam(pc, tc);
  }

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  Scalar reducedPressure(const Scalar& p) const noexcept { return p / pc_; }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature
  Scalar reducedTemperature(const Scalar& t) const noexcept { return t / tc_; }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  Scalar pressure(const Scalar& t, const Scalar& v) const noexcept {
    if constexpr (std::is_same_v<CorrectionPolicy,
                                 NoCorrectionPolicy<Scalar>>) {
      return EosPolicy::pressure(t, v, a_, b_);
    } else {
      const auto tr = this->reducedTemperature(t);
      const auto alpha = corrector_.value(tr);
      return EosPolicy::pressure(t, v, alpha * a_, b_);
    }
  }

  struct StateParams {
    Scalar pressure;
    Scalar temperature;
    Scalar reducedPressure;
    Scalar reducedTemperature;
    Scalar reducedAttractionParam;
    Scalar reducedRepulsionParam;
  };

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @return A list of Z-factors
  template <typename CubicEquationSolver>
  std::pair<std::vector<Scalar>, StateParams> zfactor(
      const Scalar& p, const Scalar& t,
      const CubicEquationSolver& solver) const noexcept {
    const auto pr = this->reducedPressure(p);
    const auto tr = this->reducedTemperature(t);
    const auto ar = this->reducedAttractionParam(pr, tr);
    const auto br = this->reducedRepulsionParam(pr, tr);
    if constexpr (std::is_same_v<CorrectionPolicy,
                                 NoCorrectionPolicy<Scalar>>) {
      return {solver(EosPolicy::zfactorCubicEq(ar, br)),
              {p, t, pr, tr, ar, br}};
    } else {
      const auto alpha = corrector_.value(tr);
      const auto ar2 = alpha * ar;
      return {solver(EosPolicy::zfactorCubicEq(ar2, br)),
              {p, t, pr, tr, ar2, br}};
    }
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  Scalar lnFugacityCoeff(const Scalar& z,
                         const StateParams& params) const noexcept {
    return EosPolicy::lnFugacityCoeff(z, params.reducedAttractionParam,
                                      params.reducedRepulsionParam);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  Scalar fugacityCoeff(const Scalar& z,
                       const StateParams& params) const noexcept {
    using std::exp;
    return exp(this->lnFugacityCoeff(z, params));
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  Scalar residualEnthalpy(const Scalar& z,
                          const StateParams& params) const noexcept {
    if constexpr (std::is_same_v<VanDerWaalsEosPolicy<Scalar>, EosPolicy>) {
      return EosPolicy::residualEnthalpy(z, params.temperature,
                                         params.reducedAttractionParam,
                                         params.reducedRepulsionParam)
    } else {
      const auto beta = corrector_.derivative(params.reducedTemperature);
      return EosPolicy::residualEnthalpy(z, params.temperature,
                                         params.reducedAttractionParam,
                                         params.reducedRepulsionParam, beta);
    }
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  Scalar residualEntropy(const Scalar& z,
                         const StateParams& params) const noexcept {
    if constexpr (std::is_same_v<EosPolicy, VanDerWaalsEosPolicy<Scalar>>) {
      return EosPolicy::residualEntropy(z, params.reducedAttractionParam,
                                        params.reducedRepulsionParam)
    } else {
      const auto beta = corrector_.derivative(params.reducedTemperature);
      return EosPolicy::residualEntropy(z, params.reducedAttractionParam,
                                        params.reducedRepulsionParam, beta);
    }
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  Scalar residualHelmholtzEnergy(const Scalar& z,
                                 const StateParams& params) const noexcept {
    return EosPolicy::residualHelmholtzEnergy(z, params.temperature,
                                              params.reducedAttractionParam,
                                              params.reducedRepulsionParam);
  }

 protected:
  // Member functions

  const CorrectionPolicy& correctionPolicy() const noexcept {
    return corrector_;
  }

  CorrectionPolicy& correctionPolicy() noexcept { return corrector_; }

 private:
  // Static functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static Scalar attractionParam(const Scalar& pc, const Scalar& tc) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return (omegaA * R * R) * tc * tc / pc;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static Scalar repulsionParam(const Scalar& pc, const Scalar& tc) noexcept {
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

  Scalar pc_;  /// Critical pressure
  Scalar tc_;  /// Critical temperature
  Scalar a_;   /// Attraction parameter
  Scalar b_;   /// Repulsion parameter

  /// Temperature correction policy for attraction parameter
  CorrectionPolicy corrector_;
};

}  // namespace eos