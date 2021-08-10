#pragma once

#include <utility>  // std::pair
#include <vector>   // std::vector

#include "eos/IdentityCorrectionFactor.hpp"  // eos::IdentityCorrectionFactor
#include "eos/ThermodynamicConstants.hpp"    // eos::gasConstant

namespace eos {

/**
 * @brief Two-parameter cubic equation of state
 *
 * @tparam Scalar scalar
 * @tparam EosPolicy EoS policy
 * @tparam CorrectionFactor correction factor for attraction parameter
 */
template <typename Scalar, typename EosPolicy, typename CorrectionFactor>
class CubicEosBase {
 public:
  /// Coefficient for attraction parameter
  static constexpr auto omegaA = EosPolicy::omegaA;
  /// Coefficient for repulsion parameter
  static constexpr auto omegaB = EosPolicy::omegaB;

  CubicEosBase() = default;
  CubicEosBase(const CubicEosBase&) = default;
  CubicEosBase(CubicEosBase&&) = default;

  /**
   * @brief Construct a new Cubic Eos Base object
   *
   * @param[in] pc critical pressure
   * @param[in] tc critical temperature
   * @param[in] alpha correction factor for attraction parameter
   */
  CubicEosBase(const Scalar& pc, const Scalar& tc,
               const CorrectionFactor& alpha) noexcept
      : pc_{pc},
        tc_{tc},
        a_{attractionParam(pc, tc)},
        b_{repulsionParam(pc, tc)},
        alpha_{alpha} {}

  CubicEosBase& operator=(const CubicEosBase&) = default;
  CubicEosBase& operator=(CubicEosBase&&) = default;

  /**
   * @brief Set the critical parameters
   *
   * @param[in] pc critical pressure
   * @param[in] tc critical temperature
   */
  void setCriticalParams(const Scalar& pc, const Scalar& tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    a_ = attractionParam(pc, tc);
    b_ = repulsionParam(pc, tc);
  }

  CorrectionFactor& correctionFactor() noexcept { return alpha_; }
  const CorrectionFactor& correctionFactor() const noexcept { return alpha_; }

  /**
   * @brief Compute reduced pressure
   *
   * @param[in] p pressure
   * @return Scalar reduced pressure
   */
  Scalar reducedPressure(const Scalar& p) const noexcept { return p / pc_; }

  /**
   * @brief Compute reduced temperature
   *
   * @param t temperature
   * @return Scalar reduced temperature
   */
  Scalar reducedTemperature(const Scalar& t) const noexcept { return t / tc_; }

  /**
   * @brief Compute pressure at given temperature and volume
   *
   * @param[in] t temperature
   * @param[in] v volume
   * @return Scalar pressure
   */
  Scalar pressure(const Scalar& t, const Scalar& v) const noexcept {
    const auto tr = this->reducedTemperature(t);
    return EosPolicy::pressure(t, v, a_, b_, alpha_.value(tr));
  }

  /**
   * @brief Thermodynamic state parameters
   *
   * Data members are used to compute thermodynamic properties.
   */
  struct StateParams {
    Scalar pressure;                ///< presure
    Scalar temperature;             ///< temperature
    Scalar reducedPressure;         ///< reduced pressure
    Scalar reducedTemperature;      ///< reduced temperature
    Scalar reducedAttractionParam;  ///< reduced attraction parameter
    Scalar reducedRepulsionParam;   ///< reduced repulsion parameter
  };

  /**
   * @brief Compute Z-factor at given pressure and temperature
   *
   * @param[in] p pressure
   * @param[in] t temperature
   * @return std::pair<std::vector<Scalar>, StateParams> a pair of Z-factors and
   *    thermodynamic state parameters
   */
  template <typename CubicEquationSolver>
  std::pair<std::vector<Scalar>, StateParams> zfactor(
      const Scalar& p, const Scalar& t,
      const CubicEquationSolver& solver) const noexcept {
    const auto pr = this->reducedPressure(p);
    const auto tr = this->reducedTemperature(t);
    const auto ar = this->reducedAttractionParam(pr, tr);
    const auto br = this->reducedRepulsionParam(pr, tr);
    if constexpr (std::is_same_v<CorrectionFactor, IdentityCorrectionFactor>) {
      return {solver(EosPolicy::zfactorCubicEq(ar, br)),
              {p, t, pr, tr, ar, br}};
    } else {
      const auto ar2 = alpha_.value(tr) * ar;
      return {solver(EosPolicy::zfactorCubicEq(ar2, br)),
              {p, t, pr, tr, ar2, br}};
    }
  }

  /**
   * @brief Compute the natural logarithm of a fugacity coefficient
   *
   * @param[in] z Z-factor
   * @param[in] params thermodynamic state parameters
   * @return Scalar the natural log of fugacity coefficient
   */
  Scalar lnFugacityCoeff(const Scalar& z,
                         const StateParams& params) const noexcept {
    return EosPolicy::lnFugacityCoeff(z, params.reducedAttractionParam,
                                      params.reducedRepulsionParam);
  }

  /**
   * @brief Compute fugacity coefficient
   *
   * @param[in] z Z-factor
   * @param[in] params thermodynamic state parameters
   * @return Scalar fugacity coefficient
   */
  Scalar fugacityCoeff(const Scalar& z,
                       const StateParams& params) const noexcept {
    using std::exp;
    return exp(this->lnFugacityCoeff(z, params));
  }

  /**
   * @brief Compute residual enthalpy
   *
   * @param[in] z Z-factor
   * @param[in] params thermodynamic state parameters
   * @return Scalar residual enthalpy
   */
  Scalar residualEnthalpy(const Scalar& z,
                          const StateParams& params) const noexcept {
    if constexpr (std::is_same_v<VanDerWaalsEosPolicy<Scalar>, EosPolicy>) {
      return EosPolicy::residualEnthalpy(z, params.temperature,
                                         params.reducedAttractionParam,
                                         params.reducedRepulsionParam)
    } else {
      const auto dlnAlpha_dlnT = alpha_.derivative(params.reducedTemperature);
      return EosPolicy::residualEnthalpy(
          z, params.temperature, params.reducedAttractionParam,
          params.reducedRepulsionParam, dlnAlpha_dlnT);
    }
  }

  /**
   * @brief Compute residual entropy
   *
   * @param[in] z Z-factor
   * @param[in] params thermodynamic state parameters
   * @return Scalar residual entropy
   */
  Scalar residualEntropy(const Scalar& z,
                         const StateParams& params) const noexcept {
    if constexpr (std::is_same_v<EosPolicy, VanDerWaalsEosPolicy<Scalar>>) {
      return EosPolicy::residualEntropy(z, params.reducedAttractionParam,
                                        params.reducedRepulsionParam)
    } else {
      const auto dlnAlpha_dlnT = alpha_.derivative(params.reducedTemperature);
      return EosPolicy::residualEntropy(z, params.reducedAttractionParam,
                                        params.reducedRepulsionParam,
                                        dlnAlpha_dlnT);
    }
  }

  /**
   * @brief Compute residual Helmholtz energy
   *
   * @param[in] z Z-factor
   * @param[in] params thermodynamic state parameters
   * @return Scalar residual Helmholtz energy
   */
  Scalar residualHelmholtzEnergy(const Scalar& z,
                                 const StateParams& params) const noexcept {
    return EosPolicy::residualHelmholtzEnergy(z, params.temperature,
                                              params.reducedAttractionParam,
                                              params.reducedRepulsionParam);
  }

 private:
  /// @name Private static functions
  //@{

  /**
   * @brief Compute attraction parameter
   *
   * @param[in] pc critical pressure
   * @param[in] tc critical temperature
   * @return Scalar attraction parameter
   */
  static Scalar attractionParam(const Scalar& pc, const Scalar& tc) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return (omegaA * R * R) * tc * tc / pc;
  }

  /**
   * @brief Compute repulsion parameter
   *
   * @param[in] pc critical pressure
   * @param[in] tc critical temperature
   * @return Scalar repulsion parameter
   */
  static Scalar repulsionParam(const Scalar& pc, const Scalar& tc) noexcept {
    constexpr auto R = gasConstant<Scalar>();
    return (omegaB * R) * tc / pc;
  }

  /**
   * @brief Compute reduced attraction parameter
   *
   * @param[in] pr reduced pressure
   * @param[in] tr reduced temperature
   * @return Scalar reduced attraction parameter
   */
  static Scalar reducedAttractionParam(const Scalar& pr,
                                       const Scalar& tr) noexcept {
    return omegaA * pr / (tr * tr);
  }

  /**
   * @brief Compute reduced repulsion parameter
   *
   * @param[in] pr reduced pressure
   * @param[in] tr reduced temperature
   * @return Scalar reduced repulsion parameter
   */
  static Scalar reducedRepulsionParam(const Scalar& pr,
                                      const Scalar& tr) noexcept {
    return omegaB * pr / tr;
  }

  Scalar pc_;               ///< Critical pressure
  Scalar tc_;               ///< Critical temperature
  Scalar a_;                ///< Attraction parameter
  Scalar b_;                ///< Repulsion parameter
  CorrectionFactor alpha_;  ///< Correction factor for attraction parameter
};

}  // namespace eos