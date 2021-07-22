#pragma once

namespace eos {

template <typename Eos>
struct CubicEosTraits;

template <typename Eos, bool UseTemperatureCorrectionFactor>
class IsobaricIsothermalState {
 public:
  using Scalar = typename CubicEosTraits<Eos>::Scalar;

  /// @param[in] t Temperature
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  /// @param[in] beta The derivative of temperature correction factor
  IsobaricIsothermalState(const Scalar& t, const Scalar& ar, const Scalar& br,
                            const Scalar& beta) noexcept
      : t_{t}, ar_{ar}, br_{br}, beta_{beta} {}

  IsobaricIsothermalState() = default;
  IsobaricIsothermalState(const IsobaricIsothermalState &) = default;
  IsobaricIsothermalState(IsobaricIsothermalState &&) = default;

  IsobaricIsothermalState &operator=(const IsobaricIsothermalState &) =
      default;
  IsobaricIsothermalState &operator=(IsobaricIsothermalState &&) = default;

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] s Isobaric-isothermal state
  /// @return A list of Z-factors
  template <typename CubicEquationSolver>
  std::vector<Scalar> zfactor(const CubicEquationSolver& solver) const noexcept {
    return solver(Eos::zfactorCubicEq(ar_, br_));
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  Scalar lnFugacityCoeff(const Scalar& z) const noexcept {
    return Eos::lnFugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  Scalar fugacityCoeff(const Scalar& z) const noexcept {
    return Eos::fugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  Scalar residualEnthalpy(const Scalar& z) const noexcept {
    return Eos::residualEnthalpy(z, t_, ar_, br_, beta_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  Scalar residualEntropy(const Scalar& z) const noexcept {
    return Eos::residualEntropy(z, ar_, br_, beta_);
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  Scalar residualHelmholtzEnergy(const Scalar& z) const noexcept {
    return Eos::residualHelmholtzEnergy(z, t_, ar_, br_);
  }

 private:
  Scalar t_;     /// Temperature
  Scalar ar_;    /// Reduced attraction parameter
  Scalar br_;    /// Reduced repulsion parameter
  Scalar beta_;  /// The derivative of temperature correction factor for
                 /// attraction parameter
};

template <typename Eos>
class IsobaricIsothermalState<Eos, false> {
 public:
  using Scalar = typename CubicEosTraits<Eos>::Scalar;
  
  /// @param[in] t Temperature
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  IsobaricIsothermalState(const Scalar& t, const Scalar& ar, const Scalar& br) noexcept
      : t_{t}, ar_{ar}, br_{br} {}

  IsobaricIsothermalState() = default;
  IsobaricIsothermalState(const IsobaricIsothermalState &) = default;
  IsobaricIsothermalState(IsobaricIsothermalState &&) = default;

  IsobaricIsothermalState &operator=(const IsobaricIsothermalState &) =
      default;
  IsobaricIsothermalState &operator=(IsobaricIsothermalState &&) = default;

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] s Isobaric-isothermal state
  /// @return A list of Z-factors
  template <typename CubicEquationSolver>
  std::vector<Scalar> zfactor(const CubicEquationSolver& solver) const noexcept {
    return solver(Eos::zfactorCubicEq(ar_, br_));
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  Scalar lnFugacityCoeff(const Scalar& z) const noexcept {
    return Eos::lnFugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  Scalar fugacityCoeff(const Scalar& z) const noexcept {
    return Eos::fugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  Scalar residualEnthalpy(const Scalar& z) const noexcept {
    return Eos::residualEnthalpy(z, t_, ar_, br_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  Scalar residualEntropy(const Scalar& z) const noexcept {
    return Eos::residualEntropy(z, ar_, br_);
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  Scalar residualHelmholtzEnergy(const Scalar& z) const noexcept {
    return Eos::residualHelmholtzEnergy(z, t_, ar_, br_);
  }

 private:
  Scalar t_;   /// Temperature
  Scalar ar_;  /// Reduced attraction parameter
  Scalar br_;  /// Reduced repulsion parameter
};

}  // namespace eos
