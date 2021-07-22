#pragma once

namespace eos {

template <typename Eos, bool UseTemperatureCorrectionFactor>
class IsobaricIsothermalState {
 public:
  /// @param[in] t Temperature
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  /// @param[in] beta The derivative of temperature correction factor
  IsobaricIsothermalState(double t, double ar, double br,
                            double beta) noexcept
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
  std::vector<double> zfactor() const noexcept {
    return Eos::zfactorCubicEq(ar_, br_).real_roots();
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  double lnFugacityCoeff(double z) const noexcept {
    return Eos::lnFugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  double fugacityCoeff(double z) const noexcept {
    return Eos::fugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  double residualEnthalpy(double z) const noexcept {
    return Eos::residualEnthalpy(z, t_, ar_, br_, beta_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  double residualEntropy(double z) const noexcept {
    return Eos::residualEntropy(z, ar_, br_, beta_);
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  double residualHelmholtzEnergy(double z) const noexcept {
    return Eos::residualHelmholtzEnergy(z, t_, ar_, br_);
  }

 private:
  double t_;     /// Temperature
  double ar_;    /// Reduced attraction parameter
  double br_;    /// Reduced repulsion parameter
  double beta_;  /// The derivative of temperature correction factor for
                 /// attraction parameter
};

template <typename Eos>
class IsobaricIsothermalState<Eos, false> {
 public:
  /// @param[in] t Temperature
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  IsobaricIsothermalState(double t, double ar, double br) noexcept
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
  std::vector<double> zfactor() const noexcept {
    return Eos::zfactorCubicEq(ar_, br_).real_roots();
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  double lnFugacityCoeff(double z) const noexcept {
    return Eos::lnFugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  double fugacityCoeff(double z) const noexcept {
    return Eos::fugacityCoeff(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  double residualEnthalpy(double z) const noexcept {
    return Eos::residualEnthalpy(z, t_, ar_, br_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  double residualEntropy(double z) const noexcept {
    return Eos::residualEntropy(z, ar_, br_);
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  double residualHelmholtzEnergy(double z) const noexcept {
    return Eos::residualHelmholtzEnergy(z, t_, ar_, br_);
  }

 private:
  double t_;   /// Temperature
  double ar_;  /// Reduced attraction parameter
  double br_;  /// Reduced repulsion parameter
};

}  // namespace eos
