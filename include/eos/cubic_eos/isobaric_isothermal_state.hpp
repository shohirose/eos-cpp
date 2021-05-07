#pragma once

namespace eos {

template <typename Eos, bool UseTemperatureCorrectionFactor>
class isobaric_isothermal_state {
 public:
  /// @param[in] t Temperature
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  /// @param[in] beta The derivative of temperature correction factor
  isobaric_isothermal_state(double t, double ar, double br,
                            double beta) noexcept
      : t_{t}, ar_{ar}, br_{br}, beta_{beta} {}

  isobaric_isothermal_state() = default;
  isobaric_isothermal_state(const isobaric_isothermal_state &) = default;
  isobaric_isothermal_state(isobaric_isothermal_state &&) = default;

  isobaric_isothermal_state &operator=(const isobaric_isothermal_state &) =
      default;
  isobaric_isothermal_state &operator=(isobaric_isothermal_state &&) = default;

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] s Isobaric-isothermal state
  /// @return A list of Z-factors
  std::vector<double> zfactor() const noexcept {
    return Eos::zfactor_cubic_eq(ar_, br_).real_roots();
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  double ln_fugacity_coeff(double z) const noexcept {
    return Eos::ln_fugacity_coeff(z, ar_, br_);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  double fugacity_coeff(double z) const noexcept {
    return Eos::fugacity_coeff(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  double residual_enthalpy(double z) const noexcept {
    return Eos::residual_enthalpy(z, t_, ar_, br_, beta_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  double residual_entropy(double z) const noexcept {
    return Eos::residual_entropy(z, ar_, br_, beta_);
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  double residual_helmholtz_energy(double z) const noexcept {
    return Eos::residual_helmholtz_energy(z, t_, ar_, br_);
  }

 private:
  double t_;     /// Temperature
  double ar_;    /// Reduced attraction parameter
  double br_;    /// Reduced repulsion parameter
  double beta_;  /// The derivative of temperature correction factor for
                 /// attraction parameter
};

template <typename Eos>
class isobaric_isothermal_state<Eos, false> {
 public:
  /// @param[in] t Temperature
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  isobaric_isothermal_state(double t, double ar, double br) noexcept
      : t_{t}, ar_{ar}, br_{br} {}

  isobaric_isothermal_state() = default;
  isobaric_isothermal_state(const isobaric_isothermal_state &) = default;
  isobaric_isothermal_state(isobaric_isothermal_state &&) = default;

  isobaric_isothermal_state &operator=(const isobaric_isothermal_state &) =
      default;
  isobaric_isothermal_state &operator=(isobaric_isothermal_state &&) = default;

  /// @brief Computes Z-factor at given pressure and temperature
  /// @param[in] s Isobaric-isothermal state
  /// @return A list of Z-factors
  std::vector<double> zfactor() const noexcept {
    return Eos::zfactor_cubic_eq(ar_, br_).real_roots();
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  double ln_fugacity_coeff(double z) const noexcept {
    return Eos::ln_fugacity_coeff(z, ar_, br_);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  double fugacity_coeff(double z) const noexcept {
    return Eos::fugacity_coeff(z, ar_, br_);
  }

  /// @brief Computes residual enthalpy
  /// @param[in] z Z-factor
  double residual_enthalpy(double z) const noexcept {
    return Eos::residual_enthalpy(z, t_, ar_, br_);
  }

  /// @brief Computes residual entropy
  /// @param[in] z Z-factor
  double residual_entropy(double z) const noexcept {
    return Eos::residual_entropy(z, ar_, br_);
  }

  /// @brief Computes residual Helmholtz energy
  /// @param[in] z Z-factor
  double residual_helmholtz_energy(double z) const noexcept {
    return Eos::residual_helmholtz_energy(z, t_, ar_, br_);
  }

 private:
  double t_;   /// Temperature
  double ar_;  /// Reduced attraction parameter
  double br_;  /// Reduced repulsion parameter
};

}  // namespace eos
