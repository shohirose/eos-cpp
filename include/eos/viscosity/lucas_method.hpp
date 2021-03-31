#pragma once

#include <cmath> // std::log, std::exp, std::fabs, std::pow

namespace eos
{

  class lucas_method
  {
  public:
    lucas_method() = default;

    /// @brief Constructs object
    /// @param[in] pc Critical pressure [Pa]
    /// @param[in] tc Critical temprature [K]
    /// @param[in] zc Critical z-factor
    /// @param[in] mw Molecular weight [kg/kmol]
    /// @param[in] dm Dipole moment [Debyes]
    /// @param[in] q Quantum parameter
    lucas_method(double pc, double tc, double zc, double mw, double dm, double q)
        : pc_{pc}, tc_{tc}, zc_{zc}, mw_{mw}, dm_{dm}, q_{q},
          dmr_{reduced_dipole_moment(dm, tc, pc)},
          xi_{inverse_viscosity(pc, tc, mw)} {}

    // Member functions

    void set_params(double pc, double tc, double zc, double mw, double dm, double q) noexcept;

    /// @brief Computes gas viscosity at low pressure
    /// @param[in] t Temperature [K]
    /// @return Viscosity [Pa-s]
    double viscosity_at_low_pressure(double t) const noexcept;

    /// @brief Computes gas viscosity at high pressure
    /// @param[in] p Pressure [Pa]
    /// @param[in] t Temperature [K]
    /// @return Viscosity [Pa-s]
    double viscosity_at_high_pressure(double p, double t) const noexcept;

  private:
    /// @brief Computes reduced dipole moment
    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    /// @param[in] dm Dipole moment
    static double reduced_dipole_moment(double pc, double tc, double dm) noexcept
    {
      return 52.46e-5 * dm * dm / (tc * tc) * pc;
    }

    /// @brief Computes reduced viscosity at low pressure.
    /// @param[in] tr Reduced temperature
    /// @param[in] fp Polarity factor at low pressure
    /// @param[in] fq Quantum factor at low pressure
    static double reduced_viscosity_at_low_pressure(double tr, double fp, double fq) noexcept
    {
      const auto z1 = (0.807 * std::pow(tr, 0.618) - 0.357 * std::exp(-0.449 * tr) + 0.340 * std::exp(-4.058 * tr) + 0.018);
      return z1 * fp * fq;
    }

    /// @brief Computes inverse viscosity.
    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    /// @param[in] mw Molecular weight
    static double inverse_viscosity(double pc, double tc, double mw) noexcept
    {
      const auto pc_bar = pc * 1e-5;
      const auto pc2 = pc_bar * pc_bar;
      return (1.0e7 * 0.176) * std::pow(tc / (mw * mw * mw * pc2 * pc2), 1.0 / 6.0);
    }

    /// @brief Computes reduced viscosity at high pressure
    /// @param[in] z1 Reduced viscosity at low pressure
    /// @param[in] pr Reduced pressure
    /// @param[in] tr Reduced temperature
    /// @return Reduced viscosity
    static double reduced_viscosity_at_high_pressure(double z1, double pr, double tr) noexcept;

    /// @brief Computes polarity factor at high pressure.
    /// @param[in] fp0 Polarity factor at low pressure
    /// @param[in] z1 Reduced viscosity at low pressure
    /// @param[in] z2 Reduced viscosity at high pressure
    static double polarity_factor_at_high_pressure(double fp0, double z1, double z2) noexcept
    {
      const auto y = z2 / z1;
      return (1.0 + (fp0 - 1.0) / (y * y * y)) / fp0;
    }

    /// @brief Computes quantum factor at high pressure.
    /// @param[in] fq0 Quantum factor at low pressure
    /// @param[in] z1 Reduced viscosity at low pressure
    /// @param[in] z2 Reduced viscosity at high pressure
    static double quantum_factor_at_high_pressure(double fq0, double z1, double z2) noexcept
    {
      const auto y = z2 / z1;
      const auto tmp = std::log(y);
      const auto tmp2 = tmp * tmp;
      return (1.0 + (fq0 - 1.0) * (1.0 / y - 0.007 * tmp2 * tmp2)) / fq0;
    }

    /// @brief Computes reduced pressure
    /// @param[in] p Pressure
    double reduced_pressure(double p) const noexcept { return p / pc_; }

    /// @brief Computes reduced temperature
    /// @param[in] t Temperature
    double reduced_temperature(double t) const noexcept { return t / tc_; }

    /// @brief Computes polarity factor at low pressure.
    /// @param[in] tr Reduced temperature
    double polarity_factor_at_low_pressure(double tr) const noexcept;

    /// @brief Computes quantum factor at low pressure.
    /// @param[in] tr Reduced temperature
    double quantum_factor_at_low_pressure(double tr) const noexcept;

    double pc_;  /// Critical pressure [Pa]
    double tc_;  /// Critical temperature [K]
    double zc_;  /// Critical Z-factor
    double mw_;  /// Molecular weight [kg/kmol]
    double dm_;  /// Dipole moment [Debyes]
    double q_;   /// Quantum parameter for H2, He, and D2
    double dmr_; /// Reduced dipole moment
    double xi_;  /// Inverse viscosity
  };

  inline lucas_method make_lucas_method(double pc, double tc, double zc, double mw, double dm, double q)
  {
    return {pc, tc, zc, mw, dm, q};
  }

} // namespace eos