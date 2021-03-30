#pragma once

#include <cmath> // std::log, std::exp, std::fabs, std::pow

namespace eos
{

  template <typename T>
  class lucas_method
  {
  public:
    using scalar_type = T;

    lucas_method() = default;

    /// @brief Constructs object
    /// @param[in] pc Critical pressure [Pa]
    /// @param[in] tc Critical temprature [K]
    /// @param[in] zc Critical z-factor
    /// @param[in] mw Molecular weight [kg/kmol]
    /// @param[in] dm Dipole moment [Debyes]
    /// @param[in] q Quantum parameter
    lucas_method(const scalar_type &pc, const scalar_type &tc, const scalar_type &zc,
                 const scalar_type &mw, const scalar_type &dm, const scalar_type &q)
        : pc_{pc},
          tc_{tc},
          zc_{zc},
          mw_{mw},
          dm_{dm},
          dmr_{reduced_dipole_moment(dm, tc, pc)},
          q_{q},
          xi_{inverse_viscosity(pc, tc, mw)} {}

    // Member functions

    void set_params(const scalar_type &pc, const scalar_type &tc, const scalar_type &zc,
                    const scalar_type &mw, const scalar_type &dm, const scalar_type &q) noexcept
    {
      pc_ = pc;
      tc_ = tc;
      zc_ = zc;
      mw_ = mw;
      dm_ = dm;
      dmr_ = reduced_dipole_moment(dm, tc, pc);
      q_ = q;
      xi_ = inverse_viscosity(pc, tc, mw)
    }

    /// @brief Computes gas viscosity at low pressure
    /// @param[in] t Temperature [K]
    /// @return Viscosity [Pa-s]
    scalar_type viscosity_at_low_pressure(const scalar_type &t) const noexcept
    {
      const auto tr = this->reduced_temperature(t);
      const auto fp = this->polarity_factor_at_low_pressure(tr);
      const auto fq = this->quantum_factor_at_low_pressure(tr);
      const auto z1 = this->reduced_viscosity_at_low_pressure(tr, fp, fq);
      return z1 / xi_;
    }

    /// @brief Computes gas viscosity at high pressure
    /// @param[in] p Pressure [Pa]
    /// @param[in] t Temperature [K]
    /// @return Viscosity [Pa-s]
    scalar_type viscosity_at_high_pressure(const scalar_type &p, const scalar_type &t) const noexcept
    {
      const auto pr = this->reduced_pressure(p);
      const auto tr = this->reduced_temperature(t);
      const auto fp0 = this->polarity_factor_at_low_pressure(tr);
      const auto fq0 = this->quantum_factor_at_low_pressure(tr);

      const auto z1 = this->reduced_viscosity_at_low_pressure(tr, fp0, fq0);
      const auto z2 = this->reduced_viscosity_at_high_pressure(z1, pr, tr);
      const auto fp = this->polarity_factor_at_high_pressure(fp0, z1, z2);
      const auto fq = this->quantum_factor_at_high_pressure(fq0, z1, z2);

      return z2 * fp * fq / xi_;
    }

  private:
    /// @brief Computes reduced dipole moment
    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    /// @param[in] dm Dipole moment
    static scalar_type reduced_dipole_moment(const scalar_type &pc, const scalar_type &tc, const scalar_type &dm) noexcept
    {
      return 52.46e-5 * dm * dm / (tc * tc) * pc;
    }

    /// @brief Computes reduced viscosity at low pressure.
    /// @param[in] tr Reduced temperature
    /// @param[in] fp Polarity factor at low pressure
    /// @param[in] fq Quantum factor at low pressure
    static scalar_type reduced_viscosity_at_low_pressure(const scalar_type &tr, const scalar_type &fp, const scalar_type &fq) noexcept
    {
      using std::exp;
      using std::pow;
      const auto z1 = (0.807 * pow(tr, 0.618) - 0.357 * exp(-0.449 * tr) + 0.340 * exp(-4.058 * tr) + 0.018);
      return z1 * fp * fq;
    }

    /// @brief Computes inverse viscosity.
    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    /// @param[in] mw Molecular weight
    static scalar_type inverse_viscosity(const scalar_type &pc, const scalar_type &tc, const scalar_type &mw) noexcept
    {
      using std::pow;
      const auto pc_bar = pc * 1e-5;
      const auto pc2 = pc_bar * pc_bar;
      return (1.0e7 * 0.176) * pow(tc / (mw * mw * mw * pc2 * pc2), 1.0 / 6.0);
    }

    /// @brief Computes reduced viscosity at high pressure
    /// @param[in] z1 Reduced viscosity at low pressure
    /// @param[in] pr Reduced pressure
    /// @param[in] tr Reduced temperature
    /// @return Reduced viscosity
    static scalar_type reduced_viscosity_at_high_pressure(const scalar_type &z1, const scalar_type &pr, const scalar_type &tr) noexcept
    {
      using std::exp;
      using std::pow;

      if (tr <= static_cast<scalar_type>(1.0))
      {
        const auto alpha = 3.262 + 14.98 * pow(pr, 5.508);
        const auto beta = 1.390 + 5.746 * pr;
        return 0.600 + 0.760 * pow(pr, alpha) + (6.990 * pow(pr, beta) - 0.6) * (1.0 - tr);
      }
      else
      {
        const auto a = 1.245e-3 / tr * exp(5.1726 * pow(tr, -0.3286));
        const auto b = a * (1.6553 * tr - 1.2723);
        const auto c = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
        const auto d = 1.7369 / tr * exp(2.2310 * pow(tr, -7.6351));
        const auto f = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
        const auto z2 = 1.0 + a * pow(pr, 1.3088) / (b * pow(pr, f) + 1.0 / (1.0 + c * pow(pr, d)));
        return z2 * z1;
      }
    }

    /// @brief Computes polarity factor at high pressure.
    /// @param[in] fp0 Polarity factor at low pressure
    /// @param[in] z1 Reduced viscosity at low pressure
    /// @param[in] z2 Reduced viscosity at high pressure
    static scalar_type polarity_factor_at_high_pressure(const scalar_type &fp0, const scalar_type &z1, const scalar_type &z2) noexcept
    {
      const auto y = z2 / z1;
      return (1.0 + (fp0 - 1.0) / (y * y * y)) / fp0;
    }

    /// @brief Computes quantum factor at high pressure.
    /// @param[in] fq0 Quantum factor at low pressure
    /// @param[in] z1 Reduced viscosity at low pressure
    /// @param[in] z2 Reduced viscosity at high pressure
    static scalar_type quantum_factor_at_high_pressure(const scalar_type &fq0, const scalar_type &z1, const scalar_type &z2) noexcept
    {
      const auto y = z2 / z1;
      using std::log;
      const auto tmp = log(y);
      const auto tmp2 = tmp * tmp;
      return (1.0 + (fq0 - 1.0) * (1.0 / y - 0.007 * tmp2 * tmp2)) / fq0;
    }

    /// @brief Computes reduced pressure
    /// @param[in] p Pressure
    scalar_type reduced_pressure(const scalar_type &p) const noexcept { return p / pc_; }

    /// @brief Computes reduced temperature
    /// @param[in] t Temperature
    scalar_type reduced_temperature(const scalar_type &t) const noexcept { return t / tc_; }

    /// @brief Computes polarity factor at low pressure.
    /// @param[in] tr Reduced temperature
    scalar_type polarity_factor_at_low_pressure(const scalar_type &tr) const noexcept
    {
      using std::fabs;
      using std::pow;
      assert(dmr_ >= 0);
      if (dmr_ < static_cast<scalar_type>(0.022))
      {
        return 1.0;
      }
      else if (dmr_ < static_cast<scalar_type>(0.075))
      {
        return 1.0 + 30.55 * pow(0.292 - zc_, 1.72);
      }
      else
      {
        return 1.0 + 30.55 * pow(0.292 - zc_, 1.72) * fabs(0.96 + 0.1 * (tr - 0.7));
      }
    }

    /// @brief Computes quantum factor at low pressure.
    /// @param[in] tr Reduced temperature
    scalar_type quantum_factor_at_low_pressure(const scalar_type &tr) const noexcept
    {
      using std::copysign;
      using std::fabs;
      using std::pow;
      constexpr auto tolerance = static_cast<scalar_type>(1e-10);
      assert(q_ >= 0);
      if (fabs(q_) < tolerance)
      {
        return 1.0;
      }
      else if (fabs(tr - 12) < tolerance)
      {
        return 1.22 * pow(q_, 0.15);
      }
      else
      {
        const auto tmp = tr - 12;
        return 1.22 * pow(q_, 0.15) * (1.0 + copysign(0.00385 * pow(tmp * tmp, 1.0 / mw_), tmp));
      }
    }

    scalar_type pc_;  /// Critical pressure [Pa]
    scalar_type tc_;  /// Critical temperature [K]
    scalar_type zc_;  /// Critical Z-factor
    scalar_type mw_;  /// Molecular weight [kg/kmol]
    scalar_type dm_;  /// Dipole moment [Debyes]
    scalar_type dmr_; /// Reduced dipole moment

    /// @brief Quantum parameter
    ///
    /// Quantum factor is required only for quantum gases: H2, He, and D2.
    /// q = 1.38 (He), q = 0.76 (H2), and q = 0.52 (D2).
    scalar_type q_;
    scalar_type xi_; /// Inverse viscosity
  };

  template <typename T>
  lucas_method<T> make_lucas_method(const T &pc, const T &tc, const T &zc, const T &mw, const T &dm, const T &q)
  {
    return {pc, tc, zc, mw, dm, q};
  }

} // namespace eos