#pragma once

#include <array> // std::array
#include <cmath> // std::exp, std::log

#include "eos/cubic_eos_base.hpp" // eos::cutic_eos_base

namespace eos
{

  template <typename T>
  class van_der_waals_eos;

  namespace detail
  {

    template <typename T>
    struct cubic_eos_traits<van_der_waals_eos<T>>
    {
      using scalar_type = T;
      static constexpr auto omega_a = static_cast<T>(0.421875);
      static constexpr auto omega_b = static_cast<T>(0.125);
    };

  } // namespace detail

  /// @brief Van der Waals Equations of State
  template <typename T>
  class van_der_waals_eos : public cutic_eos_base<van_der_waals_eos<T>>
  {
  public:
    using scalar_type = T;
    using base_type = cutic_eos_base<van_der_waals_eos<T>>;

    // Static Functions

    /// @brief Computes pressure at given temperature and volume.
    /// @param[in] t Temperature
    /// @param[in] v Volume
    /// @param[in] a Attraction parameter
    /// @param[in] b Repulsion parameter
    /// @returns Pressure
    static scalar_type pressure_impl(const scalar_type &t, const scalar_type &v, const scalar_type &a, const scalar_type &b) noexcept
    {
      return gas_constant<scalar_type>() * t / (v - b) - a / (v * v);
    }

    /// @brief Computes coefficients of cubic equation
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @returns Coefficients of the cubic equation of z-factor
    static std::array<scalar_type, 3> cubic_eq_impl(const scalar_type &a, const scalar_type &b) noexcept
    {
      return {-b - 1, a, -a * b};
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @returns Fugacity coefficient
    static scalar_type fugacity_coeff_impl(const scalar_type &z, const scalar_type &a, const scalar_type &b) noexcept
    {
      using std::exp;
      using std::log;
      return exp(-log(z - b) - a / z + z - 1);
    }

    /// @brief Computes residual enthalpy
    /// @param[in] z Z-factor
    /// @param[in] t Temperature
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @param[in] beta Temperature correction factor
    static scalar_type residual_enthalpy_impl(const scalar_type &z, const scalar_type &t, const scalar_type &a, const scalar_type &b, const scalar_type &beta) noexcept
    {
      return gas_constant<scalar_type>() * t * (z - 1 - a * (1 - beta) / z);
    }

    /// @brief Computes residual entropy
    /// @param[in] z Z-factor
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @param[in] beta Temperature correction factor
    static scalar_type residual_entropy_impl(const scalar_type &z, const scalar_type &a, const scalar_type &b, const scalar_type &beta) noexcept
    {
      using std::log;
      return gas_constant<scalar_type>() * (log(z - b) + a * beta / z);
    }

    /*
    /// @brief Computes residual molar specific heat at constant volume
    /// @param[in] z Z-factor
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @param[in] gamma Temperature correction factor
    static scalar_type residual_specific_heat_at_const_volume(const scalar_type &z, const
    scalar_type &a,
                                                    [[maybe_unused]] const scalar_type
    &b, const scalar_type &gamma) noexcept { return gas_constant * gamma * a / z;
    }
    */

    van_der_waals_eos() = default;

    van_der_waals_eos(const scalar_type &pc, const scalar_type &tc) noexcept
        : base_type{pc, tc} {}

    void set_params(const scalar_type &pc, const scalar_type &tc) noexcept
    {
      this->base_type::set_params(pc, tc);
    }

    constexpr scalar_type alpha(const scalar_type &) const noexcept { return 1.0; }

    constexpr scalar_type beta(const scalar_type &) const noexcept { return 0.0; }

    // scalar_type gamma(const scalar_type &) const noexcept { return 0.0; }
  };

  /// @brief Makes van der Waals EoS
  /// @tparam T Scalar type
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  template <typename T>
  inline van_der_waals_eos<T> make_van_der_waals_eos(const T &pc, const T &tc)
  {
    return {pc, tc};
  }

} // namespace eos