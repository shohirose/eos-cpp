#pragma once

#include <array>                              // std::array
#include <boost/math/constants/constants.hpp> // boost::math::constants::root_two
#include <cmath>                              // std::sqrt, std::exp, std::log

#include "eos/cubic_eos_base.hpp" // eos::cutic_eos_base

namespace eos
{

  template <typename scalar_type>
  class peng_robinson_eos;

  namespace detail
  {

    template <typename T>
    struct cubic_eos_traits<peng_robinson_eos<T>>
    {
      using scalar_type = T;
      static constexpr auto omega_a = static_cast<scalar_type>(0.45724);
      static constexpr auto omega_b = static_cast<scalar_type>(0.07780);
    };

  } // namespace detail

  /// @brief Peng-Robinson EoS.
  /// @tparam T Value type
  template <typename T>
  class peng_robinson_eos : public cutic_eos_base<peng_robinson_eos<T>>
  {
  public:
    using scalar_type = T;
    using base_type = cutic_eos_base<peng_robinson_eos<T>>;

    // Static functions

    /// @brief Computes pressure at given temperature and volume
    /// @param[in] t Temperature
    /// @param[in] v Volume
    /// @param[in] a Attraction parameter
    /// @param[in] b Repulsion parameter
    /// @returns Pressure
    static scalar_type pressure_impl(const scalar_type &t, const scalar_type &v, const scalar_type &a, const scalar_type &b) noexcept
    {
      constexpr auto R = gas_constant<scalar_type>();
      return R * t / (v - b) - a / (v * (v + b) + b * (v - b));
    }

    /// @brief Computes coeficients of the cubic equation of z-factor.
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @returns Coefficients of the cubic equation of z-factor.
    static std::array<scalar_type, 3> cubic_eq_impl(const scalar_type &a, const scalar_type &b) noexcept
    {
      return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
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
      return exp(z - 1 - log(z - b) - q(z, a, b));
    }

    /// @brief Computes residual enthalpy
    /// @param[in] z Z-factor
    /// @param[in] t Temperature
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @param[in] beta Temperature correction factor
    static scalar_type residual_enthalpy_impl(const scalar_type &z, const scalar_type &t, const scalar_type &a, const scalar_type &b, const scalar_type &beta) noexcept
    {
      using std::log;
      constexpr auto R = gas_constant<scalar_type>();
      return R * t * (z - 1 - (1 - beta) * q(z, a, b));
    }

    /// @brief Computes residual entropy
    /// @param[in] z Z-factor
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    /// @param[in] beta Temperature correction factor
    static scalar_type residual_entropy_impl(const scalar_type &z, const scalar_type &a, const scalar_type &b, const scalar_type &beta) noexcept
    {
      using std::log;
      constexpr auto R = gas_constant<scalar_type>();
      return R * (log(z - b) + beta * q(z, a, b));
    }

    // Constructors

    peng_robinson_eos() = default;

    /// @brief Constructs Peng-Robinson EoS
    /// @param[in] pc Critical pressrue
    /// @param[in] tc Critical temperature
    /// @param[in] omega Acentric factor
    peng_robinson_eos(const scalar_type &pc, const scalar_type &tc, const scalar_type &omega)
        : base_type{pc, tc}, omega_{omega}, m_{m(omega)} {}

    // Member functions

    /// @brief Set parameters
    /// @param[in] pc Critical pressrue
    /// @param[in] tc Critical temperature
    /// @param[in] omega Acentric factor
    void set_params(const scalar_type &pc, const scalar_type &tc, const scalar_type &omega) noexcept
    {
      this->base_type::set_params(pc, tc);
      omega_ = omega;
      m_ = m(omega);
    }

    // Member functions

    /// @brief Computes the correction factor for attraction parameter
    /// @param[in] tr Reduced temperature
    scalar_type alpha(const scalar_type &tr) const noexcept
    {
      using std::sqrt;
      const auto a = 1 + m_ * (1 - sqrt(tr));
      return a * a;
    }

    /// @brief Computes \f$ \beta = \frac{\mathrm{d} \ln \alpha}{\mathrm{d} \ln
    /// scalar_type} \f$
    /// @param[in] tr Reduced temperature
    scalar_type beta(const scalar_type &tr) const noexcept
    {
      using std::sqrt;
      const auto sqrt_tr = sqrt(tr);
      const auto a = 1 + m_ * (1 - sqrt_tr);
      return -m_ * sqrt_tr / a;
    }

  private:
    /// @brief Computes parameter \f$ m \f$ from acentric factor
    /// @param[in] omega Acentric factor
    static scalar_type m(const scalar_type &omega) noexcept
    {
      return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
    }

    /// @brief Computes a coefficient appearing in the calculation of fugacity coefficient,
    /// residual enthalpy and entropy.
    /// @param[in] z Z-factor
    /// @param[in] a Reduced attraction parameter
    /// @param[in] b Reduced repulsion parameter
    static scalar_type q(const scalar_type &z, const scalar_type &a, const scalar_type &b) noexcept
    {
      constexpr auto sqrt2 = boost::math::constants::root_two<scalar_type>();
      constexpr auto delta1 = 1 + sqrt2;
      constexpr auto delta2 = 1 - sqrt2;
      using std::log;
      return a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z + delta2 * b));
    }

    /*
    /// @brief Computes  \f[ \gamma = \frac{T_r^2}{\alpha} \cdot
    /// \frac{\mathrm{d}^2 \alpha}{\mathrm{d} T_r^2} \f]
    /// @param[in] tr Reduced temperature
    scalar_type gamma(const scalar_type &tr) const noexcept {
      using std::sqrt;
      const auto sqrt_tr = sqrt(tr);
      const auto a = 1 + m_ * (1 - sqrt_tr);
      const auto alpha = a * a;
      return m_ / (2 * alpha) * (m_ * tr + a * sqrt_tr);
    }
    */

    /// Acentric factor
    scalar_type omega_;
    /// \f$ m = 0.3796 + 1.485 \omega - 0.1644 \omega^2 + 0.01667 \omega^3 \f$
    scalar_type m_;
  };

  /// @brief Makes Peng-Robinson EoS
  /// @tparam T Scalar type
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  template <typename T>
  inline peng_robinson_eos<T> make_peng_robinson_eos(const T &pc, const T &tc, const T &omega)
  {
    return {pc, tc, omega};
  }

} // namespace eos
