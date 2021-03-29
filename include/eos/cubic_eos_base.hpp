#pragma once

#include <vector> // std::vector

#include "eos/constants.hpp" // eos::gas_constant
#include "eos/roots.hpp"     // eos::real_roots

namespace eos
{

  namespace detail
  {

    template <typename T>
    struct cubic_eos_traits
    {
    };

  } // namespace detail

  template <typename Eos>
  class isothermal_state
  {
  public:
    using scalar_type = typename detail::cubic_eos_traits<Eos>::scalar_type;

    /// @param[in] t Temperature
    /// @param[in] a Attraction parameter
    /// @param[in] b Repulsion parameter
    isothermal_state(const scalar_type &t, const scalar_type &a, const scalar_type &b) noexcept
        : t_{t}, a_{a}, b_{b} {}

    /// @brief Computes pressure at given temperature and volume
    /// @param[in] v Volume
    scalar_type pressure(const scalar_type &v) const noexcept
    {
      return Eos::pressure_impl(t_, v, a_, b_);
    }

  private:
    scalar_type t_; /// Temperature
    scalar_type a_; /// Attraction parameter
    scalar_type b_; /// Repulsion parameter
  };

  template <typename Eos>
  class isobaric_isothermal_state
  {
  public:
    using scalar_type = typename detail::cubic_eos_traits<Eos>::scalar_type;

    isobaric_isothermal_state(const scalar_type &ar, const scalar_type &br, const scalar_type &beta) noexcept
        : ar_{ar}, br_{br}, beta_{beta} {}

    /// @brief Computes Z-factor at given pressure and temperature
    /// @param[in] s Isobaric-isothermal state
    /// @return A list of Z-factors
    std::vector<scalar_type> zfactor() const noexcept
    {
      return real_roots(Eos::cubic_eq_impl(ar_, br_));
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @param[in] s Isobaric-isothermal state
    scalar_type fugacity_coeff(const scalar_type &z) const noexcept
    {
      return Eos::fugacity_coeff_impl(z, ar_, br_);
    }

    /// @brief Computes residual enthalpy
    /// @param[in] z Z-factor
    /// @param[in] s Isobaric-isothermal state
    scalar_type residual_enthalpy(const scalar_type &z) const noexcept
    {
      return Eos::residual_enthalpy_impl(z, ar_, br_, beta_);
    }

    /// @brief Computes residual entropy
    /// @param[in] z Z-factor
    /// @param[in] s Isobaric-isothermal state
    scalar_type residual_entropy(const scalar_type &z) const noexcept
    {
      return Eos::residual_entropy_impl(z, ar_, br_, beta_);
    }

  private:
    scalar_type ar_; /// Reduced attraction parameter
    scalar_type br_; /// Reduced repulsion parameter
    scalar_type beta_;
  };

  /// @brief Two-parameter cubic equation of state (EoS)
  /// @tparam Derived Concrete EoS class
  ///
  /// Derived EoS classes must have the following static functions:
  ///    - pressure_impl(t, v, a, b)
  ///    - cubic_eq_impl(ar, br)
  ///    - fugacity_coeff_impl(z, ar, br)
  ///    - residual_enthalpy_impl(z, t, ar, br, beta)
  ///    - residual_entropy_impl(z, ar, br, beta)
  ///    - alpha(tr)
  ///    - beta()
  /// where t is temperature, v is volume, a is attraction parameter, b is
  /// repulsion parameter, ar is reduced attraction parameter, br is reduced
  /// repulsion parameter, z is z-factor, and beta is a temperature correction
  /// factor.
  ///
  /// cubic_eos_traits class specialized for each concrete EoS class must be defined
  /// in the detail namespace. cubic_eos_traits class must define the following types
  /// and constants:
  ///    - scalar_type: scalar_type type
  ///    - omega_a: Constant for attraction parameter
  ///    - omega_b: Constant for repulsion parameter
  ///
  template <typename Derived>
  class cubic_eos_base
  {
  public:
    using scalar_type = typename detail::cubic_eos_traits<Derived>::scalar_type;

    static constexpr auto omega_a = detail::cubic_eos_traits<Derived>::omega_a;
    static constexpr auto omega_b = detail::cubic_eos_traits<Derived>::omega_b;

    // Constructors

    cubic_eos_base() = default;

    /// @brief Constructs cubic EoS
    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    cubic_eos_base(const scalar_type &pc, const scalar_type &tc) noexcept
        : pc_{pc},
          tc_{tc},
          ac_{this->critical_attraction_param(pc, tc)},
          bc_{this->critical_repulsion_param(pc, tc)} {}

    // Member functions

    /// @brief Set parameters
    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    void set_params(const scalar_type &pc, const scalar_type &tc) noexcept
    {
      pc_ = pc;
      tc_ = tc;
      ac_ = this->critical_attraction_param(pc, tc);
      bc_ = this->critical_repulsion_param(pc, tc);
    }

    /// @brief Creates isothermal state
    /// @param[in] t Temperature
    isothermal_state<Derived> create_isothermal_state(const scalar_type &t) const noexcept
    {
      const auto tr = this->reduced_temperature(t);
      const auto a = this->attraction_param(tr);
      const auto b = this->repulsion_param();
      return {t, a, b};
    }

    /// @brief Creates isobaric-isothermal state
    /// @param[in] p Pressure
    /// @param[in] t Temperature
    isobaric_isothermal_state<Derived> create_isobaric_isothermal_state(const scalar_type &p, const scalar_type &t) const noexcept
    {
      const auto pr = this->reduced_pressure(p);
      const auto tr = this->reduced_temperature(t);
      const auto ar = this->reduced_attraction_param(pr, tr);
      const auto br = this->reduced_repulsion_param(pr, tr);
      const auto beta = this->derived().beta(tr);
      return {ar, br, beta};
    }

    /// @brief Computes pressure at given temperature and volume
    /// @param[in] t Temperature
    /// @param[in] v Volume
    scalar_type pressure(const scalar_type &t, const scalar_type &v) const noexcept
    {
      const auto tr = this->reduced_temperature(t);
      const auto a = this->attraction_param(tr);
      const auto b = this->repulsion_param();
      return Derived::pressure_impl(t, v, a, b);
    }

    /// @brief Computes Z-factor at given pressure and temperature
    /// @param[in] p Pressure
    /// @param[in] t Temperature
    /// @return A list of Z-factors
    std::vector<scalar_type> zfactor(const scalar_type &p, const scalar_type &t) const noexcept
    {
      return this->create_isobaric_isothermal_state(p, t).zfactor();
    }

  private:
    // Static functions

    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    static scalar_type critical_attraction_param(const scalar_type &pc, const scalar_type &tc) noexcept
    {
      constexpr auto R = gas_constant<scalar_type>();
      return (omega_a * R * R) * tc * tc / pc;
    }

    /// @param[in] pc Critical pressure
    /// @param[in] tc Critical temperature
    static scalar_type critical_repulsion_param(const scalar_type &pc, const scalar_type &tc) noexcept
    {
      constexpr auto R = gas_constant<scalar_type>();
      return (omega_b * R) * tc / pc;
    }

    // Member functions

    /// @brief Computes reduced pressure
    /// @param[in] p Pressure
    scalar_type reduced_pressure(const scalar_type &p) const noexcept
    {
      return p / pc_;
    }

    /// @brief Computes reduced temperature
    /// @param[in] t Temperature
    scalar_type reduced_temperature(const scalar_type &t) const noexcept
    {
      return t / tc_;
    }

    /// @brief Returns attraction parameter at a given temperature.
    /// @param[in] tr Reduced temperature
    scalar_type attraction_param(const scalar_type &tr) const noexcept
    {
      return this->derived().alpha(tr) * ac_;
    }

    /// @brief Returns repulsion parameter.
    scalar_type repulsion_param() const noexcept { return bc_; }

    /// @brief Returns reduced attraction parameter at a given pressure and
    /// temperature.
    /// @param[in] pr Reduced pressure
    /// @param[in] tr Reduced temperature
    scalar_type reduced_attraction_param(const scalar_type &pr, const scalar_type &tr) const noexcept
    {
      return this->derived().alpha(tr) * omega_a * pr / (tr * tr);
    }

    /// @brief Returns reduced repulsion parameter at a given pressure and
    /// temperature.
    /// @param[in] pr Reduced pressure
    /// @param[in] tr Reduced temperature
    scalar_type reduced_repulsion_param(const scalar_type &pr, const scalar_type &tr) const noexcept
    {
      return omega_b * pr / tr;
    }

    // Member functions

    /// @brief Get reference to derived class object
    Derived &derived() noexcept { return static_cast<Derived &>(*this); }

    /// @brief Get const reference to derived class object
    const Derived &derived() const noexcept
    {
      return static_cast<const Derived &>(*this);
    }

    scalar_type pc_; /// Critical pressure
    scalar_type tc_; /// Critical temperature
    scalar_type ac_; /// Critical attraction parameter
    scalar_type bc_; /// Critical repulsion parameter
  };

} // namespace eos