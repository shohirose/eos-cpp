#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <type_traits>

#include "eos/cubic_eos/cubic_eos_base.hpp"

namespace eos
{

  /// @brief Estimates vapor pressure of a pure component by using Wilson
  /// equation.
  inline double estimate_vapor_pressure(double t, double pc, double tc, double omega) noexcept
  {
    assert(t <= tc);
    return pc * std::pow(10, 7.0 / 3.0 * (1 + omega) * (1 - tc / t));
  }

  enum class flash_iteration_error
  {
    success,
    not_converged,
    multiple_roots_not_found,
  };

  struct flash_iteration_result
  {
    double rsd;                  /// Relative residual
    int iter;                    /// Iteration count
    flash_iteration_error error; /// Error code
  };

  /// @brief vapor_liquid_flash calculation class
  template <typename CubicEos>
  class vapor_liquid_flash
  {
  public:
    static_assert(std::is_base_of_v<cubic_eos_base<CubicEos>, CubicEos>, "CubicEos must be derived from eos::cubic_eos_base");

    vapor_liquid_flash() = default;
    vapor_liquid_flash(const vapor_liquid_flash &) = default;
    vapor_liquid_flash(vapor_liquid_flash &&) = default;

    vapor_liquid_flash &operator=(const vapor_liquid_flash &) = default;
    vapor_liquid_flash &operator=(vapor_liquid_flash &&) = default;

    /// @brief Constructs flash object
    /// @param[in] eos EoS
    ///
    /// The default values of tolerance and maxixum iteration are 1e-6 and 100,
    /// respectively.
    vapor_liquid_flash(const CubicEos &eos) : eos_{eos}, tol_{1e-6}, maxiter_{100} {}

    /// @brief Constructs flash object
    /// @param[in] eos EoS
    /// @param[in] tol Tolerance for vapor_liquid_flash calculation convergence
    /// @param[in] maxiter Maximum iteration
    vapor_liquid_flash(const CubicEos &eos, double tol, int maxiter)
        : eos_{eos}, tol_{tol}, maxiter_{maxiter} {}

    vapor_liquid_flash &operator=(const vapor_liquid_flash &) = default;
    vapor_liquid_flash &operator=(vapor_liquid_flash &&) = default;

    /// @brief Computes vapor pressure
    /// @param[in] p_init Initial pressure
    /// @param[in] t Temperature
    /// @return A pair of vapor pressure and iteration report
    std::pair<double, flash_iteration_result> vapor_pressure(double p_init, double t) const noexcept
    {
      auto p = p_init;
      double eps = 1.0;
      int iter = 0;

      while (eps > tol_ && iter < maxiter_)
      {
        const auto state = eos_.create_isobaric_isothermal_state(p, t);
        const auto z = state.zfactor();

        if (z.size() < 2)
        {
          return {0.0, {eps, iter, flash_iteration_error::multiple_roots_not_found}};
        }

        const auto zv = *std::max_element(z.begin(), z.end());
        const auto zl = *std::min_element(z.begin(), z.end());
        const auto phiv = state.fugacity_coeff(zv);
        const auto phil = state.fugacity_coeff(zl);

        eps = std::fabs(1.0 - phil / phiv);

        // Update vapor pressure by successive substitution
        p *= phil / phiv;

        ++iter;
      }

      if (iter >= maxiter_)
      {
        return {0.0, {eps, iter, flash_iteration_error::not_converged}};
      }
      else
      {
        return {p, {eps, iter, flash_iteration_error::success}};
      }
    }

    double tolerance() const noexcept { return tol_; }
    int maxIter() const noexcept { return maxiter_; }

    void set_params(double tol, int maxiter)
    {
      tol_ = tol;
      maxiter_ = maxiter;
    }

  private:
    CubicEos eos_;
    double tol_;
    int maxiter_;
  };

  template <typename CubicEos>
  inline vapor_liquid_flash<CubicEos> make_vapor_liquid_flash(const CubicEos &eos)
  {
    return {eos};
  }

} // namespace eos