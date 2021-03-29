#pragma once

#include <algorithm> // std::max_element
#include <cassert>   // assert
#include <cmath>     // std::pow
#include <iostream>  // std::cerr
#include <type_traits>
#include <utility> // std::pair, std::forward

namespace eos
{

  /// @brief Estimates vapor pressure of a pure component by using Wilson
  /// equation.
  template <typename T>
  T estimate_vapor_pressure(const T &t, const T &pc, const T &tc, const T &omega) noexcept
  {
    assert(t <= tc);
    using std::pow;
    return pc * pow(10, 7.0 / 3.0 * (1 + omega) * (1 - tc / t));
  }

  namespace detail
  {

    template <typename Eos>
    struct cubic_eos_traits;

  } // namespace detail

  enum class FlashErrorCode
  {
    Success,               /// Calculation succeeded
    MaxIterReached,        /// Maximum iteration reached
    MultipleRootsNotFound, /// Multiple roots not found in z-factor
  };

  /// @brief vapor_liquid_flash calculation class
  template <typename Eos>
  class vapor_liquid_flash
  {
  public:
    using scalar_type = typename detail::cubic_eos_traits<Eos>::scalar_type;

    vapor_liquid_flash() = default;
    vapor_liquid_flash(const vapor_liquid_flash &) = default;
    vapor_liquid_flash(vapor_liquid_flash &&) = default;

    /// @brief Constructs flash object
    /// @param[in] eos EoS
    ///
    /// The default values of tolerance and maxixum iteration are 1e-6 and 100,
    /// respectively.
    vapor_liquid_flash(const Eos &eos) : eos_{eos}, tol_{1e-6}, maxiter_{100} {}

    /// @brief Constructs flash object
    /// @param[in] eos EoS
    /// @param[in] tol Tolerance for vapor_liquid_flash calculation convergence
    /// @param[in] maxiter Maximum iteration
    vapor_liquid_flash(const Eos &eos, const scalar_type &tol, int maxiter)
        : eos_{eos}, tol_{tol}, maxiter_{maxiter} {}

    vapor_liquid_flash &operator=(const vapor_liquid_flash &) = default;
    vapor_liquid_flash &operator=(vapor_liquid_flash &&) = default;

    enum class error_code
    {
      success,
      not_converged,
      multiple_roots_not_found,
    };

    struct result
    {
      scalar_type rsd;  /// Relative residual
      int iter;         /// Iteration count
      error_code error; /// Error code
    };

    /// @brief Computes vapor pressure
    /// @param[in] p_init Initial pressure
    /// @param[in] t Temperature
    /// @return A pair of vapor pressure and iteration report
    std::pair<scalar_type, result> vapor_pressure(const scalar_type &p_init, const scalar_type &t) const noexcept
    {
      auto p = p_init;
      scalar_type eps = 1;
      int iter = 0;

      using std::fabs;

      while (eps > tol_ && iter < maxiter_)
      {
        const auto state = eos_.create_isobaric_isothermal_state(p, t);
        const auto z = state.zfactor();

        if (z.size() < 2)
        {
          return {0, {eps, iter, error_code::multiple_roots_not_found}};
        }

        const auto zv = *std::max_element(z.begin(), z.end());
        const auto zl = *std::min_element(z.begin(), z.end());
        const auto phiv = state.fugacity_coeff(zv);
        const auto phil = state.fugacity_coeff(zl);

        eps = fabs(1 - phil / phiv);

        // Update vapor pressure by successive substitution
        p *= phil / phiv;

        ++iter;
      }

      if (iter >= maxiter_)
      {
        return {0, {eps, iter, error_code::not_converged}};
      }

      return {p, {eps, iter, error_code::success}};
    }

    scalar_type tolerance() const noexcept { return tol_; }
    int maxIter() const noexcept { return maxiter_; }

    void set_params(const scalar_type &tol, int maxiter)
    {
      tol_ = tol;
      maxiter_ = maxiter;
    }

  private:
    Eos eos_;
    scalar_type tol_;
    int maxiter_;
  };

  template <typename Eos>
  inline vapor_liquid_flash<Eos> make_vapor_liquid_flash(const Eos &eos)
  {
    return {eos};
  }

} // namespace eos