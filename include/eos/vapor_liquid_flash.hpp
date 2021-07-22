#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <utility>

namespace eos {

/// @brief Estimates vapor pressure of a pure component by using Wilson
/// equation.
inline double estimateVaporPressure(double t, double pc, double tc,
                                      double omega) noexcept {
  assert(t <= tc);
  return pc * std::pow(10, 7.0 / 3.0 * (1 + omega) * (1 - tc / t));
}

enum class FlashError {
  success,
  not_converged,
  multiple_roots_not_found,
};

struct flash_iteration_result {
  double rsd;                   /// Relative residual
  int iter;                     /// Iteration count
  FlashError error;  /// Error code
};

/// @brief VaporLiquidFlash calculation class
template <typename CubicEos>
class VaporLiquidFlash {
 public:
  VaporLiquidFlash() = default;
  VaporLiquidFlash(const VaporLiquidFlash &) = default;
  VaporLiquidFlash(VaporLiquidFlash &&) = default;

  VaporLiquidFlash &operator=(const VaporLiquidFlash &) = default;
  VaporLiquidFlash &operator=(VaporLiquidFlash &&) = default;

  /// @brief Constructs flash object
  /// @param[in] eos EoS
  ///
  /// The default values of tolerance and maxixum iteration are 1e-6 and 100,
  /// respectively.
  VaporLiquidFlash(const CubicEos &eos)
      : eos_{eos}, tol_{1e-6}, maxiter_{100} {}

  /// @brief Constructs flash object
  /// @param[in] eos EoS
  /// @param[in] tol Tolerance for VaporLiquidFlash calculation convergence
  /// @param[in] maxiter Maximum iteration
  VaporLiquidFlash(const CubicEos &eos, double tol, int maxiter)
      : eos_{eos}, tol_{tol}, maxiter_{maxiter} {}

  /// @brief Computes vapor pressure
  /// @param[in] p_init Initial pressure
  /// @param[in] t Temperature
  /// @return A pair of vapor pressure and iteration report
  std::pair<double, flash_iteration_result> vaporPressure(
      double p_init, double t) const noexcept {
    auto p = p_init;
    double eps = 1.0;
    int iter = 0;

    while (eps > tol_ && iter < maxiter_) {
      const auto state = eos_.createIsobaricIsothermalState(p, t);
      const auto z = state.zfactor();

      if (z.size() < 2) {
        return {0.0,
                {eps, iter, FlashError::multiple_roots_not_found}};
      }

      const auto zv = *std::max_element(z.begin(), z.end());
      const auto zl = *std::min_element(z.begin(), z.end());
      const auto phiv = state.fugacityCoeff(zv);
      const auto phil = state.fugacityCoeff(zl);

      eps = std::fabs(1.0 - phil / phiv);

      // Update vapor pressure by successive substitution
      p *= phil / phiv;

      ++iter;
    }

    if (iter >= maxiter_) {
      return {0.0, {eps, iter, FlashError::not_converged}};
    } else {
      return {p, {eps, iter, FlashError::success}};
    }
  }

  double tolerance() const noexcept { return tol_; }
  int maxIter() const noexcept { return maxiter_; }

  void setParams(double tol, int maxiter) {
    tol_ = tol;
    maxiter_ = maxiter;
  }

 private:
  CubicEos eos_;
  double tol_;
  int maxiter_;
};

template <typename CubicEos>
inline VaporLiquidFlash<CubicEos> makeVaporLiquidFlash(
    const CubicEos &eos) {
  return {eos};
}

}  // namespace eos