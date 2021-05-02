#pragma once

#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

#include <cstdint>
#include <gsl/gsl>
#include <stdexcept>

namespace eos {

class gsl_workspace_wrapper {
 public:
  gsl_workspace_wrapper(std::size_t n) : n_{n}, w_{nullptr} {
    if (n_ < 2) {
      throw std::invalid_argument(
          "Error: n must be equal to or larger than 2!");
    }
    w_ = gsl_poly_complex_workspace_alloc(n);
    if (!w_) {
      throw std::runtime_error(
          "Error: gsl_poly_complex_workspace_alloc failed!");
    }
  }

  gsl_workspace_wrapper(gsl_workspace_wrapper &&other)
      : n_{other.n_}, w_{other.w_} {
    other.n_ = 0;
    other.w_ = nullptr;
  }

  ~gsl_workspace_wrapper() {
    if (w_) {
      gsl_poly_complex_workspace_free(w_);
    }
  }

  gsl_workspace_wrapper &operator=(gsl_workspace_wrapper &&other) {
    n_ = other.n_;
    w_ = other.w_;
    other.n_ = 0;
    other.w_ = nullptr;
  }

  /// @brief Solve a polynomial
  /// @param[in] a Coefficients
  /// @param[out] z Complex roots
  void solve(gsl::span<const double> a, gsl::span<double> z) {
    if (a.size() != n_) {
      throw std::invalid_argument(
          "Error: the number of coefficients is incorrect!");
    }
    if (z.size() != 2 * (n_ - 1)) {
      throw std::invalid_argument(
          "Error: the number of complex roots is incorrect!");
    }

    auto *handler = gsl_set_error_handler_off();
    const auto status = gsl_poly_complex_solve(a.data(), n_, w_, z.data());
    gsl_set_error_handler(handler);

    if (status == GSL_EFAILED) {
      throw std::runtime_error("Error: gsl_poly_complex_solve failed!");
    }
  }

  void reset(std::size_t n) {
    if (n_ < 2) {
      throw std::invalid_argument(
          "Error: n must be equal to or larger than 2!");
    }

    if (n_ == n) {
      return;
    }

    n_ = n;
    if (w_) {
      gsl_poly_complex_workspace_free(w_);
    }
    w_ = gsl_poly_complex_workspace_alloc(n);
    if (!w_) {
      throw std::runtime_error(
          "Error: gsl_poly_complex_workspace_alloc failed!");
    }
  }

 private:
  std::size_t n_;
  gsl_poly_complex_workspace *w_;
};

}  // namespace eos
