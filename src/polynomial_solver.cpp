#include "eos/math/polynomial_solver.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <complex>
#include <cassert>

namespace eos {

class workspace_wrapper {
 public:
  workspace_wrapper(std::size_t n) : n_{n}, w_{nullptr} {
    if (n_ < 2) {
      throw std::invalid_argument(
          "Error: n must be equal to or larger than 2!");
    }
    w_ = gsl_poly_complex_workspace_alloc(n);
    if (!w_) {
      throw std::runtime_error("Error: Workspace allocation failed!");
    }
  }

  workspace_wrapper(workspace_wrapper &&other) : n_{other.n_}, w_{other.w_} {
    other.n_ = 0;
    other.w_ = nullptr;
  }

  ~workspace_wrapper() {
    if (w_) {
      gsl_poly_complex_workspace_free(w_);
    }
  }

  workspace_wrapper &operator=(workspace_wrapper &&other) {
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
      throw std::runtime_error("Error: Workspace allocation failed!");
    }
  }

 private:
  std::size_t n_;
  gsl_poly_complex_workspace *w_;
};

class polynomial_solver::impl {
 public:
  impl(std::size_t n)
      : workspace_{std::make_unique<workspace_wrapper>(n)},
        roots_(n - 1),
        tol_{1e-8} {}

  impl(const impl &) = delete;
  impl(impl &&other)
      : workspace_{std::move(other.workspace_)},
        roots_{std::move(other.roots_)},
        tol_{other.tol_} {}

  impl &operator=(const impl &) = delete;
  impl &operator=(impl &&other) {
    workspace_ = std::move(other.workspace_);
    roots_ = std::move(other.roots_);
    tol_ = other.tol_;
    return *this;
  };

  ~impl() {}

  std::vector<double> solve(gsl::span<const double> coeffs) {
    // Array-oriented access of the array of std::complex is guranteed.
    // Please refer to
    // https://en.cppreference.com/w/cpp/numeric/complex
    workspace_->solve(coeffs,
                      gsl::make_span(reinterpret_cast<double *>(roots_.data()),
                                     roots_.size() * 2));

    std::vector<double> x;
    x.reserve(roots_.size());

    for (const auto &z : roots_) {
      if (std::fabs(z.imag()) < tol_) {
        x.push_back(z.real());
      }
    }

    std::sort(x.begin(), x.end());
    return x;
  }

  // If n < 2, std::invalid_argument error will be thrown.
  void reset(std::size_t n) {
    assert(workspace_);
    workspace_->reset(n);
    roots_.resize(n - 1);
  }

  void set_tolerance(double tol) { tol_ = tol; }

  std::unique_ptr<impl> clone() const {
    return std::make_unique<impl>(roots_.size() + 1);
  }

 private:
  std::unique_ptr<workspace_wrapper> workspace_;
  std::vector<std::complex<double>> roots_;
  double tol_;
};

polynomial_solver::polynomial_solver() : pimpl_{} {}

polynomial_solver::polynomial_solver(size_t num_coeffs)
    : pimpl_{std::make_unique<impl>(num_coeffs)} {}

polynomial_solver::polynomial_solver(const polynomial_solver &other)
    : pimpl_{} {
  if (other.pimpl_) {
    pimpl_ = other.pimpl_->clone();
  }
}

polynomial_solver::polynomial_solver(polynomial_solver &&other)
    : pimpl_{std::move(other.pimpl_)} {}

polynomial_solver::~polynomial_solver() {}

polynomial_solver &polynomial_solver::operator=(
    const polynomial_solver &other) {
  if (other.pimpl_) {
    pimpl_ = other.pimpl_->clone();
  }
  return *this;
}

polynomial_solver &polynomial_solver::operator=(polynomial_solver &&other) {
  pimpl_ = std::move(other.pimpl_);
  return *this;
}

void polynomial_solver::reset(size_t num_coeffs) {
  if (pimpl_) {
    pimpl_->reset(num_coeffs);
  } else {
    pimpl_ = std::make_unique<impl>(num_coeffs);
  }
}

std::vector<double> polynomial_solver::solve(gsl::span<const double> a) {
  this->reset(a.size());
  return pimpl_->solve(a);
}

}  // namespace eos
