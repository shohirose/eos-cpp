#include "eos/math/polynomial_solver.hpp"

#include <cassert>
#include <complex>

#include "gsl_workspace_wrapper.hpp"

namespace eos {

class polynomial_solver::impl {
 public:
  impl(std::size_t n)
      : workspace_{std::make_unique<gsl_workspace_wrapper>(n)},
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
  std::unique_ptr<gsl_workspace_wrapper> workspace_;
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
