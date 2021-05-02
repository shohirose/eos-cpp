#include "eos/math/polynomial.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

#include <algorithm>
#include <cassert>
#include <complex>
#include <stdexcept>

namespace eos {

std::vector<double> real_roots(double a, double b, double c) noexcept {
  std::vector<double> x(3);
  const auto num_roots = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
  x.resize(num_roots);
  return x;
}

std::vector<double> real_roots(gsl::span<const double> a) {
  assert(a.size() > 0);
  std::vector<std::complex<double>> z(a.size() - 1);

  auto *w = gsl_poly_complex_workspace_alloc(a.size());
  if (!w) {
    throw std::runtime_error("Error: gsl_poly_complex_workspace_alloc failed!");
  }

  auto *handler = gsl_set_error_handler_off();
  // Array-oriented access of the array of std::complex is guranteed.
  // Please refer to
  // https://en.cppreference.com/w/cpp/numeric/complex
  const auto status = gsl_poly_complex_solve(
      a.data(), a.size(), w, reinterpret_cast<double *>(z.data()));
  gsl_poly_complex_workspace_free(w);
  gsl_set_error_handler(handler);

  if (status == GSL_EFAILED) {
    throw std::runtime_error("Error: gsl_poly_complex_solve failed!");
  }

  std::vector<double> x;
  x.reserve(z.size());
  constexpr double tol = 1e-8;
  for (const auto &zi : z) {
    if (std::fabs(zi.imag()) < tol) {
      x.push_back(zi.real());
    }
  }
  std::sort(x.begin(), x.end());
  return x;
}

}  // namespace eos
