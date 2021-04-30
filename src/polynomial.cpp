#include "eos/math/polynomial.hpp"

#include <gsl/gsl_poly.h>

#include <algorithm>
#include <cassert>
#include <complex>

namespace eos {

std::vector<double> real_roots(double a, double b, double c) noexcept {
  std::vector<double> x(3);
  const auto num_roots = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
  x.resize(num_roots);
  return x;
}

std::vector<double> real_roots(const std::vector<double> &a) noexcept {
  assert(a.size() > 0);
  std::vector<std::complex<double>> z(a.size() - 1);
  auto *w = gsl_poly_complex_workspace_alloc(a.size());
  // Array-oriented access of the array of std::complex is guranteed.
  // Please refer to
  // https://en.cppreference.com/w/cpp/numeric/complex
  gsl_poly_complex_solve(a.data(), a.size(), w,
                         reinterpret_cast<double *>(z.data()));
  gsl_poly_complex_workspace_free(w);

  std::vector<double> x;
  x.reserve(a.size() - 1);
  constexpr double tolerance = 1e-8;
  for (size_t i = 0; i < a.size() - 1; ++i) {
    if (std::fabs(z[i].imag()) < tolerance) {
      x.push_back(z[i].real());
    }
  }
  std::sort(x.begin(), x.end());
  return x;
}

}  // namespace eos
