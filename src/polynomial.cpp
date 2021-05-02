#include "eos/math/polynomial.hpp"

#include <algorithm>
#include <complex>

#include "gsl_workspace_wrapper.hpp"

namespace eos {

std::vector<double> real_roots(double a, double b, double c) noexcept {
  std::vector<double> x(3);
  const auto num_roots = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
  x.resize(num_roots);
  return x;
}

std::vector<double> real_roots(gsl::span<const double> a) {
  gsl_workspace_wrapper w(a.size());
  std::vector<std::complex<double>> z(a.size() - 1);

  // Array-oriented access of the array of std::complex is guranteed.
  // Please refer to
  // https://en.cppreference.com/w/cpp/numeric/complex
  w.solve(a,
          gsl::make_span(reinterpret_cast<double *>(z.data()), 2 * z.size()));

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
