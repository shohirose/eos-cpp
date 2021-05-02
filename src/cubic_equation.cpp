#include "eos/math/cubic_equation.hpp"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>

#include <algorithm>

namespace eos {

std::vector<double> cubic_equation::real_roots() const {
  std::vector<double> x(3);
  const auto num_roots =
      gsl_poly_solve_cubic(this->a, this->b, this->c, &x[0], &x[1], &x[2]);
  x.resize(num_roots);
  return x;
}

std::array<std::complex<double>, 3> cubic_equation::complex_roots() const {
  std::array<gsl_complex, 3> z1;
  gsl_poly_complex_solve_cubic(this->a, this->b, this->c,  //
                               &z1[0], &z1[1], &z1[2]);

  std::array<std::complex<double>, 3> z2;
  std::transform(begin(z1), end(z1), begin(z2), [](const auto& z) {
    return std::complex<double>{GSL_REAL(z), GSL_IMAG(z)};
  });
  return z2;
}

}  // namespace eos
