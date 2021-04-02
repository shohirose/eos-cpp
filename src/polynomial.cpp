#include "eos/math/polynomial.hpp"
#include <gsl/gsl_poly.h>

namespace eos
{

    std::vector<double> real_roots(double a, double b, double c) noexcept
    {
        std::vector<double> x(3);
        const auto num_roots = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
        x.resize(num_roots);
        return x;
    }

} // namespace eos
