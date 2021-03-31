#include "eos/math/polynomial.hpp"
#include <gsl/gsl_poly.h>

namespace eos
{

    std::vector<double> real_roots(const std::array<double, 3> &a) noexcept
    {
        std::vector<double> x(3);
        const auto num_roots = gsl_poly_solve_cubic(a[0], a[1], a[2], &x[0], &x[1], &x[2]);
        x.resize(num_roots);
        return x;
    }

} // namespace eos
