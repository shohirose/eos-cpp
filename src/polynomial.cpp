#include "eos/math/polynomial.hpp"
#include <complex>
#include <algorithm>
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

    template <int N>
    std::vector<double> real_roots(const std::array<double, N> &a) noexcept
    {
        // Array-oriented access of the array of std::complex is guranteed.
        // Please refer to
        // https://en.cppreference.com/w/cpp/numeric/complex
        std::array<std::complex<double>, N - 1> z;
        auto *w = gsl_poly_complex_workspace_alloc(N);
        gsl_poly_complex_solve(a.data(), N, w, reinterpret_cast<double*>(z.data()));
        gsl_poly_complex_workspace_free(w);

        std::vector<double> x;
        x.reserve(N - 1);
        constexpr double tolerance = 1e-8;
        for (int i = 0; i < N - 1; ++i)
        {
            if (std::fabs(z[i].imag()) < tolerance)
            {
                x.push_back(z[i].real());
            }
        }
        std::sort(x.begin(), x.end());
        return x;
    }

    template std::vector<double> real_roots<4>(const std::array<double, 4> &a) noexcept;
    template std::vector<double> real_roots<5>(const std::array<double, 5> &a) noexcept;
    template std::vector<double> real_roots<6>(const std::array<double, 6> &a) noexcept;

} // namespace eos
