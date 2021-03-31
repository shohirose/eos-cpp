#include "eos/viscosity/lucas_method.hpp"
#include <cassert>

namespace eos
{

    void lucas_method::set_params(double pc, double tc, double zc, double mw, double dm, double q) noexcept
    {
        pc_ = pc;
        tc_ = tc;
        zc_ = zc;
        mw_ = mw;
        dm_ = dm;
        q_ = q;
        dmr_ = reduced_dipole_moment(dm, tc, pc);
        xi_ = inverse_viscosity(pc, tc, mw);
    }

    double lucas_method::viscosity_at_low_pressure(double t) const noexcept
    {
        const auto tr = this->reduced_temperature(t);
        const auto fp = this->polarity_factor_at_low_pressure(tr);
        const auto fq = this->quantum_factor_at_low_pressure(tr);
        const auto z1 = this->reduced_viscosity_at_low_pressure(tr, fp, fq);
        return z1 / xi_;
    }

    double lucas_method::viscosity_at_high_pressure(double p, double t) const noexcept
    {
        const auto pr = this->reduced_pressure(p);
        const auto tr = this->reduced_temperature(t);
        const auto fp0 = this->polarity_factor_at_low_pressure(tr);
        const auto fq0 = this->quantum_factor_at_low_pressure(tr);

        const auto z1 = this->reduced_viscosity_at_low_pressure(tr, fp0, fq0);
        const auto z2 = this->reduced_viscosity_at_high_pressure(z1, pr, tr);
        const auto fp = this->polarity_factor_at_high_pressure(fp0, z1, z2);
        const auto fq = this->quantum_factor_at_high_pressure(fq0, z1, z2);

        return z2 * fp * fq / xi_;
    }

    double lucas_method::reduced_viscosity_at_high_pressure(double z1, double pr, double tr) noexcept
    {
        using std::exp;
        using std::pow;

        if (tr <= static_cast<double>(1.0))
        {
            const auto alpha = 3.262 + 14.98 * pow(pr, 5.508);
            const auto beta = 1.390 + 5.746 * pr;
            return 0.600 + 0.760 * pow(pr, alpha) + (6.990 * pow(pr, beta) - 0.6) * (1.0 - tr);
        }
        else
        {
            const auto a = 1.245e-3 / tr * exp(5.1726 * pow(tr, -0.3286));
            const auto b = a * (1.6553 * tr - 1.2723);
            const auto c = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
            const auto d = 1.7369 / tr * exp(2.2310 * pow(tr, -7.6351));
            const auto f = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
            const auto z2 = 1.0 + a * pow(pr, 1.3088) / (b * pow(pr, f) + 1.0 / (1.0 + c * pow(pr, d)));
            return z2 * z1;
        }
    }

    double lucas_method::polarity_factor_at_low_pressure(double tr) const noexcept
    {
        using std::fabs;
        using std::pow;
        assert(dmr_ >= 0);
        if (dmr_ < static_cast<double>(0.022))
        {
            return 1.0;
        }
        else if (dmr_ < static_cast<double>(0.075))
        {
            return 1.0 + 30.55 * pow(0.292 - zc_, 1.72);
        }
        else
        {
            return 1.0 + 30.55 * pow(0.292 - zc_, 1.72) * fabs(0.96 + 0.1 * (tr - 0.7));
        }
    }

    double lucas_method::quantum_factor_at_low_pressure(double tr) const noexcept
    {
        using std::copysign;
        using std::fabs;
        using std::pow;
        constexpr auto tolerance = static_cast<double>(1e-10);
        assert(q_ >= 0);
        if (fabs(q_) < tolerance)
        {
            return 1.0;
        }
        else if (fabs(tr - 12) < tolerance)
        {
            return 1.22 * pow(q_, 0.15);
        }
        else
        {
            const auto tmp = tr - 12;
            return 1.22 * pow(q_, 0.15) * (1.0 + copysign(0.00385 * pow(tmp * tmp, 1.0 / mw_), tmp));
        }
    }

} // namespace eos
