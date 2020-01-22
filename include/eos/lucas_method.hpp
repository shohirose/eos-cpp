// MIT License
//
// Copyright (c) 2019 Sho Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <cmath>        // std::log, std::exp, std::fabs, std::pow
#include <type_traits>  // std::conditional

#include <Eigen/Core>  // Eigen::Matrix, Eigen::Index, Eigen::Dynamic

#include "eos/constants.hpp"  // eos::gas_constant, eos::dynamic_extent

namespace eos {

/// @brief Base class for Lucas method
/// @tparam T Value type
template <typename T>
struct LucasMethodBase {
  // Static functions

  /// @brief Computes reduced dipole moment
  /// @param[in] dm Dipole moment
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical Temperature
  /// @return Reduced dipole moment
  static T reduced_dipole_moment(const T& dm, const T& pc,
                                 const T& tc) noexcept {
    const auto tmp = dm / tc;
    return 52.46e-5 * tmp * tmp * pc;
  }

  /// @brief Computes polarity factor at low pressure
  /// @param[in] dmr Reduced dipole moment
  /// @param[in] zc Critical z-factor
  /// @param[in] tr Reduced temperature
  static T low_pressure_polarity_factor(const T& dmr, const T& zc,
                                        const T& tr) noexcept {
    using std::fabs;
    using std::pow;
    assert(dmr >= 0.0);
    if (dmr >= 0.0 && dmr < 0.022) {
      return 1.0;
    } else if (dmr < 0.075) {
      return 1.0 + 30.55 * pow(0.292 - zc, 1.72);
    } else {
      return 1.0 +
             30.55 * pow(0.292 - zc, 1.72) * fabs(0.96 + 0.1 * (tr - 0.7));
    }
  }

  /// @brief Computes quantum factor at low pressure
  /// @param[in] q Quantum parameter
  /// @param[in] tr Reduced pressure
  /// @param[in] mw Molecular weight
  ///
  /// Quantum factor is required only for quantum gases: H2, He, and D2.
  /// q = 1.38 (He), q = 0.76 (H2), and q = 0.52 (D2).
  static T low_pressure_quantum_factor(const T& q, const T& tr,
                                       const T& mw) noexcept {
    using std::copysign;
    using std::pow;
    if (q == 0.0) {
      return 1.0;
    } else {
      const auto tmp = tr - 12.0;
      if (tmp == 0.0) {
        return 1.22 * pow(q, 0.15);
      } else {
        return 1.22 * pow(q, 0.15) *
               (1.0 + copysign(0.00385 * pow(tmp, 2.0 / mw), tmp));
      }
    }
  }

  /// @param[in] fp0 Polarity factor at low pressure
  /// @param[in] z1 Reduced viscosity at low pressure
  /// @param[in] z2 Reduced viscosity at high pressure
  static T polarity_factor(const T& fp0, const T& z1, const T& z2) noexcept {
    const auto y = z2 / z1;
    return (1.0 + (fp0 - 1.0) / (y * y * y)) / fp0;
  }

  /// @param[in] fq0 Quantum factor at low pressure
  /// @param[in] z1 Reduced viscosity at low pressure
  /// @param[in] z2 Reduced viscosity at high pressure
  static T quantum_factor(const T& fq0, const T& z1, const T& z2) noexcept {
    const auto y = z2 / z1;
    using std::log;
    const auto tmp = log(y);
    const auto tmp2 = tmp * tmp;
    return (1.0 + (fq0 - 1.0) * (1.0 / y - 0.007 * tmp2 * tmp2)) / fq0;
  }

  /// @param[in] tr Reduced temperature
  /// @param[in] fp Polarity factor
  /// @param[in] fq Quantum factor
  static T low_pressure_reduced_viscosity(const T& tr, const T& fp,
                                          const T& fq) noexcept {
    using std::exp;
    using std::pow;
    const auto z1 = 0.807 * pow(tr, 0.618) - 0.357 * exp(-0.449 * tr) +
                    0.340 * exp(-4.058 * tr) + 0.018;
    return z1 * fp * fq;
  }

  /// @brief Computes reduced viscosity at high pressure
  /// @param[in] z1 Reduced viscosity at low pressure
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static T high_pressure_reduced_viscosity(const T& z1, const T& pr,
                                           const T& tr) noexcept {
    return (tr <= 1.0) ? low_temperature_reduced_viscosity(pr, tr)
                       : high_temperature_reduced_viscosity(z1, pr, tr);
  }

  /// @brief Computes reduced viscosity at low temperature
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static T low_temperature_reduced_viscosity(const T& pr,
                                             const T& tr) noexcept {
    using std::pow;
    const auto alpha = 3.262 + 14.98 * pow(pr, 5.508);
    const auto beta = 1.390 + 5.746 * pr;
    return 0.600 + 0.760 * pow(pr, alpha) +
           (6.990 * pow(pr, beta) - 0.6) * (1.0 - tr);
  }

  /// @brief Computes reduced viscosity at high temperature
  /// @param[in] z1 Reduced viscosity at low pressure
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static T high_temperature_reduced_viscosity(const T& z1, const T& pr,
                                              const T& tr) noexcept {
    using std::exp;
    using std::pow;
    const auto a = 1.245e-3 / tr * exp(5.1726 * pow(tr, -0.3286));
    const auto b = a * (1.6553 * tr - 1.2723);
    const auto c = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
    const auto d = 1.7369 / tr * exp(2.2310 * pow(tr, -7.6351));
    const auto f = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
    const auto z2 = 1.0 + a * pow(pr, 1.3088) /
                              (b * pow(pr, f) + 1.0 / (1.0 + c * pow(pr, d)));
    return z2 * z1;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] mw Molecular weight
  static T inverse_reduced_viscosity(const T& pc, const T& tc,
                                     const T& mw) noexcept {
    using std::pow;
    const auto pc_bar = pc * 1e-5;
    const auto pc2 = pc_bar * pc_bar;
    return 0.176 * pow(tc / (mw * mw * mw * pc2 * pc2), 1.0 / 6.0);
  }
};

/// @brief Class to computes gas viscosity at high pressure by using Lucas
/// method.
/// @tparam T Value type
/// @tparam N Number of components
///
/// Please refer to Poling et al. (2001) "The Properties of Gases and
/// Liquids", fifth edition. McGRAW-HILL.
template <typename T, std::size_t N>
class LucasMethod : public LucasMethodBase<T> {
 public:
  using Base = LucasMethodBase<T>;
  using Vector = typename std::conditional<N == dynamic_extent,
                                           Eigen::Matrix<T, Eigen::Dynamic, 1>,
                                           Eigen::Matrix<T, N, 1>>::type;

  // Static functions

  /// @brief Computes reduced dipole moment
  /// @param[in] dm Dipole moment
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical Temperature
  /// @return Reduced dipole moment
  template <typename Derived1, typename Derived2, typename Derived3>
  static Vector reduced_dipole_moment(
      const Eigen::MatrixBase<Derived1>& dm,
      const Eigen::MatrixBase<Derived2>& pc,
      const Eigen::MatrixBase<Derived3>& tc) noexcept {
    return 52.46e-5 *
           ((dm.array() / tc.array()).square() * pc.array()).matrix();
  }

  /// @brief Computes polarity factor at low pressure
  /// @param[in] dmr Reduced dipole moment
  /// @param[in] zc Critical z-factor
  /// @param[in] tr Reduced temperature
  template <typename Derived1, typename Derived2, typename Derived3>
  static Vector low_pressure_polarity_factor(
      const Eigen::MatrixBase<Derived1>& dmr,
      const Eigen::MatrixBase<Derived2>& zc,
      const Eigen::MatrixBase<Derived3>& tr) noexcept {
    Vector fp0;
    if constexpr (N == dynamic_extent) {
      fp0.resize(dmr.size());
    }
    for (Eigen::Index i = 0; i < dmr.size(); ++i) {
      fp0[i] = Base::low_pressure_polarity_factor(dmr[i], zc[i], tr[i]);
    }
    return fp0;
  }

  /// @brief Computes quantum factor at low pressure
  /// @param[in] q Quantum parameter
  /// @param[in] tr Reduced pressure
  /// @param[in] mw Molecular weight
  ///
  /// Quantum factor is required only for quantum gases, H2, He, D2.
  /// q = 1.38 (He), q = 0.76 (H2), q = 0.52 (D2).
  /// Please refer to Poling et al. (2001) "The Properties of Gases and
  /// Liquids", fifth edition. McGRAW-HILL.
  template <typename Derived1, typename Derived2, typename Derived3>
  static Vector low_pressure_quantum_factor(
      const Eigen::MatrixBase<Derived1>& q,
      const Eigen::MatrixBase<Derived2>& tr,
      const Eigen::MatrixBase<Derived3>& mw) noexcept {
    Vector fq0;
    if constexpr (N == dynamic_extent) {
      fq0.resize(q.size());
    }
    for (Eigen::Index i = 0; i < q.size(); ++i) {
      fq0[i] = Base::low_pressure_quantum_factor(q[i], tr[i], mw[i]);
    }
    return fq0;
  }

  /// @brief Constructs object
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temprature
  /// @param[in] vc Critical volume
  /// @param[in] zc Critical z-factor
  /// @param[in] mw Molecular weight
  /// @param[in] dm Dipole moment
  /// @param[in] q Quantum parameter for quantum gases
  LucasMethod(const Eigen::Ref<const Vector>& pc,
              const Eigen::Ref<const Vector>& tc,
              const Eigen::Ref<const Vector>& vc,
              const Eigen::Ref<const Vector>& zc,
              const Eigen::Ref<const Vector>& mw,
              const Eigen::Ref<const Vector>& dm,
              const Eigen::Ref<const Vector>& q)
      : pc_{pc}, tc_{tc}, vc_{vc}, zc_{zc}, mw_{mw}, dm_{dm}, q_{q} {}

  /// @brief Computes gas viscosity
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @param[in] x Composition
  T viscosity(const T& p, const T& t, const Eigen::Ref<const Vector>& x) const
      noexcept {
    const Vector tr = (t / tc_.array()).matrix();
    const auto dmr = reduced_dipole_moment(dm_, pc_, tc_);
    const auto fp0 = low_pressure_polarity_factor(dmr, zc_, tr);
    const auto fq0 = low_pressure_quantum_factor(q_, tr, mw_);

    const auto tc_m = x.dot(tc_);
    const auto mw_m = x.dot(mw_);
    const auto fp0_m = x.dot(fp0);
    const auto fq0_m = x.dot(fq0);
    const auto zc_m = x.dot(zc_);
    const auto vc_m = x.dot(vc_);
    const auto pc_m = gas_constant * tc_m * zc_m / vc_m;
    const auto pr_m = p / pc_m;
    const auto tr_m = t / tc_m;

    const auto z1 = Base::low_pressure_reduced_viscosity(tr_m, fp0_m, fq0_m);
    const auto z2 = Base::high_pressure_reduced_viscosity(z1, pr_m, tr_m);
    const auto fp = Base::polarity_factor(fp0_m, z1, z2);
    const auto xi = Base::inverse_reduced_viscosity(pc_m, tc_m, mw_m);

    return z2 * fp / xi * 1e-7;
  }

  static constexpr bool requires_aligned_alloc =
      (N != dynamic_extent && (sizeof(Vector) % 16) == 0);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(requires_aligned_alloc)

 private:
  Vector pc_;  /// Critical pressure
  Vector tc_;  /// Critical temperature
  Vector vc_;  /// Critical volume
  Vector zc_;  /// Critical Z-factor
  Vector mw_;  /// Molecular weight
  Vector dm_;  /// Dipole moment
  Vector q_;   /// Quantum parameter
};

template <typename T>
class LucasMethod<T, 1> : public LucasMethodBase<T> {
 public:
  using Base = LucasMethodBase<T>;

  /// @brief Constructs object
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temprature
  /// @param[in] vc Critical volume
  /// @param[in] zc Critical z-factor
  /// @param[in] mw Molecular weight
  /// @param[in] dm Dipole moment
  /// @param[in] q Quantum parameter for quantum gases
  LucasMethod(const T& pc, const T& tc, const T& vc, const T& zc, const T& mw,
              const T& dm, const T& q)
      : pc_{pc}, tc_{tc}, vc_{vc}, zc_{zc}, mw_{mw}, dm_{dm}, q_{q} {}

  /// @brief Computes gas viscosity
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  T viscosity(const T& p, const T& t) const noexcept {
    const auto tr = t / tc_;
    const auto dmr = Base::reduced_dipole_moment(dm_, pc_, tc_);
    const auto fp0 = Base::low_pressure_polarity_factor(dmr, zc_, tr);
    const auto fq0 = Base::low_pressure_quantum_factor(q_, tr, mw_);
    const auto pr = p / pc_;

    const auto z1 = Base::low_pressure_reduced_viscosity(tr, fp0, fq0);
    const auto z2 = Base::high_pressure_reduced_viscosity(z1, pr, tr);
    const auto fp = Base::polarity_factor(fp0, z1, z2);
    const auto xi = Base::inverse_reduced_viscosity(pc_, tc_, mw_);

    return z2 * fp / xi * 1e-7;
  }

 private:
  T pc_;  /// Critical pressure
  T tc_;  /// Critical temperature
  T vc_;  /// Critical volume
  T zc_;  /// Critical Z-factor
  T mw_;  /// Molecular weight
  T dm_;  /// Dipole moment
  T q_;   /// Quantum parameter
};

}  // namespace eos