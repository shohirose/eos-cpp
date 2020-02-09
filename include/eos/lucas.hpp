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

#include <Eigen/Core>   // Eigen::Matrix, Eigen::Index, Eigen::Dynamic
#include <cmath>        // std::log, std::exp, std::fabs, std::pow
#include <type_traits>  // std::conditional

#include "eos/constants.hpp"  // eos::gas_constant, eos::dynamic_extent

namespace eos {

namespace lucas {

/// @brief Computes inverse viscosity.
/// @tparam T Value type
/// @param[in] pc Critical pressure [Pa]
/// @param[in] tc Critical temperature [K]
/// @param[in] mw Molecular weight [kg/kmol]
/// @return Inverse viscosity [1/(Pa-s)]
template <typename T>
T inverse_viscosity(const T& pc, const T& tc, const T& mw) noexcept {
  using std::pow;
  const auto pc_bar = pc * 1e-5;
  const auto pc2 = pc_bar * pc_bar;
  return (1.0e7 * 0.176) * pow(tc / (mw * mw * mw * pc2 * pc2), 1.0 / 6.0);
}

/// @brief Computes reduced dipole moment
/// @tparam T Value type
/// @param[in] dm Dipole moment
/// @param[in] tc Critical temperature
/// @param[in] pc Critical pressure
template <typename T>
T reduced_dipole_moment(const T& dm, const T& tc, const T& pc) noexcept {
  return 52.46e-5 * dm * dm / (tc * tc) * pc;
}

namespace low_pressure {

/// @brief Computes reduced viscosity at low pressure.
/// @tparam T Value type
/// @param[in] tr Reduced temperature
/// @param[in] fp0 Polarity factor at low pressure
/// @param[in] fq0 Quantum factor at low pressure
template <typename T>
T reduced_viscosity(const T& tr, const T& fp0, const T& fq0) noexcept {
  using std::exp;
  using std::pow;
  const auto z1 = (0.807 * pow(tr, 0.618) - 0.357 * exp(-0.449 * tr) +
                   0.340 * exp(-4.058 * tr) + 0.018);
  return z1 * fp0 * fq0;
}

/// @brief Computes polarity factor at low pressure.
/// @tparam T Value type
/// @param[in] dmr Reduced dipole momen [Debyes]
/// @param[in] zc Critical z-factor
/// @param[in] tr Reduced temperature
template <typename T>
T polarity_factor(const T& dmr, const T& zc, const T& tr) noexcept {
  using std::fabs;
  using std::pow;
  assert(dmr >= 0.0);
  if (dmr < 0.022) {
    return 1.0;
  } else if (dmr < 0.075) {
    return 1.0 + 30.55 * pow(0.292 - zc, 1.72);
  } else {
    return 1.0 + 30.55 * pow(0.292 - zc, 1.72) * fabs(0.96 + 0.1 * (tr - 0.7));
  }
}

/// @brief Computes quantum factor at low pressure.
/// @tparam T Value type
/// @param[in] q Quantum factor
/// @param[in] tr Reduced temperature [K]
/// @param[in] mw Molecular weight [kg/kmol]
///
/// Quantum factor is required only for quantum gases: H2, He, and D2.
/// q = 1.38 (He), q = 0.76 (H2), and q = 0.52 (D2).
template <typename T>
T quantum_factor(const T& q, const T& tr, const T& mw) noexcept {
  using std::copysign;
  using std::pow;
  assert(q >= 0.0);
  if (q == 0.0) {
    return 1.0;
  } else {
    const auto tmp = tr - 12.0;
    if (tmp == 0.0) {
      return 1.22 * pow(q, 0.15);
    } else {
      return 1.22 * pow(q, 0.15) *
             (1.0 + copysign(0.00385 * pow(tmp * tmp, 1.0 / mw), tmp));
    }
  }
}

/// @brief Computes gas viscosity at low pressure by using Lucas method.
/// @tparam T Value type
/// @tparam N Number of components
template <typename T, std::size_t N>
class Lucas {
 public:
  /// Vector type
  using Vector = Eigen::Matrix<T, N, 1>;

  /// @brief Constructs object
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temprature
  /// @param[in] vc Critical volume
  /// @param[in] zc Critical z-factor
  /// @param[in] mw Molecular weight
  /// @param[in] dm Dipole moment
  /// @param[in] q Quantum parameter
  Lucas(const Eigen::Ref<const Vector>& pc, const Eigen::Ref<const Vector>& tc,
        const Eigen::Ref<const Vector>& vc, const Eigen::Ref<const Vector>& zc,
        const Eigen::Ref<const Vector>& mw, const Eigen::Ref<const Vector>& dm,
        const Eigen::Ref<const Vector>& q)
      : pc_{pc}, tc_{tc}, vc_{vc}, zc_{zc}, mw_{mw}, dm_{dm}, q_{q} {}

  // Member functions

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure [Pa]
  Vector reduced_pressure(const T& p) const noexcept {
    return (p / pc_.array()).matrix();
  }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature [K]
  Vector reduced_temperature(const T& t) const noexcept {
    return (t / tc_.array()).matrix();
  }

  /// @brief Computes reduced dipole moment
  Vector reduced_dipole_moment() const noexcept {
    return 52.46e-5 *
           ((dm_.array() / tc_.array()).square() * pc_.array()).matrix();
  }

  /// @brief Computes polarity factor at low pressure
  /// @param[in] dmr Reduced dipole moment
  /// @param[in] tr Reduced temperature
  template <typename Derived1, typename Derived2>
  Vector polarity_factor(const Eigen::MatrixBase<Derived1>& dmr,
                         const Eigen::MatrixBase<Derived2>& tr) const noexcept {
    Vector fp0;
    for (Eigen::Index i = 0; i < dmr.size(); ++i) {
      fp0[i] = low_pressure::polarity_factor(dmr[i], zc_[i], tr[i]);
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
  template <typename Derived>
  Vector quantum_factor(const Eigen::MatrixBase<Derived>& tr) const noexcept {
    Vector fq0 = Vector::Ones();
    for (Eigen::Index i = 0; i < q_.size(); ++i) {
      fq0[i] = low_pressure::quantum_factor(q_[i], tr[i], mw_[i]);
      assert(!std::isnan(fq0[i]));
    }
    return fq0;
  }

  template <typename Derived1, typename Derived2, typename Derived3>
  std::tuple<T, T, T, T, T, T> mixing_properties(
      const Eigen::MatrixBase<Derived1>& x,
      const Eigen::MatrixBase<Derived2>& fp0,
      const Eigen::MatrixBase<Derived3>& fq0) const noexcept {
    Eigen::Index i;
    const auto r = mw_.maxCoeff(&i) / mw_.minCoeff();
    T a = 1.0;
    using std::pow;
    if (r > 9.0 && x[i] > 0.05 && x[i] < 0.7) {
      a -= 0.01 * pow(r, 0.87);
    }
    return {x.dot(tc_), x.dot(mw_),     x.dot(zc_),
            x.dot(vc_), x.dot(fp0) * a, x.dot(fq0)};
  }

  /// @brief Computes gas viscosity at low pressure for mixtures
  /// @param[in] t Temperature [K]
  /// @param[in] x Composition
  T viscosity(const T& t, const Eigen::Ref<const Vector>& x) const noexcept {
    const auto tr = this->reduced_temperature(t);
    const auto dmr = this->reduced_dipole_moment();
    const auto fp0 = this->polarity_factor(dmr, tr);
    const auto fq0 = this->quantum_factor(tr);

    const auto [tcm, mwm, zcm, vcm, fp0m, fq0m] =
        this->mixing_properties(x, fp0, fq0);
    const auto pcm = gas_constant * tcm * zcm / vcm;
    const auto trm = t / tcm;

    const auto z1 = low_pressure::reduced_viscosity(trm, fp0m, fq0m);
    const auto xi = lucas::inverse_viscosity(pcm, tcm, mwm);

    return z1 / xi;
  }

  static constexpr bool requires_aligned_alloc =
      (N != dynamic_extent && (sizeof(Vector) % 16) == 0);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(requires_aligned_alloc)

 protected:
  Vector pc_;  /// Critical pressure
  Vector tc_;  /// Critical temperature
  Vector vc_;  /// Critical volume
  Vector zc_;  /// Critical Z-factor
  Vector mw_;  /// Molecular weight
  Vector dm_;  /// Dipole moment
  Vector q_;   /// Quantum parameter
};

/// @brief Partial specialization of Lucas for pure components.
/// @tparam T Value type
template <typename T>
class Lucas<T, 1> {
 public:
  /// @brief Constructs object
  /// @param[in] pc Critical pressure [Pa]
  /// @param[in] tc Critical temprature [K]
  /// @param[in] zc Critical z-factor
  /// @param[in] mw Molecular weight [kg/kmol]
  /// @param[in] dm Dipole moment [Debyes]
  /// @param[in] q Quantum parameter
  Lucas(const T& pc, const T& tc, const T& zc, const T& mw, const T& dm,
        const T& q)
      : pc_{pc}, tc_{tc}, zc_{zc}, mw_{mw}, dm_{dm}, q_{q} {}

  // Member functions

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  T reduced_pressure(const T& p) const noexcept { return p / pc_; }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature
  T reduced_temperature(const T& t) const noexcept { return t / tc_; }

  /// @brief Computes reduced dipole moment
  T reduced_dipole_moment() const noexcept {
    const auto tmp = dm_ / tc_;
    return 52.46e-5 * tmp * tmp * pc_;
  }

  /// @brief Computes polarity factor at low pressure.
  /// @param[in] dmr Reduced dipole moment
  /// @param[in] tr Reduced temperature
  T polarity_factor(const T& dmr, const T& tr) const noexcept {
    return low_pressure::polarity_factor(dmr, zc_, tr);
  }

  /// @brief Computes quantum factor at low pressure.
  /// @param[in] tr Reduced temperature
  ///
  /// Quantum factor is required only for quantum gases: H2, He, and D2.
  /// q = 1.38 (He), q = 0.76 (H2), and q = 0.52 (D2).
  T quantum_factor(const T& tr) const noexcept {
    return low_pressure::quantum_factor(q_, tr, mw_);
  }

  /// @brief Computes inverse viscosity based on the corresponding states
  /// method.
  /// @return Inverse viscosity [1/(Pa-s)]
  T inverse_viscosity() const noexcept {
    return lucas::inverse_viscosity(pc_, tc_, mw_);
  }

  /// @brief Computes gas viscosity at low pressure
  /// @param[in] t Temperature [K]
  /// @return Viscosity [Pa-s]
  T viscosity(const T& t) const noexcept {
    const auto tr = this->reduced_temperature(t);
    const auto dmr = this->reduced_dipole_moment();
    const auto fp0 = this->polarity_factor(dmr, tr);
    const auto fq0 = this->quantum_factor(tr);
    const auto z1 = low_pressure::reduced_viscosity(tr, fp0, fq0);
    const auto xi = this->inverse_viscosity();
    return z1 / xi;
  }

 private:
  T pc_;  /// Critical pressure [Pa]
  T tc_;  /// Critical temperature [K]
  T zc_;  /// Critical Z-factor
  T mw_;  /// Molecular weight [kg/kmol]
  T dm_;  /// Dipole moment [Debyes]
  T q_;   /// Quantum parameter for H2, He, D2
};

}  // namespace low_pressure

namespace high_pressure {

/// @brief Computes reduced viscosity at high pressure
/// @tparam T Value type
/// @param[in] z1 Reduced viscosity at low pressure
/// @param[in] pr Reduced pressure
/// @param[in] tr Reduced temperature
/// @return Reduced viscosity
template <typename T>
T reduced_viscosity(const T& z1, const T& pr, const T& tr) noexcept {
  using std::exp;
  using std::pow;

  if (tr <= 1.0) {
    const auto alpha = 3.262 + 14.98 * pow(pr, 5.508);
    const auto beta = 1.390 + 5.746 * pr;
    return 0.600 + 0.760 * pow(pr, alpha) +
           (6.990 * pow(pr, beta) - 0.6) * (1.0 - tr);
  } else {
    const auto a = 1.245e-3 / tr * exp(5.1726 * pow(tr, -0.3286));
    const auto b = a * (1.6553 * tr - 1.2723);
    const auto c = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
    const auto d = 1.7369 / tr * exp(2.2310 * pow(tr, -7.6351));
    const auto f = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
    const auto z2 = 1.0 + a * pow(pr, 1.3088) /
                              (b * pow(pr, f) + 1.0 / (1.0 + c * pow(pr, d)));
    return z2 * z1;
  }
}

/// @brief Computes polarity factor at high pressure.
/// @tparam T Value type
/// @param[in] fp0 Polarity factor at low pressure
/// @param[in] z1 Reduced viscosity at low pressure
/// @param[in] z2 Reduced viscosity at high pressure
template <typename T>
T polarity_factor(const T& fp0, const T& z1, const T& z2) noexcept {
  const auto y = z2 / z1;
  return (1.0 + (fp0 - 1.0) / (y * y * y)) / fp0;
}

/// @brief Computes quantum factor at high pressure.
/// @tparam T Value type
/// @param[in] fq0 Quantum factor at low pressure
/// @param[in] z1 Reduced viscosity at low pressure
/// @param[in] z2 Reduced viscosity at high pressure
template <typename T>
T quantum_factor(const T& fq0, const T& z1, const T& z2) noexcept {
  const auto y = z2 / z1;
  using std::log;
  const auto tmp = log(y);
  const auto tmp2 = tmp * tmp;
  return (1.0 + (fq0 - 1.0) * (1.0 / y - 0.007 * tmp2 * tmp2)) / fq0;
}

/// @brief Computes gas viscosity at hight pressure by using Lucas method.
/// @tparam T Value type
/// @tparam N Number of components
template <typename T, std::size_t N>
class Lucas : public low_pressure::Lucas<T, N> {
 public:
  /// Base type
  using Base = low_pressure::Lucas<T, N>;
  /// Vector type
  using Vector = typename Base::Vector;

  /// @brief Constructs object
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temprature
  /// @param[in] vc Critical volume
  /// @param[in] zc Critical z-factor
  /// @param[in] mw Molecular weight
  /// @param[in] dm Dipole moment
  /// @param[in] q Quantum parameter
  Lucas(const Eigen::Ref<const Vector>& pc, const Eigen::Ref<const Vector>& tc,
        const Eigen::Ref<const Vector>& vc, const Eigen::Ref<const Vector>& zc,
        const Eigen::Ref<const Vector>& mw, const Eigen::Ref<const Vector>& dm,
        const Eigen::Ref<const Vector>& q)
      : Base{pc, tc, vc, zc, mw, dm, q} {}

  /// @brief Computes gas viscosity at low pressure for mixtures
  /// @param[in] p Pressure [Pa]
  /// @param[in] t Temperature [K]
  /// @param[in] x Composition
  T viscosity(const T& p, const T& t, const Eigen::Ref<const Vector>& x) const
      noexcept {
    const auto tr = this->reduced_temperature(t);
    const auto dmr = this->reduced_dipole_moment();
    const auto fp0 = this->polarity_factor(dmr, tr);
    const auto fq0 = this->quantum_factor(tr);

    const auto [tcm, mwm, zcm, vcm, fp0m, fq0m] =
        this->mixing_properties(x, fp0, fq0);
    const auto pcm = gas_constant * tcm * zcm / vcm;
    const auto prm = p / pcm;
    const auto trm = t / tcm;

    const auto z1 = low_pressure::reduced_viscosity(trm, fp0m, fq0m);
    const auto z2 = high_pressure::reduced_viscosity(z1, prm, trm);
    const auto fp = high_pressure::polarity_factor(fp0m, z1, z2);
    const auto fq = high_pressure::quantum_factor(fq0m, z1, z2);
    const auto xi = lucas::inverse_viscosity(pcm, tcm, mwm);

    return z2 * fp * fq / xi;
  }

  static constexpr bool requires_aligned_alloc =
      (N != dynamic_extent && (sizeof(Vector) % 16) == 0);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(requires_aligned_alloc)
};

/// @brief Partial specialization of Lucas for pure components.
/// @tparam T Value type
template <typename T>
class Lucas<T, 1> : public low_pressure::Lucas<T, 1> {
 public:
  using Base = low_pressure::Lucas<T, 1>;

  /// @brief Constructs object
  /// @param[in] pc Critical pressure [Pa]
  /// @param[in] tc Critical temprature [K]
  /// @param[in] zc Critical z-factor
  /// @param[in] mw Molecular weight [kg/kmol]
  /// @param[in] dm Dipole moment [Debyes]
  /// @param[in] q Quantum parameter
  Lucas(const T& pc, const T& tc, const T& zc, const T& mw, const T& dm,
        const T& q)
      : Base{pc, tc, zc, mw, dm, q} {}

  // Member functions

  /// @brief Computes gas viscosity at high pressure
  /// @param[in] p Pressure [Pa]
  /// @param[in] t Temperature [K]
  /// @return Viscosity [Pa-s]
  T viscosity(const T& p, const T& t) const noexcept {
    const auto pr = this->reduced_pressure(p);
    const auto tr = this->reduced_temperature(t);
    const auto dmr = this->reduced_dipole_moment();
    const auto fp0 = this->Base::polarity_factor(dmr, tr);
    const auto fq0 = this->Base::quantum_factor(tr);

    const auto z1 = low_pressure::reduced_viscosity(tr, fp0, fq0);
    const auto z2 = high_pressure::reduced_viscosity(z1, pr, tr);
    const auto fp = high_pressure::polarity_factor(fp0, z1, z2);
    const auto fq = high_pressure::quantum_factor(fq0, z1, z2);
    const auto xi = this->inverse_viscosity();

    return z2 * fp * fq / xi;
  }
};

}  // namespace high_pressure

}  // namespace lucas

}  // namespace eos