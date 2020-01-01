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

#ifndef SHIROSE_EOS_HPP
#define SHIROSE_EOS_HPP

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

#include "shirose/roots.hpp"

namespace shirose {

constexpr double gas_constant = 8.31446261815324;

template <typename T>
class van_der_waals {
 public:
  using value_type = T;

  /// Constant for attraction parameter
  static constexpr double omega_a = 0.421875;
  /// Constant for respulsion parameter
  static constexpr double omega_b = 0.125;

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  van_der_waals(const T &pc, const T &tc)
      : pc_{pc}, tc_{tc}, ac_{this->ac()}, bc_{this->bc()}, a_{}, b_{} {}

  /// @brief Sets critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void set(const T &pc, const T &tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    ac_ = this->ac();
    bc_ = this->bc();
    a_ = 0;
    b_ = 0;
  }

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    return gas_constant * t / (v - bc_) - ac_ / (v * v);
  }

  /// @brief Computes coefficients of cubic equation
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns An array of Z-factor
  std::array<T, 3> cubic_eq(const T &p, const T &t) noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    a_ = this->a(pr, tr);
    b_ = this->b(pr, tr);
    return {-b_ - 1, a_, -a_ * b_};
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] ar Reduced attraction parameter
  /// @param[in] br Reduced repulsion parameter
  /// @returns Fugacity coefficient
  ///
  /// cubic_eq() must be called in advance, and Z-factor must be a solution of
  /// the cubic equation.
  T fugacity_coeff(const T &z) const noexcept {
    using std::exp;
    using std::log;
    return exp(-log(z - b_) - a_ / z + z - 1);
  }

 private:
  /// @brief Computes attraction parameter at critical condition
  T ac() noexcept {
    return (omega_a * gas_constant * gas_constant) * tc_ * tc_ / pc_;
  }

  /// @brief Computes repulsion parameter at critical condition
  T bc() noexcept { return (omega_b * gas_constant) * tc_ / pc_; }

  /// @brief Computes the reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced attraction parameter
  T a(const T &pr, const T &tr) noexcept { return omega_a * pr / (tr * tr); }

  /// @brief Computes the reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced repulsion parameter
  T b(const T &pr, const T &tr) noexcept { return omega_b * pr / tr; }

  /// Critical pressure
  T pc_;
  /// Critical temperature
  T tc_;
  /// Attraction parameter at critical condition
  T ac_;
  /// Repulsion parameter at critical condition
  T bc_;
  /// Reduced attraction parameter
  T a_;
  /// Reduced repulsion parameter
  T b_;
};

template <typename T>
class soave_redlich_kwong {
 public:
  static constexpr double omega_a = 0.42748;
  static constexpr double omega_b = 0.08664;

  using value_type = T;

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Accentric factor
  soave_redlich_kwong(const T &pc, const T &tc, const T &omega) noexcept
      : pc_{pc},
        tc_{tc},
        omega_{omega},
        m_{this->m()},
        ac_{this->ac()},
        bc_{this->bc()},
        a_{},
        b_{} {}

  /// @brief Sets critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Accentric factor
  void set(const T &pc, const T &tc, const T &omega) noexcept {
    pc_ = pc;
    tc_ = tc;
    omega_ = omega;
    m_ = this->m();
    ac_ = this->ac();
    bc_ = this->bc();
    a_ = 0;
    b_ = 0;
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    const auto alpha = this->alpha(t / tc_);
    return gas_constant * t / (v - bc_) - alpha * ac_ / (v * (v + bc_));
  }

  /// @brief Computes the coefficients of the cubic equation of z-factor.
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns Coefficients of the cubic equation of z-factor.
  std::array<T, 3> cubic_eq(const T &p, const T &t) noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    a_ = this->a(pr, tr);
    b_ = this->b(pr, tr);
    return {-1, a_ - b_ - b_ * b_, -a_ * b_};
  }

  /// @brief Computes the fugacity coefficient
  /// @param[in] z Z-factor
  /// @returns Fugacity coefficient
  ///
  /// cubic_eq() must be called in advance, and z-factor must be one of the
  /// solutions of the cubic equation.
  T fugacity_coeff(const T &z) const noexcept {
    using std::exp;
    using std::log;
    return exp(z - 1 - log(z - b_) - a_ * a_ / b_ * log(b_ / z + 1));
  }

 private:
  /// @brief Correction factor for temperature dependency of attraction
  /// parameter
  T m() const noexcept { return 0.48 + (1.574 - 0.176 * omega_) * omega_; }

  /// @brief Computes attraction parameter
  T ac() const noexcept {
    return (omega_a * gas_constant * gas_constant) * tc_ * tc_ / pc_;
  }

  /// @brief Computes repulsion parameter
  T bc() const noexcept { return (omega_b * gas_constant) * tc_ / pc_; }

  /// @brief Computes the temperature correction factor for the attraction
  /// parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor for the attraction parameter
  T alpha(const T &tr) const noexcept {
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// @brief Computes the reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced attraction parameter
  T a(const T &pr, const T &tr) const noexcept {
    return this->alpha(tr) * omega_a * pr / (tr * tr);
  }

  /// @brief Computes the reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced repulsion parameter
  T b(const T &pr, const T &tr) const noexcept { return omega_b * pr / tr; }

  /// Critical pressure
  T pc_;
  /// Critical temperature
  T tc_;
  /// Acentric factor
  T omega_;
  /// Correction factor for temperature dependency of attraction parameter
  T m_;
  /// Attraction parameter at critical condition
  T ac_;
  /// Repulsion parameter at critical condition
  T bc_;
  /// Reduced attraction parameter
  T a_;
  /// Reduced repulsion parameter
  T b_;
};

template <typename T>
class peng_robinson {
 public:
  /// Constant for attraction parameter
  static constexpr double omega_a = 0.45724;
  /// Constant for repulsion parameter
  static constexpr double omega_b = 0.07780;

  /// Square root of 2 (two).
  static constexpr double sqrt2 = 1.4142135623730950488;
  static constexpr double delta1 = 1 + sqrt2;
  static constexpr double delta2 = 1 - sqrt2;

  using value_type = T;

  /// @brief Constructs EoS from critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  peng_robinson(const T &pc, const T &tc, const T &omega)
      : pc_{pc},
        tc_{tc},
        omega_{omega},
        m_{this->m()},
        ac_{this->ac()},
        bc_{this->bc()},
        a_{},
        b_{} {}

  /// @brief Sets critical parameters
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] omega Acentric factor
  void set(const T &pc, const T &tc, const T &omega) noexcept {
    pc_ = pc;
    tc_ = tc;
    omega_ = omega;
    m_ = this->m();
    ac_ = this->ac();
    bc_ = this->bc();
    a_ = 0;
    b_ = 0;
  }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    const auto alpha = this->alpha(t / tc_);
    return gas_constant * t / (v - bc_) -
           alpha * ac_ / ((v - bc_) * (v + bc_) + 2 * bc_ * v);
  }

  /// @brief Computes coeficients of the cubic equation of z-factor.
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns An array of coefficients of the cubic equation.
  std::array<T, 3> cubic_eq(const T &p, const T &t) noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    a_ = this->a(pr, tr);
    b_ = this->b(pr, tr);
    return {b_ - 1, a_ - (3 * b_ + 2) * b_, (-a_ + b_ + b_ * b_) * b_};
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @returns Fugacity coefficient
  ///
  /// cubic_eq() must be called in advance, and z-factor must be one of the
  /// solutions of the cubic equation.
  T fugacity_coeff(const T &z) const noexcept {
    using std::exp;
    using std::log;
    return exp(z - 1 - log(z - b_) -
               a_ / (2 * sqrt2 * b_) *
                   log((z + delta1 * b_) / (z - delta2 * b_)));
  }

 private:
  /// @brief Correction factor for temperature dependency of attraction
  /// parameter
  T m() const noexcept {
    return 0.3796 + omega_ * (1.485 - omega_ * (0.1644 - 0.01667 * omega_));
  }

  /// @brief Computes attraction parameter
  T ac() const noexcept {
    return (omega_a * gas_constant * gas_constant) * tc_ * tc_ / pc_;
  }

  /// @brief Computes repulsion parameter
  T bc() const noexcept { return (omega_b * gas_constant) * tc_ / pc_; }

  /// @brief Computes the temperature correction factor for the attraction
  /// parameter
  /// @param[in] tr Reduced temperature
  /// @returns Temperature correction factor for the attraction parameter
  T alpha(const T &tr) const noexcept {
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// @brief Computes the reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced attraction parameter
  T a(const T &pr, const T &tr) const noexcept {
    return this->alpha(tr) * omega_a * pr / (tr * tr);
  }

  /// @brief Computes the reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced repulsion parameter
  T b(const T &pr, const T &tr) const noexcept { return omega_b * pr / tr; }

  /// Critical pressure
  T pc_;
  /// Critical temperature
  T tc_;
  /// Acentric factor
  T omega_;
  /// Correction factor for temperature dependency of attraction parameter
  T m_;
  /// Attraction parameter at critical condition
  T ac_;
  /// Repulsion parameter at critical condition
  T bc_;
  /// Reduced attraction parameter
  T a_;
  /// Reduced repulsion parameter
  T b_;
};

template <typename Policy>
class equation_of_state {
 public:
  using value_type = typename Policy::value_type;

  /// @brief Constructs EoS.
  template <typename... Args>
  equation_of_state(Args &&... args) : policy_{std::forward<Args>(args)...} {}

  /// @brief Set parameters.
  template <typename... Args>
  void set(Args &&... args) noexcept {
    policy_.set(std::forward<Args>(args)...);
  }

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  value_type pressure(const value_type &t, const value_type &v) noexcept {
    return policy_.pressure(t, v);
  }

  /// @brief Computes Z-factor
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns An array of Z-factor
  std::vector<value_type> zfactor(const value_type &p,
                                  const value_type &t) noexcept {
    const auto c = policy_.cubic_eq(p, t);
    const auto x = roots(c);
    return real_roots(x);
  }

  /// @brief Computes fugacity coefficient
  /// @param[in] z Z-factor
  /// @returns Fugacity coefficient
  ///
  /// zfactor() must be called in advance, and Z-factor must be a solution of
  /// zfactor().
  value_type fugacity_coeff(const value_type &z) const noexcept {
    return policy_.fugacity_coeff(z);
  }

 private:
  Policy policy_;
};

/// @brief Van der Waals equation of state.
template <typename T>
using vdw_eos = equation_of_state<van_der_waals<T>>;

/// @brief Soave-Redlich-Kwong equation of state.
template <typename T>
using srk_eos = equation_of_state<soave_redlich_kwong<T>>;

/// @brief Peng-Robinson equation of state.
template <typename T>
using pr_eos = equation_of_state<peng_robinson<T>>;

}  // namespace shirose

#endif  // SHIROSE_EOS_HPP