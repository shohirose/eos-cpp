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
#include <limits>
#include <vector>

#include "shirose/roots.hpp"

namespace shirose {

constexpr double gas_constant = 8.31446261815324;

struct van_der_waals {
  static constexpr double omega_a = 0.421875;
  static constexpr double omega_b = 0.125;

  template <typename T>
  static T m(const T &) noexcept {
    return 0;
  }

  template <typename T>
  static std::array<T, 3> cubic_eq(const T &a, const T &b) noexcept {
    return {-b - 1, a, -a * b};
  }

  template <typename T>
  static T fugacity_coeff(const T &a, const T &b, T z) noexcept {
    using std::exp;
    using std::log;
    return exp(-log(z - b) - a / z + z - 1);
  }

  template <typename T>
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / (v * v);
  }
};

struct soave_redlich_kwong {
  static constexpr double omega_a = 0.42748;
  static constexpr double omega_b = 0.08664;

  template <typename T>
  static T m(T omega) noexcept {
    return 0.48 + 1.574 * omega - 0.176 * omega * omega;
  }

  template <typename T>
  static std::array<T, 3> cubic_eq(const T &a, const T &b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  template <typename T>
  static T fugacity_coeff(const T &a, const T &b, const T &z) noexcept {
    using std::exp;
    using std::log;
    return exp(z - 1 - log(z - b) - a * a / b * log(b / z + 1));
  }

  template <typename T>
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / (v * (v + b));
  }
};

struct peng_robinson {
  static constexpr double omega_a = 0.45724;
  static constexpr double omega_b = 0.07780;

  template <typename T>
  static T m(const T &omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  template <typename T>
  static std::array<T, 3> cubic_eq(const T &a, const T &b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  template <typename T>
  static T fugacity_coeff(const T &a, const T &b, const T &z) noexcept {
    using std::exp;
    using std::log;
    constexpr double sqrt2 = 1.4142135623730950488;
    constexpr double delta1 = 1 + sqrt2;
    constexpr double delta2 = 1 - sqrt2;
    return exp(z - 1 - log(z - b) -
               a / (2 * sqrt2 * b) * log((z + delta1 * b) / (z - delta2 * b)));
  }

  template <typename T>
  static T pressure(const T &t, const T &v, const T &a, const T &b) noexcept {
    return gas_constant * t / (v - b) - a / ((v - b) * (v + b) + 2 * b * v);
  }
};

constexpr auto dynamic_extent = std::numeric_limits<std::size_t>::max();

template <typename Policy, typename T, std::size_t N>
class eos_state {};

template <typename Policy, typename T>
class eos_state<Policy, T, 1> {
 public:
  eos_state(const T &ar, const T &br) : ar_{ar}, br_{br} {}

  std::vector<T> zfactor() const noexcept {
    const auto p = Policy::cubic_eq(ar_, br_);
    const auto x = roots(p);

    std::vector<T> z;
    z.reserve(3);
    using std::fabs;
    constexpr double eps = 1e-10;
    for (auto &&xi : x)
      if (fabs(xi.imag()) < eps) z.push_back(xi.real());
    return z;
  }

  T fugacity_coeff(const T &z) const noexcept {
    return Policy::fugacity_coeff(ar_, br_, z);
  }

 private:
  T ar_;
  T br_;
};

template <typename Policy, typename T, std::size_t N>
class pvt_relation {};

template <typename Policy, typename T>
class pvt_relation<Policy, T, 1> {
 public:
  pvt_relation(const T &a, const T &b) : a_{a}, b_{b} {}

  T pressure(const T &t, const T &v) const noexcept {
    return Policy::pressure(t, v, a_, b_);
  }

 private:
  T a_;
  T b_;
};

template <typename Policy, typename T, std::size_t N>
class equation_of_state {};

template <typename Policy, typename T>
class equation_of_state<Policy, T, 1> {
 public:
  static constexpr auto omega_a = Policy::omega_a;
  static constexpr auto omega_b = Policy::omega_b;

  equation_of_state() = default;

  equation_of_state(const T &pc, const T &tc, const T &omega)
      : pc_{pc}, tc_{tc}, omega_{omega}, m_{Policy::m(omega)} {}

  void set(const T &pc, const T &tc, const T &omega) noexcept {
    pc_ = pc;
    tc_ = tc;
    omega_ = omega;
    m_ = Policy::m(omega);
  }

  pvt_relation<Policy, T, 1> pvt() const noexcept {
    const auto a = (omega_a * gas_constant * gas_constant) * tc_ * tc_ / pc_;
    const auto b = (omega_b * gas_constant) * tc_ / pc_;
    return {a, b};
  }

  eos_state<Policy, T, 1> state(const T &p, const T &t) noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    using std::pow;
    using std::sqrt;
    const auto alpha = pow(1 + m_ * (1 - sqrt(tr)), 2);
    const auto ar = omega_a * alpha * pr / pow(tr, 2);
    const auto br = omega_b * pr / tr;
    return {ar, br};
  }

 private:
  T pc_;
  T tc_;
  T omega_;
  T m_;
};

template <typename T, std::size_t N>
using vdw_eos = equation_of_state<van_der_waals, T, N>;

template <typename T, std::size_t N>
using srk_eos = equation_of_state<soave_redlich_kwong, T, N>;

template <typename T, std::size_t N>
using pr_eos = equation_of_state<peng_robinson, T, N>;

}  // namespace shirose

#endif  // SHIROSE_EOS_HPP