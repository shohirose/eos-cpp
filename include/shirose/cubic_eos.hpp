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

#ifndef SHIROSE_CUBIC_EOS_HPP
#define SHIROSE_CUBIC_EOS_HPP

#include <vector>  // std::vector

#include "shirose/constants.hpp"  // shirose::gas_constant
#include "shirose/roots.hpp"      // shirose::roots

namespace shirose {

/// @brief Two-parameter cubic equation of state (EoS)
/// @tparam Eos EoS policy
/// @tparam Alpha Temperature correction policy for attraction parameter
///
/// Eos must have the following static functions:
///    - pressure(t, v, a, b)
///    - cubic_eq(ar, br)
///    - fugacity_coeff(z, ar, br)
/// where t is temperature, v is volume, a is attraction parameter, b is
/// repulsion parameter, ar is reduced attraction parameter, br is reduced
/// repulsion parameter, and z is z-factor.
///
/// Alpha must have the following memeber functions:
///    - value(tr), computes alpha
///    - derivative(tr), computes the derivative of alpha
///    - second_derivative(tr), computes the second derivative of alpha
/// where tr is reduced temperature.
template <typename T, typename Eos, typename Alpha>
class cubic_eos {
 public:
  // Static functions

  /// @brief Computes attraction parameter
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @returns Attraction parameter at critical condition
  static T attraction_param(const T &pc, const T &tc) noexcept {
    return (Eos::omega_a * gas_constant * gas_constant) * tc * tc / pc;
  }

  /// @brief Computes repulsion parameter
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @returns Repulsion parameter at critical condition
  static T repulsion_param(const T &pc, const T &tc) noexcept {
    return (Eos::omega_b * gas_constant) * tc / pc;
  }

  /// @brief Computes reduced attraction parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced attraction parameter
  static T reduced_attraction_param(const T &pr, const T &tr) noexcept {
    return Eos::omega_a * pr / (tr * tr);
  }

  /// @brief Computes reduced repulsion parameter
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  /// @returns Reduced repulsion parameter
  static T reduced_repulsion_param(const T &pr, const T &tr) noexcept {
    return Eos::omega_b * pr / tr;
  }

  // Constructors

  /// @brief Constructs EoS.
  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  /// @param[in] alpha Temperature correction policy for attraction parameter
  cubic_eos(const T &pc, const T &tc, const Alpha &alpha)
      : pc_{pc},
        tc_{tc},
        ac_{this->attraction_param(pc, tc)},
        bc_{this->repulsion_param(pc, tc)},
        alpha_{alpha} {}

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @returns Pressure
  T pressure(const T &t, const T &v) noexcept {
    const auto tr = t / tc_;
    const auto a = alpha_.value(tr) * ac_;
    const auto b = bc_;
    return Eos::pressure(t, v, a, b);
  }

  /// @brief A state at given pressure and temperature expressed by this
  /// equation of state
  class pt_state {
   public:
    /// @brief Constructs a state
    /// @param[in] a Reduced attration parameter
    /// @param[in] b Reduced repulsion parameter
    pt_state(const T &a, const T &b) : a_{a}, b_{b} {}

    /// @brief Computes z-factor
    /// @returns An array of z-factors
    std::vector<T> zfactor() const noexcept {
      const auto p = Eos::cubic_eq(a_, b_);
      const auto x = roots(p);
      return real_roots(x);
    }

    /// @brief Computes fugacity coefficient
    /// @param[in] z Z-factor
    /// @returns Fugacity coefficient
    ///
    /// Z-factor must be one of solutions of zfactor() function.
    T fugacity_coeff(const T &z) const noexcept {
      return Eos::fugacity_coeff(z, a_, b_);
    }

   private:
    /// Reduced attraction parameter
    T a_;
    /// Reduced repulsion parameter
    T b_;
  };

  /// @brief Creates a state of given pressure and temperature
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  /// @returns State
  pt_state state(const T &p, const T &t) noexcept {
    const auto pr = p / pc_;
    const auto tr = t / tc_;
    const auto a = alpha_.value(tr) * this->reduced_attraction_param(pr, tr);
    const auto b = this->reduced_repulsion_param(pr, tr);
    return {a, b};
  }

 private:
  /// Critical pressure
  T pc_;
  /// Critical temperature
  T tc_;
  /// Attraction parameter at critical condition
  T ac_;
  /// Repulsion parameter at critical condition
  T bc_;
  /// Temperature correction policy for attraction parameter
  Alpha alpha_;
};

}  // namespace shirose

#endif  // SHIROSE_CUBIC_EOS_HPP