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

#ifndef SHIROSE_CORRECTION_POLICY_HPP
#define SHIROSE_CORRECTION_POLICY_HPP

#include <cmath>

namespace shirose {

namespace policy {

/// @brief No temperature correction for attraction parameter
/// @tparam T Value type
template <typename T>
struct no_correction {
  /// @brief Computes temperature correction factor for attraction parameter
  T alpha(const T&) const noexcept { return 1; }

  /// @brief Computes temperature correction factor for attraction parameter
  ///
  /// \f[ \beta = \frac{T_r}{\alpha} \frac{\mathrm{d}\alpha}{\mathrm{d}T_r} \f]
  T beta(const T&) const noexcept { return 0; }

  /// @brief Computes temperature correction factor for attraction parameter
  ///
  /// \f[ \gamma = \frac{T_r^2}{\alpha}
  /// \frac{\mathrm{d}^2\alpha}{\mathrm{d}T_r^2} \f]
  T gamma(const T&) const noexcept { return 0; }
};

/// @brief Base class for temperature correction proposed by Soave (1972)
/// @tparam T Value type
template <typename T>
class soave_1972_base {
 public:
  /// @brief Constructs temperature correction policy
  /// @param[in] m Correction parameter
  soave_1972_base(const T& m) : m_{m} {}

  /// @brief Computes temperature correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  T alpha(const T& tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// @brief Computes temperature correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  ///
  /// \f[ \beta = \frac{T_r}{\alpha} \frac{\mathrm{d}\alpha}{\mathrm{d}T_r} \f]
  T beta(const T& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -tr * m_ / (a * sqrt_tr);
  }

  /// @brief Computes temperature correction factor for attraction parameter
  /// @param[in] tr Reduced temperature
  ///
  /// \f[ \gamma = \frac{T_r^2}{\alpha}
  /// \frac{\mathrm{d}^2\alpha}{\mathrm{d}T_r^2} \f]
  T gamma(const T& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    const auto alpha = a * a;
    return m_ / (2 * alpha) * (m_ * tr + a * sqrt_tr);
  }

 private:
  T m_;
};

/// @brief Temperature correction proposed by Soave (1972)
/// @tparam T Value type
template <typename T>
class soave_1972 : public soave_1972_base<T> {
 public:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static T m(const T& omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  /// @brief Constructs temperature correction policy
  /// @param[in] omega Acentric factor
  soave_1972(const T& omega) : soave_1972_base<T>{this->m(omega)} {}
};

/// @brief Temperature correction proposed by Peng and Robinson (1976)
/// @tparam T Value type
template <typename T>
class peng_robinson_1976 : public soave_1972_base<T> {
 public:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static T m(const T& omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  /// @brief Constructs temperature correction policy
  /// @param[in] omega Acentric factor
  peng_robinson_1976(const T& omega) : soave_1972_base<T>{this->m(omega)} {}
};

}  // namespace policy

}  // namespace shirose

#endif  // SHIROSE_CORRECTION_POLICY_HPP