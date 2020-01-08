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

#ifndef SHIROSE_ALPHA_FUNCTIONS_HPP
#define SHIROSE_ALPHA_FUNCTIONS_HPP

#include <cmath>

namespace shirose {

namespace alpha {

template <typename T>
struct no_correction {
  T value(const T&) const noexcept { return 1; }
  T derivative(const T&) const noexcept { return 0; }
  T second_derivative(const T&) const noexcept { return 0; }
};

template <typename T>
class soave_1972 {
 public:
  static T m(const T& omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  soave_1972(const T& omega) : m_{this->m(omega)} {}

  T value(const T& tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  T derivative(const T& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -a * m_ / sqrt_tr;
  }

  T second_derivative(const T& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return a * m_ * m_ / tr + a * m_ / (2 * tr * sqrt_tr);
  }

 private:
  T m_;
};

template <typename T>
class peng_robinson_1976 {
 public:
  static T m(const T& omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  peng_robinson_1976(const T& omega) : m_{this->m(omega)} {}

  T value(const T& tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  T derivative(const T& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -a * m_ / sqrt_tr;
  }

  T second_derivative(const T& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return a * m_ * m_ / tr + a * m_ / (2 * tr * sqrt_tr);
  }

 private:
  T m_;
};

}  // namespace alpha

}  // namespace shirose

#endif  // SHIROSE_ALPHA_FUNCTIONS_HPP