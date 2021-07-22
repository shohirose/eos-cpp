#pragma once

#include <cmath>  // std::sqrt

namespace eos {

template <typename Scalar>
class SoaveCorrectionPolicy {
 public:
  SoaveCorrectionPolicy() = default;
  SoaveCorrectionPolicy(const Scalar& m) : m_{m} {}

  SoaveCorrectionPolicy(const SoaveCorrectionPolicy&) = default;
  SoaveCorrectionPolicy(SoaveCorrectionPolicy&&) = default;

  SoaveCorrectionPolicy& operator=(const SoaveCorrectionPolicy&) = default;
  SoaveCorrectionPolicy& operator=(SoaveCorrectionPolicy&&) = default;

  const Scalar& m() const noexcept { return m_; }
  Scalar& m() noexcept { return m_; }

  /// @param[in] tr Reduced temperature
  Scalar value(const Scalar& tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /// @param[in] tr Reduced temperature
  Scalar derivative(const Scalar& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    return -m_ * sqrt_tr * (1 + m_ * (1 - sqrt_tr));
  }

 private:
  Scalar m_;
};

}  // namespace eos
