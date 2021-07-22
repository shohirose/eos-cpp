#pragma once

namespace eos {

template <typename Scalar>
class NoCorrectionPolicy {
 public:
  constexpr Scalar value(const Scalar&) const noexcept { return 1; }
  constexpr Scalar derivative(const Scalar&) const noexcept { return 0; }
};

}  // namespace eos
