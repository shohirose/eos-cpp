#pragma once

namespace eos {

/**
 * @brief Correction policy for cubic equations of state
 *
 * This correction policy does not apply any modification to attraction
 * parameter.
 *
 * @tparam Scalar
 */
template <typename Scalar>
class NoCorrectionPolicy {
 public:
  constexpr Scalar value(const Scalar&) const noexcept { return 1; }
  constexpr Scalar derivative(const Scalar&) const noexcept { return 0; }
};

}  // namespace eos
