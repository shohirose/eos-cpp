#pragma once

namespace eos {

/**
 * @brief Identity correction factor for attraction parameter
 *
 * This correction factor applies no temperature correction to attraction
 * parameter, which means \f$\alpha=1\f$.
 */
template <typename Scalar>
struct IdentityCorrectionFactor {
  constexpr Scalar value(const Scalar&) const noexcept { return 1; }
  constexpr Scalar derivative(const Scalar&) const noexcept { return 0; }
};

}  // namespace eos
