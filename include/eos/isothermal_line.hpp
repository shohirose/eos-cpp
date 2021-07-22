#pragma once

namespace eos {

template <typename Eos>
class CubicEosTraits;

template <typename Eos>
class IsothermalLine {
 public:
  using Scalar = typename CubicEosTraits<Eos>::Scalar;

  /// @param[in] t Temperature
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  IsothermalLine(const Scalar &t, const Scalar &a, const Scalar &b) noexcept
      : t_{t}, a_{a}, b_{b} {}

  IsothermalLine() = default;
  IsothermalLine(const IsothermalLine &) = default;
  IsothermalLine(IsothermalLine &&) = default;

  IsothermalLine &operator=(const IsothermalLine &) = default;
  IsothermalLine &operator=(IsothermalLine &&) = default;

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] v Volume
  Scalar pressure(const Scalar &v) const noexcept {
    return Eos::pressure(t_, v, a_, b_);
  }

 private:
  Scalar t_;  /// Temperature
  Scalar a_;  /// Attraction parameter
  Scalar b_;  /// Repulsion parameter
};

}  // namespace eos
