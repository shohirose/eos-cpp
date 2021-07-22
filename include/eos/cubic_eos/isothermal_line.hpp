#pragma once

namespace eos {

template <typename Eos>
class IsothermalLine {
 public:
  /// @param[in] t Temperature
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  IsothermalLine(double t, double a, double b) noexcept
      : t_{t}, a_{a}, b_{b} {}

  IsothermalLine() = default;
  IsothermalLine(const IsothermalLine &) = default;
  IsothermalLine(IsothermalLine &&) = default;

  IsothermalLine &operator=(const IsothermalLine &) = default;
  IsothermalLine &operator=(IsothermalLine &&) = default;

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] v Volume
  double pressure(double v) const noexcept {
    return Eos::pressure(t_, v, a_, b_);
  }

 private:
  double t_;  /// Temperature
  double a_;  /// Attraction parameter
  double b_;  /// Repulsion parameter
};

}  // namespace eos
