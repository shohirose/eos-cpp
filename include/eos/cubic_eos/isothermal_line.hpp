#pragma once

namespace eos {

template <typename Eos>
class isothermal_line {
 public:
  /// @param[in] t Temperature
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  isothermal_line(double t, double a, double b) noexcept
      : t_{t}, a_{a}, b_{b} {}

  isothermal_line() = default;
  isothermal_line(const isothermal_line &) = default;
  isothermal_line(isothermal_line &&) = default;

  isothermal_line &operator=(const isothermal_line &) = default;
  isothermal_line &operator=(isothermal_line &&) = default;

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
