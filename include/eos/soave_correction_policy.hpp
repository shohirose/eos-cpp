#pragma once

#include <cmath>  // std::sqrt

namespace eos {

/**
 * @brief Correction factor proposed by Soave (1972).
 *
 * @tparam Scalar scalar
 */
template <typename Scalar>
class SoaveCorrectionPolicy {
 public:
  SoaveCorrectionPolicy() = default;

  /**
   * @brief Construct a new Soave Correction Policy object
   *
   * @param m coefficient used to calculate a correction factor
   */
  SoaveCorrectionPolicy(const Scalar& m) : m_{m} {}

  SoaveCorrectionPolicy(const SoaveCorrectionPolicy&) = default;
  SoaveCorrectionPolicy(SoaveCorrectionPolicy&&) = default;

  SoaveCorrectionPolicy& operator=(const SoaveCorrectionPolicy&) = default;
  SoaveCorrectionPolicy& operator=(SoaveCorrectionPolicy&&) = default;

  const Scalar& m() const noexcept { return m_; }
  Scalar& m() noexcept { return m_; }

  /**
   * @brief Compute the correction factor at given temperature
   *
   * @param tr reduced temperature
   * @return Scalar correction factor for attraction parameter
   */
  Scalar value(const Scalar& tr) const noexcept {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /**
   * @brief Compute the logarithmic derivative of correction factor
   *
   * This function calculates the following derivative:
   * \f[ \frac{d \ln \alpha}{d \ln T} \f]
   * where \f$\alpha\f$ is the correction factor and \f$T\f$ is temperature.
   *
   * @param tr reduced temperature
   * @return Scalar logarithmic derivative of correction factor
   */
  Scalar derivative(const Scalar& tr) const noexcept {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    return -m_ * sqrt_tr * (1 + m_ * (1 - sqrt_tr));
  }

 private:
  Scalar m_;
};

}  // namespace eos
