#pragma once

#include <cmath>  // std::sqrt

namespace eos {

/**
 * @brief Correction factor for attraction parameter proposed by Soave (1972).
 *
 * Soave (1972) proposed the following correction factor:
 * \f[
 *    \alpha(T_r) := \left[ 1 + m \left( 1 - \sqrt{T_r} \right) \right]^2
 * \f]
 * where \f$T_r\f$ is reduced temperature and \f$m\f$ is a coefficient.
 * Then, its logarithmic derivative in terms of temperature becomes
 * \f[
 *    \frac{d \ln \alpha}{d \ln T} = - m \cdot \sqrt{\frac{T_r}{\alpha}}
 * \f]
 * 
 * @tparam Scalar scalar
 */
template <typename Scalar>
class SoaveCorrectionFactor {
 public:
  SoaveCorrectionFactor() = default;

  /**
   * @brief Construct a new SoaveCorrectionFactor object
   *
   * @param m coefficient used to calculate a correction factor
   */
  SoaveCorrectionFactor(const Scalar& m) : m_{m} {}

  SoaveCorrectionFactor(const SoaveCorrectionFactor&) = default;
  SoaveCorrectionFactor(SoaveCorrectionFactor&&) = default;

  SoaveCorrectionFactor& operator=(const SoaveCorrectionFactor&) = default;
  SoaveCorrectionFactor& operator=(SoaveCorrectionFactor&&) = default;

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
    return -m_ * sqrt_tr / (1 + m_ * (1 - sqrt_tr));
  }

 private:
  Scalar m_;
};

}  // namespace eos
