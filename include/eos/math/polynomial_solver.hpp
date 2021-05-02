#pragma once

#include <gsl/gsl>  // gsl::span
#include <memory>   // std::unique_ptr
#include <vector>   // std::vector

namespace eos {

/// @brief Solver for polynomials with real coefficients
class polynomial_solver {
 public:
  polynomial_solver();
  polynomial_solver(std::size_t num_coeffs);
  polynomial_solver(const polynomial_solver&);
  polynomial_solver(polynomial_solver&&);

  ~polynomial_solver();

  polynomial_solver& operator=(const polynomial_solver&);
  polynomial_solver& operator=(polynomial_solver&&);

  void reset(std::size_t num_coeffs);

  /// @brief Solve a polynomial for real roots
  /// @param[in] a Coefficients of a polynomial
  /// @returns Real roots in the ascending order
  ///
  /// The polynomial takes the form of:
  /// \f[
  ///  a[0] + a[1] x + a[2] x^2 + ... + a[N-1] x^{N-1} = 0
  /// \f]
  ///
  /// The number of coefficients are automatically modified if it is not
  /// set or different from that of the given coefficients.
  std::vector<double> solve(gsl::span<const double> a);

 private:
  class impl;
  std::unique_ptr<impl> pimpl_;
};

}  // namespace eos
