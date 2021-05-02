#include "eos/math/polynomial.hpp"
#include "eos/math/polynomial_solver.hpp"

#include <gtest/gtest.h>

#include <vector>

TEST(RealRootsTest, CubicEquationTest) {
  // x^3 -x = (x - 1)(x + 1)x =0
  {
    const auto x = eos::real_roots(0.0, -1.0, 0.0);
    ASSERT_EQ(x.size(), 3);
    EXPECT_NEAR(x[0], -1.0, 1e-6);
    EXPECT_NEAR(x[1], 0.0, 1e-6);
    EXPECT_NEAR(x[2], 1.0, 1e-6);
  }

  // x^3 - 1 = (x - 1)(x^2 + x + 1) =0
  {
    const auto x = eos::real_roots(0.0, 0.0, -1.0);
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }

  // x^3 - x^2 + x - 1 = (x - 1)^3 = 0
  {
    const auto x = eos::real_roots(-1.0, 1.0, -1.0);
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }
}

TEST(RealRootsTest, GeneralPolynomialEquationTest) {
  // x^4 + 2x^3 - x^2 - 2x = (x + 2)(x - 1)(x + 1)x = 0
  {
    const std::vector<double> a = {0.0, -2.0, -1.0, 2.0, 1.0};
    const auto x = eos::real_roots(a);
    ASSERT_EQ(x.size(), 4);
    EXPECT_NEAR(x[0], -2.0, 1e-6);
    EXPECT_NEAR(x[1], -1.0, 1e-6);
    EXPECT_NEAR(x[2], 0.0, 1e-6);
    EXPECT_NEAR(x[3], 1.0, 1e-6);
  }

  // x^3 -x = (x - 1)(x + 1)x = 0
  {
    const std::vector<double> a = {0.0, -1.0, 0.0, 1.0};
    const auto x = eos::real_roots(a);
    ASSERT_EQ(x.size(), 3);
    EXPECT_NEAR(x[0], -1.0, 1e-6);
    EXPECT_NEAR(x[1], 0.0, 1e-6);
    EXPECT_NEAR(x[2], 1.0, 1e-6);
  }

  // x^3 - 1 = (x - 1)(x^2 + x + 1) =0
  {
    const std::vector<double> a = {-1.0, 0.0, 0.0, 1.0};
    const auto x = eos::real_roots(a);
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }

  // x^3 - x^2 + x - 1 = (x - 1)^3 = 0
  {
    const std::vector<double> a = {-1.0, 1.0, -1.0, 1.0};
    const auto x = eos::real_roots(a);
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }
}

TEST(PolynomialSolverTest, GeneralPolynomialEquationTest)

{
  eos::polynomial_solver solver;

  // x^4 + 2x^3 - x^2 - 2x = (x + 2)(x - 1)(x + 1)x = 0
  {
    const std::vector<double> a = {0.0, -2.0, -1.0, 2.0, 1.0};
    const auto x = solver.solve(a);
    ASSERT_EQ(x.size(), 4);
    EXPECT_NEAR(x[0], -2.0, 1e-6);
    EXPECT_NEAR(x[1], -1.0, 1e-6);
    EXPECT_NEAR(x[2], 0.0, 1e-6);
    EXPECT_NEAR(x[3], 1.0, 1e-6);
  }

  // x^3 -x = (x - 1)(x + 1)x = 0
  {
    const std::vector<double> a = {0.0, -1.0, 0.0, 1.0};
    const auto x = solver.solve(a);
    ASSERT_EQ(x.size(), 3);
    EXPECT_NEAR(x[0], -1.0, 1e-6);
    EXPECT_NEAR(x[1], 0.0, 1e-6);
    EXPECT_NEAR(x[2], 1.0, 1e-6);
  }

  // x^3 - 1 = (x - 1)(x^2 + x + 1) =0
  {
    const std::vector<double> a = {-1.0, 0.0, 0.0, 1.0};
    const auto x = solver.solve(a);
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }

  // x^3 - x^2 + x - 1 = (x - 1)^3 = 0
  {
    const std::vector<double> a = {-1.0, 1.0, -1.0, 1.0};
    const auto x = solver.solve(a);
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }
}
