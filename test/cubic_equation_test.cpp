#include "eos/math/cubic_equation.hpp"

#include <gtest/gtest.h>

TEST(CubicEquationTest, RealRootsTest) {
  // x^3 -x = (x - 1)(x + 1)x =0
  {
    eos::cubic_equation eq(0.0, -1.0, 0.0);
    const auto x = eq.real_roots();
    ASSERT_EQ(x.size(), 3);
    EXPECT_NEAR(x[0], -1.0, 1e-6);
    EXPECT_NEAR(x[1], 0.0, 1e-6);
    EXPECT_NEAR(x[2], 1.0, 1e-6);
  }

  // x^3 - 1 = (x - 1)(x^2 + x + 1) =0
  {
    eos::cubic_equation eq(0.0, 0.0, -1.0);
    const auto x = eq.real_roots();
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }

  // x^3 - x^2 + x - 1 = (x - 1)^3 = 0
  {
    eos::cubic_equation eq(-1.0, 1.0, -1.0);
    const auto x = eq.real_roots();
    ASSERT_EQ(x.size(), 1);
    EXPECT_NEAR(x[0], 1.0, 1e-6);
  }
}