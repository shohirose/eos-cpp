#include <gsl/gsl_poly.h>
#include <gtest/gtest.h>

#include <iostream>

#include "eos/peng_robinson_eos.hpp"
#include "eos/soave_redlich_kwong_eos.hpp"
#include "eos/van_der_waals_eos.hpp"

struct CubicEquationSolver {
  std::vector<double> operator()(const std::array<double, 3>& a) const noexcept {
    std::vector<double> x(3);
    const auto n = gsl_poly_solve_cubic(a[0], a[1], a[2], &x[0], &x[1], &x[2]);
    x.resize(n);
    return x;
  }
};

TEST(CubicEosTest, VanDerWaalsEosTest) {
  // Methane
  // Taken from "Properties of Gases and Liquids, 5th edition"
  const double pc = 4.6e6;  // Critical pressure [Pa]
  const double tc = 190.6;  // Critical temperature [K]

  auto eos = eos::makeVanDerWaalsEos(pc, tc);
  const double p = 3.5e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  const auto [z, params] = eos.zfactor(p, t, CubicEquationSolver{});
  ASSERT_EQ(z.size(), std::size_t{3});
  EXPECT_NEAR(z[0], 0.208227, 1e-6);
  EXPECT_NEAR(z[1], 0.287939, 1e-6);
  EXPECT_NEAR(z[2], 0.604543, 1e-6);
  EXPECT_NEAR(eos.fugacityCoeff(z[0], params), 0.748170, 1e-6);
  EXPECT_NEAR(eos.fugacityCoeff(z[1], params), 0.750815, 1e-6);
  EXPECT_NEAR(eos.fugacityCoeff(z[2], params), 0.736909, 1e-6);

  EXPECT_NEAR(eos.pressure(t, 0.001), 1.333628e6, 1.0);
  EXPECT_NEAR(eos.pressure(t, 0.01), 1.480043e5, 0.1);
  EXPECT_NEAR(eos.pressure(t, 0.1), 1.494944e4, 0.01);
}

TEST(CubicEosTest, SoaveRedlichKwongEosTest) {
  // Methane
  // Taken from "Properties of Gases and Liquids, 5th edition"
  const double pc = 4.6e6;     // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.011;  // Acentric factor

  auto eos = eos::makeSoaveRedlichKwongEos(pc, tc, omega);

  const double p = 3.5e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  const auto [z, params] = eos.zfactor(p, t, CubicEquationSolver{});
  ASSERT_EQ(z.size(), std::size_t{3});
  EXPECT_NEAR(z[0], 0.153854, 1e-6);
  EXPECT_NEAR(z[1], 0.328967, 1e-6);
  EXPECT_NEAR(z[2], 0.517177, 1e-6);
  EXPECT_NEAR(eos.fugacityCoeff(z[0], params), 0.68414, 1e-5);
  EXPECT_NEAR(eos.fugacityCoeff(z[1], params), 0.70153, 1e-5);
  EXPECT_NEAR(eos.fugacityCoeff(z[2], params), 0.69864, 1e-5);

  EXPECT_NEAR(eos.pressure(t, 0.001), 1.309626e6, 1.0);
  EXPECT_NEAR(eos.pressure(t, 0.01), 1.477157e5, 0.1);
  EXPECT_NEAR(eos.pressure(t, 0.1), 1.494651e4, 0.01);
}

TEST(CubicEosTest, PengRobinsonEosTest) {
  // Methane
  // Taken from "Properties of Gases and Liquids, 5th edition"
  const double pc = 4.6e6;     // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.011;  // Acentric factor

  auto eos = eos::makePengRobinsonEos(pc, tc, omega);
  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  const auto [z, params] = eos.zfactor(p, t, CubicEquationSolver{});
  ASSERT_EQ(z.size(), std::size_t{3});
  EXPECT_NEAR(z[0], 0.123511, 1e-6);
  EXPECT_NEAR(z[1], 0.198808, 1e-6);
  EXPECT_NEAR(z[2], 0.623953, 1e-6);
  EXPECT_NEAR(eos.fugacityCoeff(z[0], params), 0.75906, 1e-5);
  EXPECT_NEAR(eos.fugacityCoeff(z[1], params), 0.76541, 1e-5);
  EXPECT_NEAR(eos.fugacityCoeff(z[2], params), 0.72557, 1e-5);

  EXPECT_NEAR(eos.pressure(t, 0.001), 1.295462e6, 1.0);
  EXPECT_NEAR(eos.pressure(t, 0.01), 1.475243e5, 0.1);
  EXPECT_NEAR(eos.pressure(t, 0.1), 1.494454e4, 0.01);
}
