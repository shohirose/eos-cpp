// MIT License
//
// Copyright (c) 2019 Sho Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "eos/eos.hpp"

#include <gtest/gtest.h>

#include <iostream>

using namespace eos;

TEST(EosTest, VanDerWaalsEosTest) {
  // Methane
  const double pc = 4e6;    // Critical pressure [Pa]
  const double tc = 190.6;  // Critical temperature [K]

  auto eos = make_vdw_eos(pc, tc);
  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]
  const auto state = eos.state(p, t);
  const auto z = state.zfactor();

  ASSERT_EQ(z.size(), 3);

  EXPECT_NEAR(z[0], 0.616434, 1e-6);
  EXPECT_NEAR(z[1], 0.207498, 1e-6);
  EXPECT_NEAR(z[2], 0.275339, 1e-6);

  EXPECT_NEAR(state.fugacity_coeff(z[0]), 0.741050, 1e-6);
  EXPECT_NEAR(state.fugacity_coeff(z[1]), 0.756747, 1e-6);
  EXPECT_NEAR(state.fugacity_coeff(z[2]), 0.758617, 1e-6);

  EXPECT_NEAR(eos.pressure(t, 0.001), 1.309708e6, 1.0);
  EXPECT_NEAR(eos.pressure(t, 0.01), 1.477564e5, 0.1);
  EXPECT_NEAR(eos.pressure(t, 0.1), 1.494696e4, 0.01);
}

TEST(EosTest, SoaveRedlichKwongEosTest) {
  // Methane
  const double pc = 4e6;       // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.008;  // Acentric factor

  auto eos = make_srk_eos(pc, tc, omega);

  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]
  const auto state = eos.state(p, t);
  const auto z = state.zfactor();

  ASSERT_EQ(z.size(), 3);

  EXPECT_NEAR(z[0], 0.536884, 1e-6);
  EXPECT_NEAR(z[1], 0.152443, 1e-6);
  EXPECT_NEAR(z[2], 0.310673, 1e-6);

  EXPECT_NEAR(state.fugacity_coeff(z[0]), 0.70353, 1e-5);
  EXPECT_NEAR(state.fugacity_coeff(z[1]), 0.69289, 1e-5);
  EXPECT_NEAR(state.fugacity_coeff(z[2]), 0.70862, 1e-5);

  EXPECT_NEAR(eos.pressure(t, 0.001), 1.283055e6, 1.0);
  EXPECT_NEAR(eos.pressure(t, 0.01), 1.474262e5, 0.1);
  EXPECT_NEAR(eos.pressure(t, 0.1), 1.494359e4, 0.01);
}

TEST(EosTest, PengRobinsonEosTest) {
  // Methane
  const double pc = 4e6;       // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.008;  // Acentric factor

  auto eos = make_pr_eos(pc, tc, omega);
  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]
  const auto state = eos.state(p, t);
  const auto z = state.zfactor();

  ASSERT_EQ(z.size(), 3);

  EXPECT_NEAR(z[0], 0.510231, 1e-6);
  EXPECT_NEAR(z[1], 0.135628, 1e-6);
  EXPECT_NEAR(z[2], 0.292355, 1e-6);

  EXPECT_NEAR(state.fugacity_coeff(z[0]), 0.68362, 1e-5);
  EXPECT_NEAR(state.fugacity_coeff(z[1]), 0.67210, 1e-5);
  EXPECT_NEAR(state.fugacity_coeff(z[2]), 0.68819, 1e-5);

  EXPECT_NEAR(eos.pressure(t, 0.001), 1.267541e6, 1.0);
  EXPECT_NEAR(eos.pressure(t, 0.01), 1.472064e5, 0.1);
  EXPECT_NEAR(eos.pressure(t, 0.1), 1.494132e4, 0.01);
}

TEST(FlashTest, VaporPressureTest) {
  // Methane
  const double pc = 4e6;       // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.008;  // Acentric factor

  // Temperature
  const double t = 180;

  const auto p_init = estimate_vapor_pressure(t, pc, tc, omega);
  EXPECT_NEAR(p_init, 2.90772e6, 10);

  auto eos = make_pr_eos(pc, tc, omega);
  auto flash = make_flash(eos);
  auto [pvap, report] = flash.vapor_pressure(p_init, t);

  EXPECT_NEAR(pvap, 2.87515e6, 10);
  EXPECT_EQ(static_cast<int>(report.error), 0);
}

// Unit test using Example 9-10 in Poling et al. 2001. "The Properties of
// Gases and Liquids", fifth edition. McGRAW-HILL.
TEST(LucasMethodTest, PureComponentTest) {
  // Ammonia
  const double pc = 113.53 * 1e5;  // Critical pressure [Pa]
  const double tc = 405.5;         // Critical temperature [K]
  const double mw = 17.031;        // Molecular weight [g/mol]
  const double zc = 0.244;         // Critical z-factor
  const double dm = 1.47;          // Dipole moment [debyes]
  const double q = 0.0;            // Quantum factor

  const lucas::high_pressure::Lucas<double, 1> lucas(pc, tc, zc, mw, dm, q);

  const double p = 300.0 * 1e5;             // Pressure [Pa]
  const double t = 420.0;                   // Temperature [K]
  const auto visc = lucas.viscosity(p, t);  // Viscosity [Pa-s]

  EXPECT_NEAR(visc, 602.0 * 1e-6 * 0.1, 1e-7);
}

TEST(LucasMethodTest, MultiComponentsTest) {
  using Eigen::Vector3d;

  // CH4, N2, CO2
  const Vector3d pc = {4.599e6, 3.398e6, 7.374e6};  // Critical pressure [Pa]
  const Vector3d tc = {190.56, 126.20, 304.12};     // Critical temperature [K]
  const Vector3d vc = {98.06e-6, 90.10e-6, 94.07e-6};  // Critical volume [m3]
  const Vector3d zc = {0.286, 0.289, 0.274};           // Critical z-factor
  const Vector3d mw = {16.043, 28.014, 44.010};  // Molecular weight [g/mol]
  const Vector3d dm = {0.0, 0.0, 0.0};           // Dipole moment [debyes]
  const Vector3d q = {0.0, 0.0, 0.0};            // Quantum factor

  const lucas::high_pressure::Lucas<double, 3> lucas(pc, tc, vc, zc, mw, dm, q);

  const auto p = 7e6;                          // pressure [Pa]
  const auto t = 300.0;                        // temperature [K]
  const Vector3d x = {0.8, 0.15, 0.05};        // Composition [CH4, N2, CO2]
  const auto visc = lucas.viscosity(p, t, x);  // Viscosity [Pa-s]

  EXPECT_NEAR(visc, 1.41055e-5, 1e-10);
}