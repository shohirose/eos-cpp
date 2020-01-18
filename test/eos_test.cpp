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
  auto flash = eos::flash<pr_eos<double>>(eos);
  const auto result = flash.vapor_pressure(p_init, t);
  const auto pvap = result.first;
  const auto report = result.second;

  EXPECT_NEAR(pvap, 2.87515e6, 10);
  EXPECT_EQ(static_cast<int>(report.error), 0);
}