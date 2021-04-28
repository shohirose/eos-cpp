#include <gtest/gtest.h>

#include <iostream>

#include "eos/cubic_eos/peng_robinson_eos.hpp"
#include "eos/cubic_eos/soave_redlich_kwong_eos.hpp"
#include "eos/cubic_eos/van_der_waals_eos.hpp"

TEST(CubicEosTest, VanDerWaalsEosTest) {
  // Methane
  const double pc = 4e6;    // Critical pressure [Pa]
  const double tc = 190.6;  // Critical temperature [K]

  auto eos = eos::make_van_der_waals_eos(pc, tc);
  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  {
    const auto state = eos.create_isobaric_isothermal_state(p, t);
    const auto z = state.zfactor();

    ASSERT_EQ(z.size(), 3);

    EXPECT_NEAR(z[0], 0.207498, 1e-6);
    EXPECT_NEAR(z[1], 0.275339, 1e-6);
    EXPECT_NEAR(z[2], 0.616434, 1e-6);

    EXPECT_NEAR(state.fugacity_coeff(z[0]), 0.756747, 1e-6);
    EXPECT_NEAR(state.fugacity_coeff(z[1]), 0.758617, 1e-6);
    EXPECT_NEAR(state.fugacity_coeff(z[2]), 0.741050, 1e-6);
  }

  {
    const auto state = eos.create_isothermal_state(t);
    EXPECT_NEAR(state.pressure(0.001), 1.309708e6, 1.0);
    EXPECT_NEAR(state.pressure(0.01), 1.477564e5, 0.1);
    EXPECT_NEAR(state.pressure(0.1), 1.494696e4, 0.01);
  }
}

TEST(CubicEosTest, SoaveRedlichKwongEosTest) {
  // Methane
  const double pc = 4e6;       // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.008;  // Acentric factor

  auto eos = eos::make_soave_redlich_kwong_eos(pc, tc, omega);

  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  {
    const auto state = eos.create_isobaric_isothermal_state(p, t);
    const auto z = state.zfactor();

    ASSERT_EQ(z.size(), 3);

    EXPECT_NEAR(z[0], 0.152443, 1e-6);
    EXPECT_NEAR(z[1], 0.310673, 1e-6);
    EXPECT_NEAR(z[2], 0.536884, 1e-6);

    EXPECT_NEAR(state.fugacity_coeff(z[0]), 0.69289, 1e-5);
    EXPECT_NEAR(state.fugacity_coeff(z[1]), 0.70862, 1e-5);
    EXPECT_NEAR(state.fugacity_coeff(z[2]), 0.70353, 1e-5);
  }

  {
    const auto state = eos.create_isothermal_state(t);
    EXPECT_NEAR(state.pressure(0.001), 1.283055e6, 1.0);
    EXPECT_NEAR(state.pressure(0.01), 1.474262e5, 0.1);
    EXPECT_NEAR(state.pressure(0.1), 1.494359e4, 0.01);
  }
}

TEST(CubicEosTest, PengRobinsonEosTest) {
  // Methane
  const double pc = 4e6;       // Critical pressure [Pa]
  const double tc = 190.6;     // Critical temperature [K]
  const double omega = 0.008;  // Acentric factor

  auto eos = eos::make_peng_robinson_eos(pc, tc, omega);
  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  {
    const auto state = eos.create_isobaric_isothermal_state(p, t);
    const auto z = state.zfactor();

    ASSERT_EQ(z.size(), 3);

    EXPECT_NEAR(z[0], 0.135628, 1e-6);
    EXPECT_NEAR(z[1], 0.292355, 1e-6);
    EXPECT_NEAR(z[2], 0.510231, 1e-6);

    EXPECT_NEAR(state.fugacity_coeff(z[0]), 0.67210, 1e-5);
    EXPECT_NEAR(state.fugacity_coeff(z[1]), 0.68819, 1e-5);
    EXPECT_NEAR(state.fugacity_coeff(z[2]), 0.68362, 1e-5);
  }

  {
    const auto state = eos.create_isothermal_state(t);
    EXPECT_NEAR(state.pressure(0.001), 1.267541e6, 1.0);
    EXPECT_NEAR(state.pressure(0.01), 1.472064e5, 0.1);
    EXPECT_NEAR(state.pressure(0.1), 1.494132e4, 0.01);
  }
}
