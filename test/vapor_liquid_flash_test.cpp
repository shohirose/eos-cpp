#include <gtest/gtest.h>

#include "eos/vapor_liquid_flash.hpp"
#include "eos/peng_robinson_eos.hpp"

TEST(VaporLiquidFlashTest, VaporPressureTest)
{
  using namespace eos;
  // Methane
  const double pc = 4e6;      // Critical pressure [Pa]
  const double tc = 190.6;    // Critical temperature [K]
  const double omega = 0.008; // Acentric factor

  // Temperature
  const double t = 180;

  const auto p0 = eos::estimate_vapor_pressure(t, pc, tc, omega);
  EXPECT_NEAR(p0, 2.90772e6, 10);

  const auto eos = eos::make_peng_robinson_eos(pc, tc, omega);
  auto flash = eos::make_vapor_liquid_flash(eos);
  auto [pvap, result] = flash.vapor_pressure(p0, t);

  EXPECT_NEAR(pvap, 2.87515e6, 10);
  EXPECT_EQ(result.error, decltype(flash)::error_code::success);
}