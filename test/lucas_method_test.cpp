#include <gtest/gtest.h>

#include "eos/lucas_method.hpp"

// The following unit tests are taken from examples in:
// Poling et al. 2001. "The Properties of Gases and Liquids", fifth edition,
// McGRAW-HILL.

// Example 9-4
TEST(LucasMethodTest, LowPressureViscosityTest) {
  // Methanol
  const double pc = 8.097e6;  // Critical pressure [Pa]
  const double tc = 512.64;   // Critical temperature [K]
  const double mw = 32.042;   // Molecular weight [kg/kmol]
  const double zc = 0.224;    // Critical Z-factor
  const double dm = 1.7;      // Dipole moment [debyes]
  const double q = 0.0;       // Quantum parameter
  const double t = 550.0;     // Temperature [K]

  const auto lucas = eos::make_lucas_method(pc, tc, zc, mw, dm, q);
  const auto visc = lucas.viscosity_at_low_pressure(t);
  EXPECT_NEAR(visc, 178.0e-7, 1.0e-7);
}

// Example 9-10
TEST(LucasMethodTest, HighPressureViscosityTest) {
  // Ammonia
  const double pc = 113.53 * 1e5;  // Critical pressure [Pa]
  const double tc = 405.5;         // Critical temperature [K]
  const double mw = 17.031;        // Molecular weight [g/mol]
  const double zc = 0.244;         // Critical z-factor
  const double dm = 1.47;          // Dipole moment [debyes]
  const double q = 0.0;            // Quantum parameter
  const double p = 300.0 * 1e5;    // Pressure [Pa]
  const double t = 420.0;          // Temperature [K]
  
  const auto lucas = eos::make_lucas_method(pc, tc, zc, mw, dm, q);
  const auto visc = lucas.viscosity_at_high_pressure(p, t);
  EXPECT_NEAR(visc, 602.0 * 1e-6 * 0.1, 1e-7);
}

/*
// Example 9-7
TEST(LowPressureTest, MultiComponentsTest) {
  using Eigen::Vector2d;

  // Ammonia, Hydrogen
  const Vector2d pc = {113.5e5, 13.0e5};   // Critical pressure [Pa]
  const Vector2d tc = {405.5, 33.2};       // Critical temperature [K]
  const Vector2d mw = {17.031, 2.016};     // Molecular weight [g/mol]
  const Vector2d zc = {0.244, 0.306};      // Critical Z-factor
  const Vector2d vc = {72.5e-6, 64.3e-6};  // Critical volume [m3]
  const Vector2d dm = {1.47, 0.0};         // Dipole moment [debyes]
  const Vector2d q = {0.0, 0.76};          // Quantum parameter
  const double t = 33 + 273.15;            // Temperature [K]

  const lucas::low_pressure::Lucas<double, 2> lucas(pc, tc, vc, zc, mw, dm, q);
  const Vector2d x = {67.7e-2, 32.3e-2};    // Composition
  const auto visc = lucas.viscosity(t, x);  // [Pa-s]
  EXPECT_NEAR(visc, 116.1e-7, 0.2e-7);
}

// Original test case
TEST(HighPressureTest, MultiComponentsTest) {
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
*/