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

#include <gtest/gtest.h>

#include "eos/eos.hpp"

using namespace eos;

// The following unit tests are taken from examples in:
// Poling et al. 2001. "The Properties of Gases and Liquids", fifth edition,
// McGRAW-HILL.

// Example 9-4
TEST(LowPressureTest, PureComponentTest) {
  // Methanol
  const double pc = 8.097e6;  // Critical pressure [Pa]
  const double tc = 512.64;   // Critical temperature [K]
  const double mw = 32.042;   // Molecular weight [kg/kmol]
  const double zc = 0.224;    // Critical Z-factor
  const double dm = 1.7;      // Dipole moment [debyes]
  const double q = 0.0;       // Quantum parameter
  const double t = 550.0;     // Temperature [K]
  const auto tr = t / tc;     // Reduced temperature

  const auto dmr = lucas::reduced_dipole_moment(dm, tc, pc);
  EXPECT_NEAR(dmr, 4.67e-2, 1e-2);

  const auto fp0 = lucas::low_pressure::polarity_factor(dmr, zc, tr);
  EXPECT_NEAR(fp0, 1.30, 1e-2);

  const auto xi = lucas::inverse_viscosity(pc, tc, mw);
  EXPECT_NEAR(xi, 4.70e4, 0.01e3);

  const auto mur = lucas::low_pressure::reduced_viscosity(tr, fp0, 1.0);
  // mur is 0.836 in the original example but it is actually 0.838.
  EXPECT_NEAR(mur, 0.838, 1e-3);

  const lucas::low_pressure::Lucas<double, 1> lucas(pc, tc, zc, mw, dm, q);
  const auto visc = lucas.viscosity(t);
  EXPECT_NEAR(visc, 178.0e-7, 1.0e-7);
}

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

// Example 9-10
TEST(HighPressureTest, PureComponentTest) {
  // Ammonia
  const double pc = 113.53 * 1e5;  // Critical pressure [Pa]
  const double tc = 405.5;         // Critical temperature [K]
  const double mw = 17.031;        // Molecular weight [g/mol]
  const double zc = 0.244;         // Critical z-factor
  const double dm = 1.47;          // Dipole moment [debyes]
  const double q = 0.0;            // Quantum factor
  const double p = 300.0 * 1e5;    // Pressure [Pa]
  const double t = 420.0;          // Temperature [K]
  const auto pr = p / pc;
  const auto tr = t / tc;

  const auto dmr = lucas::reduced_dipole_moment(dm, tc, pc);
  EXPECT_NEAR(dmr, 7.827e-2, 0.001e-2);

  const auto fp0 = lucas::low_pressure::polarity_factor(dmr, zc, tr);
  EXPECT_NEAR(fp0, 1.164, 1e-3);

  const auto xi = lucas::inverse_viscosity(pc, tc, mw);
  EXPECT_NEAR(xi, 4.95e4, 0.01e4);

  const auto z1 = lucas::low_pressure::reduced_viscosity(tr, fp0, 1.0);
  // z1 is 0.7259 in the original example but it is actually 0.7256.
  EXPECT_NEAR(z1, 0.7259, 3e-4);

  const auto z2 = lucas::high_pressure::reduced_viscosity(z1, pr, tr);
  // z2 is 0.3466 in the original example but it is actually 3.4649.
  EXPECT_NEAR(z2, 3.466, 2e-3);

  const auto fp = lucas::high_pressure::polarity_factor(fp0, z1, z2);
  EXPECT_NEAR(fp, 0.860, 1e-3);

  const lucas::high_pressure::Lucas<double, 1> lucas(pc, tc, zc, mw, dm, q);

  const auto visc = lucas.viscosity(p, t);  // Viscosity [Pa-s]

  EXPECT_NEAR(visc, 602.0 * 1e-6 * 0.1, 1e-7);
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