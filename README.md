# eos-cpp

## Overview

eos-cpp is a header-only library for cubic equations of state (EoS) written in c++.

This library provides templated classes for three types of cubic EoS: Van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- `eos::VanDerWaalsEos<Scalar>`
- `eos::SoaveRedlichKwongEos<Scalar>`
- `eos::PengRobinsonEos<Scalar>`

You can customize and create a new cubic EoS by inheriting `CubicEosBase<Scalar, EosPolicy, CorrectionFactor>`. `CubicEosBase` is a template base class for general two-parameter cubic EoS. `CubicEosBase` requires `Scalar`, `EosPoicy`, and `CorrectionFactor` types. While `EosPolicy` class implements functions to calculate thermodynamic properties such as pressure and compressibility, `CorrectionFactor` class provides the calculation method of the correction factor for attraction parameter.

An `EosPolicy` class must provide the following static functions:

- `pressure` : Compute pressure.
- `zfactorCubicEq` : Compute coefficients of the cubic equation of Z-factor.
- `lnFugacityCoeff` : Compute the natural log of fugacity coefficient.
- `residualEnthalpy` : Compute residual enthalpy.
- `residualEntropy` : Compute residual entropy.
- `residualHelmholzEnergy` : Compute residual Helmholtz energy.

A `CorrectionFactor` class must provide the following member functions:

- `value` : Compute the correction factor.
- `derivative` : Compute the derivative of the correction factor.

Helper functions are defined for each EoS to easily create EoS instances:

- `eos::makeVanDerWaalsEos`
- `eos::makeSoaveRedlichKwongEos`
- `eos::makePengRobinsonEos`

## Dependencies

This library only depends on the STL library.

Unit tests use [GNU Scientific Library](https://www.gnu.org/software/gsl/) and [Googletest](https://github.com/google/googletest). Googletest is included as a git submodule under `third-party` directory.

## An Example of Usage

```cpp
#include <gsl/gsl_poly.h> // gsl_poly_solve_cubic

#include "eos/peng_robinson_eos.hpp"

// A solver for cubic equation using GSL library
struct CubicEquationSolver {
  std::vector<double> operator()(const std::array<double, 3>& a) {
    std::vector<double> x(3);
    const auto n = gsl_poly_solve_cubic(a[0], a[1], a[2], &x[0], &x[1], &x[2]);
    x.resize(n);
    return x;
  }
};

int main() {
  // Methane
  const double pc = 4.6e6;    // Critical pressure [Pa]
  const double tc = 190.6;    // Critical temperature [K]
  const double omega = 0.008; // Acentric factor
  
  // Create an instance of Peng-Robinson EoS
  const auto eos = eos::makePengRobinsonEos(pc, tc, omega);
  
  {
    const double p = 3e6;    // Pressure [Pa]
    const double t = 180.0;  // Temperature [K]

    // Compute z-factor at the pressure and temperature.
    // Please note that there can be multile values for Z-factor.
    // `params` contains parameters required for calculating thermodynamic properties,
    // e.g. fugacity coeffcients.
    const auto [z, params] = eos.zfactor(p, t, CubicEquationSolver{});

    // Compute the fugacity coefficient of a Z-factor.
    const auto phi = eos.fugacityCoeff(z[0], params);
  }

  {
    const double t = 180.0;  // Temperature [K]
    const double v = 0.001;  // Volume [m3]

    // Compute pressure at given temperature and volume.
    const auto p = eos.pressure(t, v);
  }
}
```