# eoscpp

This is c++ code for equation of states (EoSs).

This module provides templated classes for three types of cubic EoSs: van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- vdw_eos
- srk_eos
- pr_eos

## Example of Usage

```cpp
// Critical parameters
const double pc = 4e6;      // Critical pressure [Pa]
const double tc = 190.6;    // Critical temperature [K]
const double omega = 0.008; // Acentric factor

// Creates EoS object
const auto eos = pr_eos<double>(pc, tc, omega);

// Computes z-factor and fugacity coefficient
{
  const double p = 3e6;    // Pressure [Pa]
  const double t = 180.0;  // Temperature [K]

  // Computes z-factor and fugacity coefficient.
  // Please note that there can be multile values for z-factor.
  const auto z = eos.zfactor(p, t);
  const auto phi = eos.fugacity_coeff(z[0]);
}

// Computes pressure at given temperature and volume
{
  const double t = 180.0;  // Temperature [K]
  const double v = 0.001;  // Volume [m3]
  const auto p = eos.pressure(t, v);
}
```