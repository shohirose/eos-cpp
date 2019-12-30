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
const double omega = 0.008; // Accentric factor

// Creates EoS object
const auto eos = pr_eos<double, 1>(pc, tc, omega);

// Computes z-factor and fugacity coefficient
{
  const double p = 3e6;
  const double t = 180.0;

  // Creates a state at a given pressure and temperature.
  const auto state = eos.state(p, t);

  // Computes z-factor and fugacity coefficient.
  // Please note that there can be multile values for z-factor.
  const auto z = state.zfactor();
  const auto phi = state.fugacity_coeff(z[0]);
}

// PVT line calculation
{
  // Creates PVT line
  const auto pvt = eos.pvt();
  // Volume [m3]
  const auto v = 0.001;
  // Computes pressure at a given temperature and volume
  const auto p = pvt.pressure(t, v);
}
```