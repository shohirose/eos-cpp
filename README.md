# eoscpp

This is c++ code for equation of states (EoSs).

This module provides templated classes for three types of cubic EoSs: van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- vdw_eos
- srk_eos
- pr_eos

## Example of Usage

First, let's create an eos:

```cpp
// Critical parameters of methane
const double pc = 4e6;      // Critical pressure [Pa]
const double tc = 190.6;    // Critical temperature [K]
const double omega = 0.008; // Acentric factor

// Creates EOS
const auto eos = make_pr_eos(pc, tc, omega);
```

Z-factor and fugacity coefficients at given pressure and temperature can be computed from a state:

```cpp
const double p = 3e6;    // Pressure [Pa]
const double t = 180.0;  // Temperature [K]

// Creates a state at given pressure and temperature  
const auto state = eos.state(p, t);

// Computes z-factor at the pressure and temperature
// Please note that there can be multile values for z-factor.
const auto z = state.zfactor();

// Computes fugacity coefficient from a corresponding z-factor
const auto phi = state.fugacity_coeff(z[0]);
```

Pressure at given temperature and volume can be computed:

```cpp
const double t = 180.0;  // Temperature [K]
const double v = 0.001;  // Volume [m3]
const auto p = eos.pressure(t, v);
```

Pressure along a given temperature can be computed from an isothermal state:

```cpp
const double t = 180.0; // Temperature [K]

// Creates isothermal state
const auto state = eos.state(t);

const std::size_t n = 100; // Number of samples
std::vector<double> v(n); // Array of volume [m3]
std::vector<double> p(n); // Array of pressure [Pa]

// Initialize arrays ...

// Computes presure along an isothermal line
for (std::size_t i = 0; i < n; ++i) {
  p[i] = state.pressure(v[i]);
}
```

Vapor pressure can be computed by flash calculation:

```cpp
// Estimate initial pressure for flash calculation by using Wilson equation
const auto p_init = estimate_vapor_pressure(t, pc, tc, omega);

// Creates flash object
const auto flash = make_flash(make_pr_eos(pc, tc, omega));

// Computes vapor pressure at a given temperature
const auto result = flash.vapor_pressure(p_init, t);

const auto p_vap = result.first;   // Vapor pressure
const auto report = result.second; // Iteration report
```