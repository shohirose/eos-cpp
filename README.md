# eoscpp

This is c++ code for equation of states (EoSs).

This module provides templated classes for three types of cubic EoS: van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- `eos::van_der_waals_eos`
- `eos::soave_redlich_kwong_eos`
- `eos::peng_robinson_eos`

In addition, you can customize and create a new cubic EoS by deriving from `cubic_eos_base`. `cubic_eos_base` is a template class for general two-parameter cubic EoS. `cubic_eos_base` requires a derived class which satisfies the following conditions:

- Defines the following static functions:
    - `pressure`
    - `zfactor_cubic_eq`
    - `fugacity_coeff`
    - `residual_enthalpy`
    - `residual_entropy`
- Defines the following member functions
    - `alpha`
    - `beta`

An example of defining a custom EoS class:

```cpp
class my_cubic_eos : public cubic_eos_base<my_cubic_eos> {
 public:
  using base_type = cubic_eos_base<my_cubic_eos>;

  static double pressure(double a, double b) noexcept;
  static std::array<double, 3> zfactor_cubic_eq(double a, double b) noexcept;
  static double fugacity_coeff(double z, double a, double b) noexcept;
  static double residual_enthalpy(double z, double t, double a, double b, double beta) noexcept;
  static double residual_entropy(double z, double a, double b, double beta) noexcept;

  my_cubic_eos() = default;
  my_cubic_eos(double pc, double tc, /* ... */) noexcept;
  void set_params(double pc, double tc, /* ... */) noexcept;
  double alpha(double tr) const noexcept;
  double beta(double tr) const noexcept;
};

```

In addition to the custom eos class, a custom `cubic_eos_traits` class must be defined. The trait class must define the following constants:

- `omega_a`
- `omega_b`

Helper functions are defined for each EoS to easily create EoS objects:

- `eos::make_van_der_waals_eos`
- `eos::make_soave_redlich_kwong_eos`
- `eos::make_peng_robinson_eos`

## Dependencies

This library depends on the [GNU Scientific Library](https://www.gnu.org/software/gsl/). The [Googletest](https://github.com/google/googletest) is used for unit testing but it is included as a git submodule under the `third-party` directory. Please make sure to run `git submodule update`.

## Example of Usage

First, let's create an EoS:

```cpp
// Critical parameters of methane
const double pc = 4e6;      // Critical pressure [Pa]
const double tc = 190.6;    // Critical temperature [K]
const double omega = 0.008; // Acentric factor

// Creates EoS
const auto eos = eos::make_peng_robinson_eos(pc, tc, omega);
```

Z-factor and fugacity coefficients at given pressure and temperature can be computed from a state:

```cpp
const double p = 3e6;    // Pressure [Pa]
const double t = 180.0;  // Temperature [K]

// Creates a state at given pressure and temperature  
const auto state = eos.create_isobaric_isothermal_state(p, t);

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

const auto state = eos.create_isothermal_state(t);

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
const auto p_init = eos::estimate_vapor_pressure(t, pc, tc, omega);

const auto eos = eos::make_peng_robinson_eos(pc, tc, omega);
const auto flash = eos::make_vapor_liquid_flash(eos);

// Computes vapor pressure at a given temperature
const auto [p_vap, result] = flash.vapor_pressure(p_init, t);
```