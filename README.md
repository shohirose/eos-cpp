# eoscpp

This is c++ code for equation of states (EoSs).

This module provides templated classes for three types of cubic EoS: van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- `eos::van_der_waals_eos`
- `eos::soave_redlich_kwong_eos`
- `eos::peng_robinson_eos`

In addition, you can customize and create a new cubic EoS by deriving from `cubic_eos_base`. `cubic_eos_base` is a template class for general two-parameter cubic EoS. `cubic_eos_base` requires a policy class which satisfies the following conditions:

- Declares the following types:
    - `derived_type`
    - `value_type`
- Defines the following constants:
    - `omega_a`
    - `omega_b`
- Defines the following static functions:
    - `pressure`
    - `cubic_eq`
    - `fugacity_coeff`
    - `residual_enthalpy`
    - `residual_entropy`

An example of defining a custom EoS class:

```cpp
// Forward declaration of my custom EoS class.
template <typename T>
class my_cubic_eos;

// Policy class for my custom EoS class.
template <typename T>
struct my_cubic_eos_policy {
  using derived_type = my_cubic_eos;
  using value_type = T;

  static constexpr double omega_a = // ...
  static constexpr double omega_b = // ...

  static T pressure(const T& a, const T& b) noexcept {
    // ...
  }

  static std::array<T, 3> cubic_eq(cosnt T& a, const T& b) noexcept {
    // ...
  }

  // ...
};

template <typename T>
class my_cubic_eos : public cubic_eos_base<my_cubic_eos_policy<T>> {
 public:
   using base_type = cubic_eos_base<my_cubic_eos_policy<T>>;

   my_cubic_eos() = default;

   my_cubic_eos(const T& pc, const T& tc, /* ... */) noexcept
     : base_type{pc, tc},
       // ....
     {}

  void set_params(const T& pc, const T& tc, /* ... */) noexcept {
    this->base_type::set_params(pc, tc);
    // ...
  }
  
  T alpha(const T& tr) const noexcept {
    // ...
  }

  T beta(const T& tr) const noexcept {
    // ...
  }
};

```

Helper functions are defined for each EoS to easily create EoS objects:

- `eos::make_vdw_eos`
- `eos::make_sdk_eos`
- `eos::make_pr_eos`

## Example of Usage

First, let's create an EoS:

```cpp
// Critical parameters of methane
const double pc = 4e6;      // Critical pressure [Pa]
const double tc = 190.6;    // Critical temperature [K]
const double omega = 0.008; // Acentric factor

// Creates EoS
const auto eos = eos::make_pr_eos(pc, tc, omega);
```

Z-factor and fugacity coefficients at given pressure and temperature can be computed from a state:

```cpp
const double p = 3e6;    // Pressure [Pa]
const double t = 180.0;  // Temperature [K]

// Creates a state at given pressure and temperature  
const auto state = eos.state(p, t);

// Computes z-factor at the pressure and temperature
// Please note that there can be multile values for z-factor.
const std::vector<double> z = state.zfactor();

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
const auto p_init = eos::estimate_vapor_pressure(t, pc, tc, omega);

// Creates flash object
const auto flash = eos::make_flash(make_pr_eos(pc, tc, omega));

// Computes vapor pressure at a given temperature
auto [p_vap, report] = flash.vapor_pressure(p_init, t);
```