#pragma once

#include <cstddef>

namespace eos {

/// @{
/// @name Physical constants
///
/// The following definition is based on the International System of Units(SI)
/// 9th edition 2019.

/// Gas constant [J/mol-K]
static constexpr double gas_constant = 8.314'462'618'153'24;

/// Avogadro constant [1/mol]
static constexpr double avogadro_constant = 6.022'140'76e23;

/// Boltzman constant [J/K]
static constexpr double boltzman_constant = 1.380'649e-23;

/// Planck constant [J-s]
static constexpr double planck_constant = 6.626'070'15e-34;

/// @}

/// @{
/// @name Other constants

/// Size of vectors
static constexpr std::size_t dynamic_extent = SIZE_MAX;

/// @}

}  // namespace eos