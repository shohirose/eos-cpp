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