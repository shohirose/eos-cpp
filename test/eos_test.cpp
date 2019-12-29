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

#include "shirose/eos.hpp"
#include <gtest/gtest.h>
#include <iostream>

using namespace shirose;

TEST(EosTest, VanDerWaalsEosTest) {
  const double pc = 4e6;
  const double tc = 190.6;
  const double omega = 0.008;

  auto eos = vdw_eos<double, 1>(pc, tc, omega);
  const double p = 3e6;
  const double t = 180.0;
  eos.set_pt(p, t);
  const auto z = eos.zfactor();

  for (size_t i = 0; i < z.size(); ++i) {
    std::cout << "z = " << z[i] << ' ' << "phi = " << eos.fugacity_coeff(z[i])
              << '\n';
  }
}