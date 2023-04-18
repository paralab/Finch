/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#pragma once

#include <stabilizertest.h>

class Box1DLinearStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "1D Box (Linear)"; }
  ElemType elm_type() override { return kElem1d; }
  GridType grid_type() override { return kGrid1d; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int nsd() override { return 1; }

  ZEROPTV node_position(int i) override {
    double A0[2][1] = {};
    A0[0][0] = -7.0/5.0;
    A0[1][0] = 2.0/5.0;

    assert(i >= 0 && i < 2);
    return ZEROPTV(A0[i][0]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[2][2] = {};
    A0[0][0] = -1.0/2.0;
    A0[0][1] = 1.0/2.0;
    A0[1][0] = -1.0/2.0;
    A0[1][1] = 1.0/2.0;

    assert(itg_pt >= 0 && itg_pt < 2);
    assert(bf >= 0 && bf < 2);
    return A0[itg_pt][bf];
  }
};

