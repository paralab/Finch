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

class Tri2DLinearStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "2D Triangle (Linear)"; }
  ElemType elm_type() override { return kElem2dTriangle; }
  GridType grid_type() override { return kGrid2dTriangle; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int nsd() override { return 2; }

  ZEROPTV node_position(int i) override {
    double A0[3][2] = {};
    A0[0][0] = -3.0;
    A0[0][1] = -3.0;
    A0[1][0] = 2.0;
    A0[2][0] = 1.0;
    A0[2][1] = 2.0;

    assert(i >= 0 && i < 3);
    return ZEROPTV(A0[i][0], A0[i][1]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[3][3] = {};
    A0[0][0] = -1.0/5.0;
    A0[0][1] = 1.0/2.0;
    A0[0][2] = -3.0/1.0E1;
    A0[1][0] = -1.0/5.0;
    A0[1][1] = 1.0/2.0;
    A0[1][2] = -3.0/1.0E1;
    A0[2][0] = -1.0/5.0;
    A0[2][1] = 1.0/2.0;
    A0[2][2] = -3.0/1.0E1;

    assert(itg_pt >= 0 && itg_pt < 3);
    assert(bf >= 0 && bf < 3);
    return A0[itg_pt][bf];
  }
};

