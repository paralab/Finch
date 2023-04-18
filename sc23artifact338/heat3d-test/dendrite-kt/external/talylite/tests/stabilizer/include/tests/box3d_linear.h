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

class Box3DLinearStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "3D Box (Linear)"; }
  ElemType elm_type() override { return kElem3dHexahedral; }
  GridType grid_type() override { return kGrid3dBox; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    double A0[8][3] = {};
    A0[0][0] = -7.0/5.0;
    A0[0][1] = -4.0/5.0;
    A0[0][2] = -1.0/5.0;
    A0[1][0] = 2.0/5.0;
    A0[1][1] = -4.0/5.0;
    A0[1][2] = -1.0/5.0;
    A0[2][0] = 2.0/5.0;
    A0[2][1] = 3.0/5.0;
    A0[2][2] = -1.0/5.0;
    A0[3][0] = -7.0/5.0;
    A0[3][1] = 3.0/5.0;
    A0[3][2] = -1.0/5.0;
    A0[4][0] = -7.0/5.0;
    A0[4][1] = -4.0/5.0;
    A0[4][2] = 2.0/5.0;
    A0[5][0] = 2.0/5.0;
    A0[5][1] = -4.0/5.0;
    A0[5][2] = 2.0/5.0;
    A0[6][0] = 2.0/5.0;
    A0[6][1] = 3.0/5.0;
    A0[6][2] = 2.0/5.0;
    A0[7][0] = -7.0/5.0;
    A0[7][1] = 3.0/5.0;
    A0[7][2] = 2.0/5.0;

    assert(i >= 0 && i < 8);
    return ZEROPTV(A0[i][0], A0[i][1], A0[i][2]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[8][8] = {};
    double t2 = sqrt(3.0);
    double t3 = t2*(1.0/1.2E1);
    double t4 = t3+1.0/6.0;
    double t5 = -t3+1.0/6.0;
    double t6 = t3-1.0/6.0;
    double t7 = -t3-1.0/6.0;
    A0[0][0] = t2*(-1.0/1.2E1)-1.0/6.0;
    A0[0][1] = t4;
    A0[0][2] = 1.0/1.2E1;
    A0[0][3] = -1.0/1.2E1;
    A0[0][4] = -1.0/1.2E1;
    A0[0][5] = 1.0/1.2E1;
    A0[0][6] = t5;
    A0[0][7] = t6;
    A0[1][0] = t7;
    A0[1][1] = t4;
    A0[1][2] = 1.0/1.2E1;
    A0[1][3] = -1.0/1.2E1;
    A0[1][4] = -1.0/1.2E1;
    A0[1][5] = 1.0/1.2E1;
    A0[1][6] = t5;
    A0[1][7] = t6;
    A0[2][0] = -1.0/1.2E1;
    A0[2][1] = 1.0/1.2E1;
    A0[2][2] = t4;
    A0[2][3] = t7;
    A0[2][4] = t6;
    A0[2][5] = t5;
    A0[2][6] = 1.0/1.2E1;
    A0[2][7] = -1.0/1.2E1;
    A0[3][0] = -1.0/1.2E1;
    A0[3][1] = 1.0/1.2E1;
    A0[3][2] = t4;
    A0[3][3] = t7;
    A0[3][4] = t6;
    A0[3][5] = t5;
    A0[3][6] = 1.0/1.2E1;
    A0[3][7] = -1.0/1.2E1;
    A0[4][0] = -1.0/1.2E1;
    A0[4][1] = 1.0/1.2E1;
    A0[4][2] = t5;
    A0[4][3] = t6;
    A0[4][4] = t7;
    A0[4][5] = t4;
    A0[4][6] = 1.0/1.2E1;
    A0[4][7] = -1.0/1.2E1;
    A0[5][0] = -1.0/1.2E1;
    A0[5][1] = 1.0/1.2E1;
    A0[5][2] = t5;
    A0[5][3] = t6;
    A0[5][4] = t7;
    A0[5][5] = t4;
    A0[5][6] = 1.0/1.2E1;
    A0[5][7] = -1.0/1.2E1;
    A0[6][0] = t6;
    A0[6][1] = t5;
    A0[6][2] = 1.0/1.2E1;
    A0[6][3] = -1.0/1.2E1;
    A0[6][4] = -1.0/1.2E1;
    A0[6][5] = 1.0/1.2E1;
    A0[6][6] = t4;
    A0[6][7] = t7;
    A0[7][0] = t6;
    A0[7][1] = t5;
    A0[7][2] = 1.0/1.2E1;
    A0[7][3] = -1.0/1.2E1;
    A0[7][4] = -1.0/1.2E1;
    A0[7][5] = 1.0/1.2E1;
    A0[7][6] = t4;
    A0[7][7] = t7;

    assert(itg_pt >= 0 && itg_pt < 8);
    assert(bf >= 0 && bf < 8);
    return A0[itg_pt][bf];
  }
};

