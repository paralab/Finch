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

class Box2DQuadraticStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "2D Box (Quadratic)"; }
  ElemType elm_type() override { return kElem2dBox; }
  GridType grid_type() override { return kGrid2dBox; }
  kBasisFunction basis_function() override { return BASIS_QUADRATIC; }
  int nsd() override { return 2; }

  ZEROPTV node_position(int i) override {
    double A0[9][2] = {};
    A0[0][0] = -1.0;
    A0[1][0] = 1.0;
    A0[2][0] = 1.0;
    A0[2][1] = 2.0;
    A0[3][0] = -1.0;
    A0[3][1] = 2.0;
    A0[5][0] = 1.0;
    A0[5][1] = 1.0;
    A0[6][1] = 2.0;
    A0[7][0] = -1.0;
    A0[7][1] = 1.0;
    A0[8][1] = 1.0;

    assert(i >= 0 && i < 9);
    return ZEROPTV(A0[i][0], A0[i][1]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[9][9] = {};
    double t2 = sqrt(1.5E1);
    double t3 = t2*(1.0/4.4E1);
    double t4 = t3+9.0/4.4E1;
    double t5 = t2*(5.0/4.4E1);
    double t6 = t5-2.1E1/4.4E1;
    double t7 = t2*(1.0/1.1E1);
    double t8 = t2*(1.0/1.76E2);
    double t9 = -t3-9.0/4.4E1;
    double t10 = -t5+2.1E1/4.4E1;
    double t11 = -t7+2.0/1.1E1;
    double t12 = t2*(1.0/2.4E1);
    double t13 = t2*(3.0/1.76E2);
    double t14 = t2*(1.7E1/1.76E2);
    double t15 = t8-1.3E1/1.76E2;
    double t16 = t2*(5.0/1.76E2);
    double t17 = t2*(2.0/3.3E1);
    double t18 = t2*(1.0/3.3E1);
    double t19 = t7-2.0/1.1E1;
    double t20 = t2+2.0;
    double t21 = 1.0/t20;
    double t22 = t14-6.7E1/1.76E2;
    double t23 = t13-1.7E1/1.76E2;
    double t24 = t16+2.3E1/1.76E2;
    double t25 = -t8+1.3E1/1.76E2;
    double t26 = t18+1.0/4.4E1;
    double t27 = t17-9.0/4.4E1;
    A0[0][0] = t2*(-5.0/1.76E2)-2.3E1/1.76E2;
    A0[0][1] = t15;
    A0[0][2] = t2*(-1.7E1/1.76E2)+6.7E1/1.76E2;
    A0[0][3] = t2*(-3.0/1.76E2)+1.7E1/1.76E2;
    A0[0][4] = t4;
    A0[0][5] = t2*(-2.0/3.3E1)+9.0/4.4E1;
    A0[0][6] = t6;
    A0[0][7] = t2*(-1.0/3.3E1)-1.0/4.4E1;
    A0[0][8] = t19;
    A0[1][0] = t9;
    A0[1][1] = t4;
    A0[1][2] = t6;
    A0[1][3] = t10;
    A0[1][5] = t21;
    A0[1][7] = t11;
    A0[2][0] = t25;
    A0[2][1] = t24;
    A0[2][2] = t23;
    A0[2][3] = t22;
    A0[2][4] = t9;
    A0[2][5] = t26;
    A0[2][6] = t10;
    A0[2][7] = t27;
    A0[2][8] = t11;
    A0[3][5] = t12-1.0/4.0;
    A0[3][7] = -t12-1.0/4.0;
    A0[3][8] = 1.0/2.0;
    A0[4][5] = 1.0/2.0;
    A0[4][7] = -1.0/2.0;
    A0[5][5] = t12+1.0/4.0;
    A0[5][7] = -t12+1.0/4.0;
    A0[5][8] = -1.0/2.0;
    A0[6][0] = -t13+1.7E1/1.76E2;
    A0[6][1] = -t14+6.7E1/1.76E2;
    A0[6][2] = t15;
    A0[6][3] = -t16-2.3E1/1.76E2;
    A0[6][4] = t6;
    A0[6][5] = -t17+9.0/4.4E1;
    A0[6][6] = t4;
    A0[6][7] = -t18-1.0/4.4E1;
    A0[6][8] = t19;
    A0[7][0] = t10;
    A0[7][1] = t6;
    A0[7][2] = t4;
    A0[7][3] = t9;
    A0[7][5] = t21;
    A0[7][7] = t11;
    A0[8][0] = t22;
    A0[8][1] = t23;
    A0[8][2] = t24;
    A0[8][3] = t25;
    A0[8][4] = t10;
    A0[8][5] = t26;
    A0[8][6] = t9;
    A0[8][7] = t27;
    A0[8][8] = t11;

    assert(itg_pt >= 0 && itg_pt < 9);
    assert(bf >= 0 && bf < 9);
    return A0[itg_pt][bf];
  }
};

