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

class Box1DCubicStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "1D Box (Cubic)"; }
  ElemType elm_type() override { return kElem1d; }
  GridType grid_type() override { return kGrid1d; }
  kBasisFunction basis_function() override { return BASIS_CUBIC; }
  int nsd() override { return 1; }

  ZEROPTV node_position(int i) override {
    double A0[4][1] = {};
    A0[0][0] = -1.0;
    A0[1][0] = 1.0;
    A0[2][0] = -1.0/3.0;
    A0[3][0] = 1.0/3.0;

    assert(i >= 0 && i < 4);
    return ZEROPTV(A0[i][0]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[4][4] = {};
    double t2 = sqrt(5.0);
    double t3 = sqrt(6.0);
    double t4 = t2*t3*(2.7E1/3.5E1);
    double t5 = t4+1.6E1/7.0;
    double t6 = 1.0/t5;
    double t7 = t2*t3*(2.7E1/2.8E2);
    double t8 = t2*t3*(2.0/3.5E1);
    double t9 = t8+3.0/7.0;
    double t10 = sqrt(t9);
    double t11 = t10*(9.0/8.0);
    double t12 = t2*t3*(8.1E1/2.8E2);
    double t13 = sqrt(3.0E1);
    double t17 = t13*2.0;
    double t14 = -t17+1.5E1;
    double t15 = sqrt(t14);
    double t16 = sqrt(3.5E1);
    double t18 = sqrt(4.2E1);
    double t19 = t13*3.32843857072932E-2;
    double t20 = t15*t16*5.594014404587092E-4;
    double t21 = t15*t16*(3.0/4.54E2);
    double t22 = t15*t18*(1.1E1/4.54E2);
    double t23 = t13*(3.0/2.27E2);
    double t24 = t7-t11+3.7E1/5.6E1;
    double t25 = t6*t24;
    double t26 = t7+t11+3.7E1/5.6E1;
    double t27 = -t11+t12+2.7E1/5.6E1;
    double t28 = t11+t12+2.7E1/5.6E1;
    double t29 = t6*t28;
    A0[0][0] = -t6*t26;
    A0[0][1] = t25;
    A0[0][2] = t29;
    A0[0][3] = -t6*t27;
    A0[1][0] = t13*(-4.65002447381302E-2)+t20-t15*t18*1.382187725800061E-2+2.547723935389134E-1;
    A0[1][1] = t19-t15*t16*7.167330955877211E-3-t15*t18*1.040719763186724E-2-3.1E1/6.81E2;
    A0[1][2] = t21+t22+t23-1.61E2/2.27E2;
    A0[1][3] = 1.0/2.0;
    A0[2][0] = -t19+t15*t16*7.167330955877211E-3+t15*t18*1.040719763186724E-2+3.1E1/6.81E2;
    A0[2][1] = t13*4.65002447381302E-2-t20+t15*t18*1.382187725800061E-2-2.547723935389134E-1;
    A0[2][2] = -1.0/2.0;
    A0[2][3] = -t21-t22-t23+1.61E2/2.27E2;
    A0[3][0] = -t25;
    A0[3][1] = t6*t26;
    A0[3][2] = t6*t27;
    A0[3][3] = -t29;

    assert(itg_pt >= 0 && itg_pt < 4);
    assert(bf >= 0 && bf < 4);
    return A0[itg_pt][bf];
  }
};

