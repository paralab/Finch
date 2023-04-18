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

class Tet3DLinearStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "3D Tetrahedron (Linear)"; }
  ElemType elm_type() override { return kElem3dTetrahedral; }
  GridType grid_type() override { return kGrid3dTet; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    double A0[4][3] = {};
    A0[0][0] = -3.88;
    A0[0][1] = 2.06;
    A0[0][2] = -1.1;
    A0[1][0] = -1.47;
    A0[1][1] = -8.3E-1;
    A0[1][2] = -7.0E-1;
    A0[2][0] = 5.1E-1;
    A0[2][1] = -1.4;
    A0[2][2] = 1.3;
    A0[3][0] = -9.5E-1;
    A0[3][1] = 2.7E-1;
    A0[3][2] = 4.8;

    assert(i >= 0 && i < 4);
    return ZEROPTV(A0[i][0], A0[i][1], A0[i][2]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[4][4] = {};
    A0[0][0] = 1.230964467005076E-1;
    A0[0][1] = -3.718966312874942E-1;
    A0[0][2] = 3.769035532994924E-1;
    A0[0][3] = -1.281033687125058E-1;
    A0[1][0] = 1.230964467005076E-1;
    A0[1][1] = -3.718966312874942E-1;
    A0[1][2] = 3.769035532994924E-1;
    A0[1][3] = -1.281033687125058E-1;
    A0[2][0] = 1.230964467005076E-1;
    A0[2][1] = -3.718966312874942E-1;
    A0[2][2] = 3.769035532994924E-1;
    A0[2][3] = -1.281033687125058E-1;
    A0[3][0] = 1.230964467005076E-1;
    A0[3][1] = -3.718966312874942E-1;
    A0[3][2] = 3.769035532994924E-1;
    A0[3][3] = -1.281033687125058E-1;

    assert(itg_pt >= 0 && itg_pt < 4);
    assert(bf >= 0 && bf < 4);
    return A0[itg_pt][bf];
  }
};

