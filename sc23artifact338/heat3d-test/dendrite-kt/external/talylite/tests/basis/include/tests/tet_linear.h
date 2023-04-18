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

#include <basistest.h>

// Yes, this should be split into a .cpp/.h pair.
// But this is only ever #included in one file (main.cpp) and I'm lazy.

// This test does not have analytic values because we don't know the
// analytic values of the tetrahedron gauss points.

class Tet3DLinearBasisTest : public BasisTest
{
public:
  std::string name() override { return "3D Tetrahedron (Linear)"; }
  ElemType elm_type() override { return kElem3dTetrahedral; }
  GridType grid_type() override { return kGrid3dTet; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-3.88, 2.06, -1.1),
      ZEROPTV(-1.47, -0.83, -0.7),
      ZEROPTV(0.51, -1.4, 1.3),
      ZEROPTV(-0.95, 0.27, 4.8)
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(0.58541020, 0.13819660, 0.13819660),
      ZEROPTV(0.13819660, 0.58541020, 0.13819660),
      ZEROPTV(0.13819660, 0.13819660, 0.58541020),
      ZEROPTV(0.13819660, 0.13819660, 0.13819660)
    };
    return pts[i];
  }

  double weight(int i) override {
    return 1.0/4.0 / 6.0;
  }

  double N(int itg_pt, int bf) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();
    double zeta_Var = pt.z();

    const double n[4] = {
      -zeta_Var-eta-xi+1.0,
      xi,
      eta,
      zeta_Var
    };

    return n[bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[4] = {
      ZEROPTV(-1, -1, -1),
      ZEROPTV(1, 0, 0),
      ZEROPTV(0, 1, 0),
      ZEROPTV(0, 0, 1),
    };
    return dnde[bf](axis);
  }

  double dXde(int itg_pt, int i, int j) override {
    static const double dxde[3][3] = {
      {2.41, 4.39, 2.93},
      {-2.89, -3.46, -1.79},
      {4.0E-1, 2.4, 5.9}
    };
    return dxde[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    static const double cof_values[3][3] = {
      {-1.6118E1, 1.6335E1, -6.94E2/1.25E2},
      {-1.8869E1, 1.3047E1, -4.028},
      {2.2797, -4.1538, 4.3485}
    };
    return cof_values[i][j];
  }

  double jacobian_det(int itg_pt) override {
    return 16.59891;
  }

  double dN(int itg_pt, int i, int axis) override {
    double dn[4][3] = {};
    dn[0][0] = 3.214066465810104E-1;
    dn[0][1] = 5.934124590108628E-1;
    dn[0][2] = -1.490700292971044E-1;
    dn[1][0] = -9.710276156687397E-1;
    dn[1][1] = -1.136761389753906;
    dn[1][2] = 1.373403434321892E-1;
    dn[2][0] = 9.841007632428876E-1;
    dn[2][1] = 7.860154672806829E-1;
    dn[2][2] = -2.502453474354641E-1;
    dn[3][0] = -3.344797941551584E-1;
    dn[3][1] = -2.426665365376401E-1;
    dn[3][2] = 2.619750333003794E-1;

    return dn[i][axis];
  }
};

