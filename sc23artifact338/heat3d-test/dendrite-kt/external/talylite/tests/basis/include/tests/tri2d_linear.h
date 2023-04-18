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

class Tri2DLinearBasisTest : public BasisTest
{
public:
  std::string name() override { return "2D Triangle (Linear)"; }
  ElemType elm_type() override { return kElem2dTriangle; }
  GridType grid_type() override { return kGrid2dTriangle; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 2; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[3] = {
      ZEROPTV(-3, -3),
      ZEROPTV(2, 0),
      ZEROPTV(1, 2),
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[3] = {
      ZEROPTV(0.5, 0),
      ZEROPTV(0, 0.5),
      ZEROPTV(0.5, 0.5),
    };
    return pts[i];
  }

  double weight(int i) override {
    return 1.0/3.0 / 2.0;
  }

  double N(int itg_pt, int bf) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    const double n[3] = {
      -eta-xi+1.0,
      xi,
      eta,
    };

    return n[bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[3] = {
      ZEROPTV(-1, -1),
      ZEROPTV(1, 0),
      ZEROPTV(0, 1),
    };
    return dnde[bf](axis);
  }

  double dXde(int itg_pt, int i, int j) override {
    static const double dxde[2][2] = {
      { 5, 4 },
      { 3, 5 },
    };
    return dxde[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    static const double cof_values[4] = { 5, -4, -3, 5 };
    return cof_values[2*j+i];
  }

  double jacobian_det(int itg_pt) override {
    return 13.0;
  }

  double dN(int itg_pt, int i, int axis) override {
    double dn[3][2];
    dn[0][0] = -2.0/1.3E1;
    dn[0][1] = -1.0/1.3E1;
    dn[1][0] = 5.0/1.3E1;
    dn[1][1] = -4.0/1.3E1;
    dn[2][0] = -3.0/1.3E1;
    dn[2][1] = 5.0/1.3E1;

    return dn[i][axis];
  }
};

