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

class Box1DLinearBasisTest : public BasisTest
{
public:
  std::string name() override { return "1D (Linear)"; }
  ElemType elm_type() override { return kElem1d; }
  GridType grid_type() override { return kGrid1d; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 1; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[2] = {
      ZEROPTV(-7.0/5.0),
      ZEROPTV(2.0/5.0)
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[2] = {
      ZEROPTV(-sqrt(3)/3.0),
      ZEROPTV(sqrt(3)/3.0)
    };
    return pts[i];
  }

  double weight(int i) override {
    return 1;
  }

  double N(int itg_pt, int bf) override {
    static const double n[2][2] = {
      { pow(3.0, 1.0/2) / 6.0 + 1.0/2,    1.0/2 - pow(3.0, 1.0/2) / 6.0 },
      { 1.0/2 - pow(3, 1.0/2)/6.0,        pow(3.0, 1.0/2) / 6.0 + 1.0/2 }
    };
    return n[itg_pt][bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[2][2] = {
      { ZEROPTV(-1.0/2), ZEROPTV(1.0/2) },
      { ZEROPTV(-1.0/2), ZEROPTV(1.0/2) }
    };
    return dnde[itg_pt][bf](axis);
  }

  double dXde(int itg_pt, int i, int j) override {
    return 9.0/10;
  }

  double cof(int itg_pt, int i, int j) override {
    return 1;
  }

  double jacobian_det(int itg_pt) override {
    return 9.0/10;
  }

  double dN(int itg_pt, int i, int axis) override {
    static const double dn[2][2][1] = {
      { {-5.0/9}, {5.0/9} },
      { {-5.0/9}, {5.0/9} }
    };
    return dn[itg_pt][i][axis];
  }
};

