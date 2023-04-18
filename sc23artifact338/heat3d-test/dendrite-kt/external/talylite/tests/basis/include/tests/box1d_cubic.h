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

class Box1DCubicBasisTest : public BasisTest
{
public:
  std::string name() override { return "1D (Cubic)"; }
  ElemType elm_type() override { return kElem1d; }
  GridType grid_type() override { return kGrid1d; }
  kBasisFunction basis_function() override { return BASIS_CUBIC; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 1; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-7.0/5.0),
      ZEROPTV(2.0/5.0),
      ZEROPTV(-4.0/5.0),
      ZEROPTV(-1.0/5.0),
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(sqrt(35) * sqrt(sqrt(30) * 2 + 15) * (-1.0 / 35)),
      ZEROPTV(sqrt(35) * sqrt(sqrt(30) * -2 + 15) * (-1.0 / 35)),
      ZEROPTV(sqrt(35) * sqrt(sqrt(30) * -2 + 15) * (1.0 / 35)),
      ZEROPTV(sqrt(35) * sqrt(sqrt(30) * 2.0 + 15) * (1.0 / 35))
    };
    return pts[i];
  }

  double weight(int i) override {
    static const double wts[4] = {
      sqrt(30) * (-1.0 / 36) + 0.5,
      sqrt(30) * (1.0 / 36) + 0.5,
      sqrt(30) * (1.0 / 36) + 0.5,
      sqrt(30) * (-1.0 / 36) + 0.5
    };
    return wts[i];
  }

  double N(int itg_pt, int bf) override {
    double xi = gauss_point(itg_pt).x();

    double t2 = xi*xi;
    double t3 = t2*(9.0/1.6E1);
    double t4 = t2*xi*(2.7E1/1.6E1);
    double t5 = xi*(1.0/1.6E1);
    double n[4] = {
      n[0] = t3+t5-t2*xi*(9.0/1.6E1)-1.0/1.6E1,
      n[1] = t3-t5+t2*xi*(9.0/1.6E1)-1.0/1.6E1,
      n[2] = -t3+t4-xi*(2.7E1/1.6E1)+9.0/1.6E1,
      n[3] = -t3-t4+xi*(2.7E1/1.6E1)+9.0/1.6E1,
    };

    return n[bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    double xi = gauss_point(itg_pt).x();
    double t2 = xi*(9.0/8.0);
    double t3 = xi*xi;
    double t4 = t3*(8.1E1/1.6E1);

    double dnde[4] = {
      t2-t3*(2.7E1/1.6E1)+1.0/1.6E1,
      t2+t3*(2.7E1/1.6E1)-1.0/1.6E1,
      -t2+t4-2.7E1/1.6E1,
      -t2-t4+2.7E1/1.6E1,
    };

    return dnde[bf];
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
    double xi = gauss_point(itg_pt).x();
    double t2 = xi*(5.0/4.0);
    double t3 = xi*xi;
    double t4 = t3*(4.5E1/8.0);

    double dn[4] = {
      t2-t3*(1.5E1/8.0)+5.0/7.2E1,
      t2+t3*(1.5E1/8.0)-5.0/7.2E1,
      -t2+t4-1.5E1/8.0,
      -t2-t4+1.5E1/8.0,
    };

    return dn[i];
  }
};

