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

class Box2DLinearBasisTest : public BasisTest
{
public:
  std::string name() override { return "2D Box (Linear)"; }
  ElemType elm_type() override { return kElem2dBox; }
  GridType grid_type() override { return kGrid2dBox; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 2; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-3.88, -1.4),
      ZEROPTV(0.51, -1.4),
      ZEROPTV(0.51, 0.27),
      ZEROPTV(-3.88, 0.27)
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-1 / sqrt(3), -1/sqrt(3)),
      ZEROPTV(1/sqrt(3), -1/sqrt(3)),
      ZEROPTV(-1/sqrt(3), 1/sqrt(3)),
      ZEROPTV(1/sqrt(3), 1/sqrt(3))
    };
    return pts[i];
  }

  double weight(int i) override {
    return 1;
  }

  double N(int itg_pt, int bf) override {
    static const double n[4][4] = {
      { (2.0 + sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 - sqrt(3)) / 6.0,  1.0 / 6.0 },
      { 1.0 / 6.0,  (2.0 + sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 - sqrt(3)) / 6.0 },
      { 1.0 / 6.0,  (2.0 - sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 + sqrt(3)) / 6.0 },
      { (2.0 - sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 + sqrt(3)) / 6.0,  1.0 / 6.0 },
    };
    return n[itg_pt][bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[4][4] = {
      { ZEROPTV((-3 - sqrt(3)) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (-3 + sqrt(3)) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0),
        ZEROPTV((sqrt(3) - 3) / 12.0, (3 + sqrt(3)) / 12.0)
      },
      { ZEROPTV((-3 - sqrt(3)) / 12.0, (sqrt(3) - 3) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (3 + sqrt(3)) / 12.0),
        ZEROPTV((-3 + sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0)
      },
      { ZEROPTV((sqrt(3) - 3) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (sqrt(3) - 3) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0),
        ZEROPTV((-3 - sqrt(3)) / 12.0, (3 + sqrt(3)) / 12.0)
      },
      { ZEROPTV((sqrt(3) - 3) / 12.0, (sqrt(3) - 3) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (3 + sqrt(3)) / 12.0),
        ZEROPTV((-3 - sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0)
      }
    };
    return dnde[itg_pt][bf](axis);
  }

  double dXde(int itg_pt, int i, int j) override {
    static const double dxde[4][2][2] = {
      { { 2.195, 0 },
        { 0, 0.835 },
      },
      { { 2.195, 0 },
        { 0, 0.835 },
      },
      { { 2.195, 0 },
        { 0, 0.835 },
      },
      { { 2.195, 0 },
        { 0, 0.835 },
      }
    };
    return dxde[itg_pt][i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    static const double cof_values[4] = { 0.835, 0, 0, 2.195 };
    return cof_values[2*j+i];
  }

  double jacobian_det(int itg_pt) override {
    return 2.195 * 0.835;
  }

  double dN(int itg_pt, int i, int axis) override {
    static const double dn[4][4][2] = {
      {                            // GP0
        { -(3 + sqrt(3)) / (6 * 4.39),  -(3 + sqrt(3)) / (6 * 1.67) },  // dN(0, 0) and dN(0, 1)
        {  (3 + sqrt(3)) / (6 * 4.39),  -(3 - sqrt(3)) / (6 * 1.67) },  // dN(1, 0) and dN(1, 1)
        {  (3 - sqrt(3)) / (6 * 4.39),   (3 - sqrt(3)) / (6 * 1.67) },  // dN(2, 0) and dN(2, 1)
        { -(3 - sqrt(3)) / (6 * 4.39),   (3 + sqrt(3)) / (6 * 1.67) }   // dN(3, 0) and dN(3, 1)
      },
      {                            // GP1
        { -(3 + sqrt(3)) / (6 * 4.39),  -(3 - sqrt(3)) / (6 * 1.67) },  // dN(0, 0) and dN(0, 1)
        {  (3 + sqrt(3)) / (6 * 4.39),  -(3 + sqrt(3)) / (6 * 1.67) },  // dN(1, 0) and dN(1, 1)
        {  (3 - sqrt(3)) / (6 * 4.39),   (3 + sqrt(3)) / (6 * 1.67) },  // dN(2, 0) and dN(2, 1)
        { -(3 - sqrt(3)) / (6 * 4.39),   (3 - sqrt(3)) / (6 * 1.67) }   // dN(3, 0) and dN(3, 1)
      },
      {                            // GP2
        { -(3 - sqrt(3)) / (6 * 4.39),  -(3 + sqrt(3)) / (6 * 1.67) },  // dN(0, 0) and dN(0, 1)
        {  (3 - sqrt(3)) / (6 * 4.39),  -(3 - sqrt(3)) / (6 * 1.67) },  // dN(1, 0) and dN(1, 1)
        {  (3 + sqrt(3)) / (6 * 4.39),   (3 - sqrt(3)) / (6 * 1.67) },  // dN(2, 0) and dN(2, 1)
        { -(3 + sqrt(3)) / (6 * 4.39),   (3 + sqrt(3)) / (6 * 1.67) }   // dN(3, 0) and dN(3, 1)
      },
      {                            // GP3
        { -(3 - sqrt(3)) / (6 * 4.39),  -(3 - sqrt(3)) / (6 * 1.67) },  // dN(0, 0) and dN(0, 1)
        {  (3 - sqrt(3)) / (6 * 4.39),  -(3 + sqrt(3)) / (6 * 1.67) },  // dN(1, 0) and dN(1, 1)
        {  (3 + sqrt(3)) / (6 * 4.39),   (3 + sqrt(3)) / (6 * 1.67) },  // dN(2, 0) and dN(2, 1)
        { -(3 + sqrt(3)) / (6 * 4.39),   (3 - sqrt(3)) / (6 * 1.67) }   // dN(3, 0) and dN(3, 1)
      }
    };
    return dn[itg_pt][i][axis];
  }
};

