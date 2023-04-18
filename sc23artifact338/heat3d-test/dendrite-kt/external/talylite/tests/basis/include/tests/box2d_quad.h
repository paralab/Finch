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

class Box2DQuadBasisTest : public BasisTest
{
public:
  std::string name() override { return "2D Box (Quadratic)"; }
  ElemType elm_type() override { return kElem2dBox; }
  GridType grid_type() override { return kGrid2dBox; }
  kBasisFunction basis_function() override { return BASIS_QUADRATIC; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 2; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[9] = {
      ZEROPTV(-7.0/5, -4.0/5), ZEROPTV(2.0/5, -4.0/5), ZEROPTV(2.0/5, 3.0/5), ZEROPTV(-7.0/5, 3.0/5),
      ZEROPTV(-1.0/2, -4.0/5), ZEROPTV(2.0/5, -1.0/10), ZEROPTV(-1.0/2, 3.0/5), ZEROPTV(-7.0/5, -1.0/10), ZEROPTV(-1.0/2, -1.0/10)
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[9] = {
      ZEROPTV(-(sqrt(3.0)*sqrt(5.0))/5.0, -(sqrt(3.0)*sqrt(5.0))/5.0),  // 0
      ZEROPTV(0, -(sqrt(3.0)*sqrt(5.0))/5.0),  // 4
      ZEROPTV((sqrt(3.0)*sqrt(5.0))/5.0, -(sqrt(3.0)*sqrt(5.0))/5.0),  // 1
      ZEROPTV(-(sqrt(3.0)*sqrt(5.0))/5.0, 0),  // 7
      ZEROPTV(0, 0),  // 8
      ZEROPTV((sqrt(3.0)*sqrt(5.0))/5.0, 0),  // 5
      ZEROPTV(-(sqrt(3.0)*sqrt(5.0))/5.0,  (sqrt(3.0)*sqrt(5.0))/5.0),  // 3
      ZEROPTV(0, (sqrt(3.0)*sqrt(5.0))/5.0),  // 6
      ZEROPTV((sqrt(3.0)*sqrt(5.0))/5.0,  (sqrt(3.0)*sqrt(5.0))/5.0),  // 2
    };
    return pts[i];
  }

  double weight(int i) override {
    static const double wts[9] = {
      25.0/81,  // 0
      40.0/81,  // 4
      25.0/81,  // 1
      40.0/81,  // 7
      64.0/81,  // 8
      40.0/81,  // 5
      25.0/81,  // 3
      40.0/81,  // 6
      25.0/81,  // 2
    };
    return wts[i];
  }

  double N(int itg_pt, int bf) override {
    static const double n[9][9] = {
      { pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0,2.0)*(3.0/2.0E1),
        (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0,2.0)*(3.0/2.0E1),
        (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        4.0/2.5E1 },  // 0

      { 0, 0, 0, 0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.0E1),
        0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.0E1),
        0,
        2.0/5.0 },  // 4

      { (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0,2.0)*(3.0/2.0E1),
        (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0,2.0)*(3.0/2.0E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        4.0/2.5E1 },  // 1

      { 0, 0, 0, 0, 0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.0E1),
        0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.0E1),
        2.0/5.0 },  // 7

      { 0, 0, 0, 0, 0, 0, 0, 0,
        1.0 },  // 8

      { 0, 0, 0, 0, 0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.0E1),
        0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.0E1),
        2.0/5.0 },  // 5

      { (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0,2.0)*(3.0/2.0E1),
        (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0,2.0)*(3.0/2.0E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        4.0/2.5E1 },  // 3

      { 0, 0, 0, 0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.0E1),
        0,
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.0E1),
        0,
        2.0/5.0 },  // 6

      { pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0,2.0)*(3.0/2.0E1),
        (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        pow(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0,2.0)*(3.0/2.0E1),
        (sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(3.0/2.0E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.5E1),
        sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.5E1),
        4.0/2.5E1 },  // 2
    };
    return n[itg_pt][bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[9][9] = {
      { ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)+3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)-1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)+1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)-3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)-3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)+1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)-1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)+3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(4.0/2.5E1),
                sqrt(3.0)*sqrt(5.0)*(4.0/2.5E1)) },  // 0

      { ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/2.0E1), 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.0E1), 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.0E1), 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/2.0E1), 0),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(-1.0/5.0)-1.0/2.0),
        ZEROPTV(1.0/5.0, 0),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(-1.0/5.0)+1.0/2.0),
        ZEROPTV(-1.0/5.0, 0),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(2.0/5.0)) },  // 4

      { ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)-3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)-1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)+1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)+3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)+3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)+1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)-1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)-3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-4.0/2.5E1),
                sqrt(3.0)*sqrt(5.0)*(4.0/2.5E1)) },  // 1

      { ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/2.0E1)),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/2.0E1)),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.0E1)),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.0E1)),
        ZEROPTV(0, -1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-1.0/5.0)+1.0/2.0, 0),
        ZEROPTV(0, 1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-1.0/5.0)-1.0/2.0, 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(2.0/5.0), 0) },  // 7

      { ZEROPTV(0, 0),
        ZEROPTV(0, 0),
        ZEROPTV(0, 0),
        ZEROPTV(0, 0),
        ZEROPTV(0, -1.0/2.0),
        ZEROPTV(1.0/2.0, 0),
        ZEROPTV(0, 1.0/2.0),
        ZEROPTV(-1.0/2.0, 0),
        ZEROPTV(0, 0) },  // 8

      { ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/2.0E1)),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/2.0E1)),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.0E1)),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.0E1)),
        ZEROPTV(0, -1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0, 0),
        ZEROPTV(0, 1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0, 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-2.0/5.0), 0) },  // 5

      { ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)-3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)-1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)+1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)+3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(3.0/2.5E1)+3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)+1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-2.0/2.5E1)-1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)-3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(4.0/2.5E1),
                sqrt(3.0)*sqrt(5.0)*(-4.0/2.5E1)) },  // 3

      { ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/2.0E1), 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/2.0E1), 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/2.0E1), 0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/2.0E1), 0),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0),
        ZEROPTV(1.0/5.0, 0),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0),
        ZEROPTV(-1.0/5.0, 0),
        ZEROPTV(0, sqrt(3.0)*sqrt(5.0)*(-2.0/5.0)) },  // 6

      { ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/1.0E1),
                sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/1.0E1)),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)+3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)-1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)+1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)-3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)-3.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)+1.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(2.0/2.5E1)-1.0/5.0,
                sqrt(3.0)*sqrt(5.0)*(-3.0/2.5E1)+3.0/5.0),
        ZEROPTV(sqrt(3.0)*sqrt(5.0)*(-4.0/2.5E1),
                sqrt(3.0)*sqrt(5.0)*(-4.0/2.5E1)) },  // 2
    };
    return dnde[itg_pt][bf](axis);
  }

  double dXde(int itg_pt, int i, int j) override {
    static const double dxde[2][2] = {
      { 9.0/10.0, 0 },
      { 0, 7.0/10.0 },
    };

    // same for all itg_pts
    return dxde[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    static const double cof_values[2][2] = {
      {7.0/10.0, 0},
      {0, 9.0/10.0}
    };
    return cof_values[i][j];
  }

  double jacobian_det(int itg_pt) override {
    return 63.0/100;
  }

  double dN(int itg_pt, int i, int axis) override {
    static const double dn[9][9][2] = {
     { { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(2.0/1.5E1)+2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(-4.0/3.5E1)-2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-4.0/4.5E1)+2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(6.0/3.5E1)-6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(2.0/1.5E1)-2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(-4.0/3.5E1)+2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-4.0/4.5E1)-2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(6.0/3.5E1)+6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(8.0/4.5E1),
         sqrt(3.0)*sqrt(5.0)*(8.0/3.5E1) } },  // 0

      { { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/1.8E1), 0 },
        { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.8E1), 0 },
        { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.8E1), 0 },
        { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/1.8E1), 0 },
        { 0, sqrt(3.0)*sqrt(5.0)*(-2.0/7.0)-5.0/7.0 },
        { 2.0/9.0, 0 },
        { 0, sqrt(3.0)*sqrt(5.0)*(-2.0/7.0)+5.0/7.0 },
        { -2.0/9.0, 0 },
        { 0, sqrt(3.0)*sqrt(5.0)*(4.0/7.0) } },  // 4

     { { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(-2.0/1.5E1)-2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(-4.0/3.5E1)-2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(4.0/4.5E1)+2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(6.0/3.5E1)+6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-2.0/1.5E1)+2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(-4.0/3.5E1)+2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(4.0/4.5E1)-2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(6.0/3.5E1)-6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-8.0/4.5E1),
         sqrt(3.0)*sqrt(5.0)*(8.0/3.5E1) } },  // 1

      { { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/1.4E1) },
        { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/1.4E1) },
        { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.4E1) },
        { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.4E1) },
        { 0, -2.0/7.0 },
        { sqrt(3.0)*sqrt(5.0)*(-2.0/9.0)+5.0/9.0, 0 },
        { 0, 2.0/7.0 },
        { sqrt(3.0)*sqrt(5.0)*(-2.0/9.0)-5.0/9.0, 0 },
        { sqrt(3.0)*sqrt(5.0)*(4.0/9.0), 0 } },  // 7

      { { 0, 0 },
        { 0, 0 },
        { 0, 0 },
        { 0, 0 },
        { 0, -5.0/7.0 },
        { 5.0/9.0, 0 },
        { 0, 5.0/7.0 },
        { -5.0/9.0, 0 },
        { 0, 0 } },  // 8

      { { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/1.4E1) },
        { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/1.4E1) },
        { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.4E1) },
        { 0, sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.4E1) },
        { 0, -2.0/7.0 },
        { sqrt(3.0)*sqrt(5.0)*(2.0/9.0)+5.0/9.0, 0 },
        { 0, 2.0/7.0 },
        { sqrt(3.0)*sqrt(5.0)*(2.0/9.0)-5.0/9.0, 0 },
        { sqrt(3.0)*sqrt(5.0)*(-4.0/9.0), 0 } },  // 5

     { { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(-1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(2.0/1.5E1)-2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(4.0/3.5E1)-2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-4.0/4.5E1)+2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(-6.0/3.5E1)+6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(2.0/1.5E1)+2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(4.0/3.5E1)+2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-4.0/4.5E1)-2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(-6.0/3.5E1)-6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(8.0/4.5E1),
         sqrt(3.0)*sqrt(5.0)*(-8.0/3.5E1) } },  // 3

      { { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(-1.0/1.8E1), 0 },
        { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(1.0/1.8E1), 0 },
        { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(1.0/1.8E1), 0 },
        { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(-1.0/1.8E1), 0 },
        { 0, sqrt(3.0)*sqrt(5.0)*(2.0/7.0)-5.0/7.0 },
        { 2.0/9.0, 0 },
        { 0, sqrt(3.0)*sqrt(5.0)*(2.0/7.0)+5.0/7.0 },
        { -2.0/9.0, 0 },
        { 0, sqrt(3.0)*sqrt(5.0)*(-4.0/7.0) } },  // 6

     { { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0/2.0)*(1.0/9.0),
         sqrt(3.0)*sqrt(5.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)-1.0)*(sqrt(3.0)*sqrt(5.0)*(1.0/5.0)+1.0/2.0)*(1.0/7.0) },
       { sqrt(3.0)*sqrt(5.0)*(-2.0/1.5E1)+2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(4.0/3.5E1)-2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(4.0/4.5E1)+2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(-6.0/3.5E1)-6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-2.0/1.5E1)-2.0/3.0,
         sqrt(3.0)*sqrt(5.0)*(4.0/3.5E1)+2.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(4.0/4.5E1)-2.0/9.0,
         sqrt(3.0)*sqrt(5.0)*(-6.0/3.5E1)+6.0/7.0 },
       { sqrt(3.0)*sqrt(5.0)*(-8.0/4.5E1),
         sqrt(3.0)*sqrt(5.0)*(-8.0/3.5E1) } },  // 2
    };
    return dn[itg_pt][i][axis];
  }


  double d2Nde(int itg_pt, int i, int axis) override {
    const ZEROPTV gp = gauss_point(itg_pt);
    double ksi = gp.x();
    double eta = gp.y();

    double t2 = eta-1.0;
    double t3 = eta*t2*(1.0/2.0);
    double t4 = eta*2.0;
    double t5 = t4-1.0;
    double t6 = ksi*2.0;
    double t7 = t6+1.0;
    double t8 = ksi+1.0;
    double t9 = ksi*t8*(1.0/2.0);
    double t10 = eta+1.0;
    double t11 = eta*t10*(1.0/2.0);
    double t12 = t4+1.0;
    double t13 = t6-1.0;
    double t14 = ksi-1.0;
    double t15 = ksi*t14*(1.0/2.0);
    double t16 = ksi*ksi;
    double t17 = -t16+1.0;
    double t18 = eta*eta;
    double t19 = -t18+1.0;

    double A0[9][3] = {};
    A0[0][0] = t3;
    A0[0][1] = t5*t13*(1.0/4.0);
    A0[0][2] = t15;
    A0[1][0] = t3;
    A0[1][1] = t5*t7*(1.0/4.0);
    A0[1][2] = t9;
    A0[2][0] = t11;
    A0[2][1] = t7*t12*(1.0/4.0);
    A0[2][2] = t9;
    A0[3][0] = t11;
    A0[3][1] = t12*t13*(1.0/4.0);
    A0[3][2] = t15;
    A0[4][0] = -eta*t2;
    A0[4][1] = -ksi*t5;
    A0[4][2] = t17;
    A0[5][0] = t19;
    A0[5][1] = -eta*t7;
    A0[5][2] = -ksi*t8;
    A0[6][0] = -eta*t10;
    A0[6][1] = -ksi*t12;
    A0[6][2] = t17;
    A0[7][0] = t19;
    A0[7][1] = -eta*t13;
    A0[7][2] = -ksi*t14;
    A0[8][0] = t18*2.0-2.0;
    A0[8][1] = eta*ksi*4.0;
    A0[8][2] = t16*2.0-2.0;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 3);
    return A0[i][axis];
  }

  double d2N(int itg_pt, int i, int axis) override {
    const ZEROPTV gp = gauss_point(itg_pt);
    double ksi = gp.x();
    double eta = gp.y();

  double t2 = eta-1.0;
  double t3 = eta*t2*(5.0E1/8.1E1);
  double t4 = eta*2.0;
  double t5 = t4-1.0;
  double t6 = ksi*2.0;
  double t7 = t6+1.0;
  double t8 = ksi+1.0;
  double t9 = ksi*t8*(5.0E1/4.9E1);
  double t10 = eta+1.0;
  double t11 = eta*t10*(5.0E1/8.1E1);
  double t12 = t4+1.0;
  double t13 = t6-1.0;
  double t14 = ksi-1.0;
  double t15 = ksi*t14*(5.0E1/4.9E1);
  double t16 = ksi*ksi;
  double t17 = t16*(-1.0E2/4.9E1)+1.0E2/4.9E1;
  double t18 = eta*eta;
  double t19 = t18*(-1.0E2/8.1E1)+1.0E2/8.1E1;

    double A0[9][3] = {};
  A0[0][0] = t3;
  A0[0][1] = t5*t13*(2.5E1/6.3E1);
  A0[0][2] = t15;
  A0[1][0] = t3;
  A0[1][1] = t5*t7*(2.5E1/6.3E1);
  A0[1][2] = t9;
  A0[2][0] = t11;
  A0[2][1] = t7*t12*(2.5E1/6.3E1);
  A0[2][2] = t9;
  A0[3][0] = t11;
  A0[3][1] = t12*t13*(2.5E1/6.3E1);
  A0[3][2] = t15;
  A0[4][0] = eta*t2*(-1.0E2/8.1E1);
  A0[4][1] = ksi*t5*(-1.0E2/6.3E1);
  A0[4][2] = t17;
  A0[5][0] = t19;
  A0[5][1] = eta*t7*(-1.0E2/6.3E1);
  A0[5][2] = ksi*t8*(-1.0E2/4.9E1);
  A0[6][0] = eta*t10*(-1.0E2/8.1E1);
  A0[6][1] = ksi*t12*(-1.0E2/6.3E1);
  A0[6][2] = t17;
  A0[7][0] = t19;
  A0[7][1] = eta*t13*(-1.0E2/6.3E1);
  A0[7][2] = ksi*t14*(-1.0E2/4.9E1);
  A0[8][0] = t18*(2.0E2/8.1E1)-2.0E2/8.1E1;
  A0[8][1] = eta*ksi*(4.0E2/6.3E1);
  A0[8][2] = t16*(2.0E2/4.9E1)-2.0E2/4.9E1;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 3);
    return A0[i][axis];
  }
};

