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

#include <element_test.h>

class Elem3DTetTest : public ElementTest {
 public:
  const char* name() override { return "Tetrahedron Element (3D)"; }
  ElemType elm_type() override { return kElem3dTetrahedral; }
  int basis_order() override { return 1; }
  int nsd() override { return 3; }

  int nodes_per_surface() override { return 3; }
  int surface_count() override { return 4; }

  const ZEROPTV* node_positions() override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-3.88, 2.06, -1.1),   // 4
      ZEROPTV(-1.47, -0.83, -0.7),  // 1
      ZEROPTV(0.51, -1.4, 1.3),     // 2
      ZEROPTV(-0.95, 0.27, 4.8),    // 3
    };
    return pts;
  }

  ZEROPTV normal(int surf_id) override {
    switch (surf_id) {
      case 1: return ZEROPTV(-5.335 / sqrt(131.60738036), -9.85 / sqrt(131.60738036), 2.4744 / sqrt(131.60738036));
      case 2: return ZEROPTV(16.118 / sqrt(621.02611709), 18.869 / sqrt(621.02611709), -2.2797 / sqrt(621.02611709));
      case 3: return ZEROPTV(-16.335 / sqrt(454.31048844), -13.047 / sqrt(454.31048844), 4.1538 / sqrt(454.31048844));
      case 4: return ZEROPTV(5.552 / sqrt(65.95894025), 4.028 / sqrt(65.95894025), -4.3485 / sqrt(65.95894025));
      default: return ZEROPTV();
    }
  }

  double measure() override {
    return 2.7664850219327195013;  // according to paraview
  }

  std::vector< std::pair<ZEROPTV, bool> > inner_points() override {
    return std::vector< std::pair<ZEROPTV, bool> > {{
      { ZEROPTV(-1.4475, 0.025, 1.075), true },  // midpoint

      { ZEROPTV(-1.47, -0.83, -0.7), true },  // vertex 0
      { ZEROPTV(0.51, -1.4, 1.3), true },  // vertex 1
      { ZEROPTV(-0.95, 0.27, 4.8), true },  // vertex 2
      { ZEROPTV(-3.88, 2.06, -1.1), true },  // vertex 3
    }};
  }
};
