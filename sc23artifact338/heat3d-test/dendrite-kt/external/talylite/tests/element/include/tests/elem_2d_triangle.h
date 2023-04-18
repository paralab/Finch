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

class Elem2DTriangleTest : public ElementTest {
 public:
  const char* name() override { return "Triangle Element (2D)"; }
  ElemType elm_type() override { return kElem2dTriangle; }
  int basis_order() override { return 1; }
  int nsd() override { return 2; }
  int surface_count() override { return 3; }
  int nodes_per_surface() override { return 2; }

  const ZEROPTV* node_positions() override {
    static const ZEROPTV pts[3] = {
      ZEROPTV(-2, -1, 0),
      ZEROPTV(1, 0, 0),
      ZEROPTV(-2, 2, 0),
    };
    return pts;
  }

  ZEROPTV normal(int surf_id) override {
    switch (surf_id) {
      case -1: return ZEROPTV(-1, 0, 0);
      case -2: return ZEROPTV(1/sqrt(10), -3/sqrt(10));
      case 1: return ZEROPTV(2/sqrt(13), 3/sqrt(13));
      default: return ZEROPTV();
    }
  }

  std::vector< std::pair<ZEROPTV, bool> > inner_points() override {
    return std::vector< std::pair<ZEROPTV, bool> > {{
      { ZEROPTV(0, 0, 0), true },  // center
      { ZEROPTV(-2, 0, 0), true },  // left edge
      { ZEROPTV(-2.001, 0, 0), false },  // a little past left edge
      { ZEROPTV(-1.99, -0.99, 0), true },  // just inside bottom left corner
      { ZEROPTV(-1, -1 + 1.0 / 3.0, 0), true },  // bottom edge
      { ZEROPTV(-1, 2 - 2.0 / 3.0, 0), true },  // top edge
      { ZEROPTV(-1, 2 - 2.0 / 3.0 + 0.001, 0), false },  // just past top edge
      { ZEROPTV(1.001, 0, 0), false },  // to the right of vertex 1

      { ZEROPTV(-2, -1, 0), true },  // vertex 0
      { ZEROPTV(1, 0, 0), true },  // vertex 1
      { ZEROPTV(-2, 2, 0), true },  // vertex 2
    }};
  }
};
