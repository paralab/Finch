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

class Elem2DBoxTest : public ElementTest {
 public:
  const char* name() override { return "Box Element (2D)"; }
  ElemType elm_type() override { return kElem2dBox; }
  int basis_order() override { return 1; }
  int nsd() override { return 2; }
  int nodes_per_surface() override { return 2; }
  int surface_count() override { return 4; }

  const ZEROPTV* node_positions() override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-1.0, -1.0),
      ZEROPTV(2.0, -1.0),
      ZEROPTV(1.0, 2.0),
      ZEROPTV(-1.0, 2.0),
    };
    return pts;
  }

  ZEROPTV normal(int surface_id) override {
    if (surface_id == -1) {
      return ZEROPTV(-1, 0, 0);
    } else if (surface_id == +1) {
      return ZEROPTV(3 / sqrt(10), 1 / sqrt(10));
    } else if (surface_id == -2) {
      return ZEROPTV(0, -1, 0);
    } else if (surface_id == +2) {
      return ZEROPTV(0, +1, 0);
    }

    return ZEROPTV(0, 0, 0);
  }

  std::vector< std::pair<ZEROPTV, bool> > inner_points() override {
    return std::vector< std::pair<ZEROPTV, bool> > {{
      { ZEROPTV(0, 0, 0), true },  // center
      { ZEROPTV(2, 0, 0), false },  // above bottom right vertex
      { ZEROPTV(-1, 0, 0), true },  // left edge
      { ZEROPTV(0, 2, 0), true },  // top edge
      { ZEROPTV(1 + 2.0 / 3.0, 0, 0), true },  // right edge
      { ZEROPTV(1 + 2.0 / 3.0 + 0.0001, 0, 0), false },  // slightly past right edge

      { ZEROPTV(-1, -1, 0), true },  // vertex 0
      { ZEROPTV(2, -1, 0), true },  // vertex 1
      { ZEROPTV(1, 2, 0), true },  // vertex 2
      { ZEROPTV(-1, 2, 0), true },  // vertex 3
    }};
  }
};
