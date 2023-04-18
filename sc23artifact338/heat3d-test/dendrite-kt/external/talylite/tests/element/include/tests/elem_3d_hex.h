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

class Elem3DHexTest : public ElementTest {
 public:
  const char* name() override { return "Hexahedral Element (3D)"; }
  ElemType elm_type() override { return kElem3dHexahedral; }
  int basis_order() override { return 1; }
  int nsd() override { return 3; }

  int surface_count() override { return 6; }
  int nodes_per_surface() override { return 4; }

  const ZEROPTV* node_positions() override {
    static const ZEROPTV pts[] = {
      ZEROPTV(-1.5, -1, -2),
      ZEROPTV(+1, -1, -1),
      ZEROPTV(+1, +1, -1),
      ZEROPTV(-1.5, +1, -2),

      ZEROPTV(-1, -1, +1),
      ZEROPTV(+1, -1, +1),
      ZEROPTV(+1, +1, +1),
      ZEROPTV(-1, +1, +1),
    };
    return pts;
  }

  ZEROPTV normal(int surface_id) override {
    switch (surface_id) {
      case -1: return ZEROPTV(-6 / sqrt(37), 0, 1 / sqrt(37));
      case +1: return ZEROPTV(1, 0, 0);
      case -2: return ZEROPTV(0, -1, 0);
      case +2: return ZEROPTV(0, +1, 0);
      case -3: return ZEROPTV(2 / sqrt(29), 0, -5 / sqrt(29));
      case +3: return ZEROPTV(0, 0, +1);
      default: return ZEROPTV();
    }
  }

  std::vector< std::pair<ZEROPTV, bool> > inner_points() override {
    return std::vector< std::pair<ZEROPTV, bool> > {{
      { ZEROPTV(0, 0, 0), true },  // center

      // random points
      { ZEROPTV(-1.5, -1, -3), false },
      { ZEROPTV(-1.6, -1, -2.1), false },
      { ZEROPTV(+1, -1.1, -1), false },
      { ZEROPTV(+0.9, -0.9, -1), true },
      { ZEROPTV(-0.9, 0.8, 0.9), true },
      { ZEROPTV(-0.9, 0.8, 1.2), false },

      // vertices
      { ZEROPTV(-1.5, -1, -2), true },
      { ZEROPTV(+1, -1, -1), true },
      { ZEROPTV(+1, +1, -1), true },
      { ZEROPTV(-1.5, +1, -2), true },
      { ZEROPTV(-1, -1, +1), true },
      { ZEROPTV(+1, -1, +1), true },
      { ZEROPTV(+1, +1, +1), true },
      { ZEROPTV(-1, +1, +1), true },
    }};
  }
};
