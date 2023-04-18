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
#include <talyfem/grid/elem_types/elem2dtriangle.h>

#include <algorithm>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>

#include <talyfem/utils/utils.h>

namespace TALYFEMLIB {

/**
 * Node ordering for triangle elements (taken from Gmsh docs)
 *
 * Triangle:               Triangle6:          Triangle9/10:
 *
 *  v
 *  ^
 *  |
 *  2                       2                    2
 *  |`\                     |`\                  | \
 *  |  `\                   |  `\                7   6
 *  |    `\                 5    `4              |     \
 *  |      `\               |      `\            8  (9)  5
 *  |        `\             |        `\          |         \
 *  0----------1 --> u      0-----3----1         0---3---4---1
 * @return
 */

const int* ELEM2dTriangle::GetSurfaceCheckArray() const {
  static const int T3n2DCheckArray[] = {
  // {surface_id, node_id1, node_id2}
    -1, 2, 0,
    +1, 1, 2,
    -2, 0, 1
  };

  static const int T6n2DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3}
      -1, 2, 0, 5,
      +1, 1, 2, 4,
      -2, 0, 1, 3
  };

  // 9 and 10 node triangles have the same surfaces
  // (the extra node is on the interior)
  static const int T9n2DCheckArray[] = {
      // {surface_id, node_id1, node_id2}
      -1, 2, 0, 7, 8,
      +1, 1, 2, 5, 6,
      -2, 0, 1, 3, 4
  };

  switch (n_nodes()) {
    case 3: return T3n2DCheckArray;
    case 6: return T6n2DCheckArray;
    case 9: return T9n2DCheckArray;
    case 10: return T9n2DCheckArray;  // intentionally the same
    default: throw NotImplementedException();
  }
}

int ELEM2dTriangle::GetNodesPerSurface() const {
  switch (n_nodes()) {
    case 3: return 2;
    case 6: return 3;
    case 9: return 4;
    case 10: return 4;
    default: throw NotImplementedException();
  }
}

// This implementation assumes counter-clockwise connectivity
// and surface array node order. If this not the case, the normal
// will point inwards (bad!).
ZEROPTV ELEM2dTriangle::CalculateNormal(const GRID* grid, int surf_id) const {
  const int* surf_nodes = GetNodesInSurface(surf_id);
  const int nodes[] = { surf_nodes[0], surf_nodes[1], (surf_nodes[1] + 1) % 3 };

  const ZEROPTV& p1 = grid->GetNode(ElemToLocalNodeID(nodes[0]))->location();
  const ZEROPTV& p2 = grid->GetNode(ElemToLocalNodeID(nodes[1]))->location();
  const ZEROPTV& p3 = grid->GetNode(ElemToLocalNodeID(nodes[2]))->location();
  ZEROPTV elem_normal;
  elem_normal.crossProduct(p2 - p1, p3 - p1);

  ZEROPTV edge_normal;
  edge_normal.crossProduct(p2 - p1, elem_normal);
  edge_normal.Normalize();

  return edge_normal;
}

void ELEM2dTriangle::Validate(const GRID* grid) {
  /*const ZEROPTV edge1 = grid->GetNode(ElemToLocalNodeID(1))->location()
                      - grid->GetNode(ElemToLocalNodeID(0))->location();
  const ZEROPTV edge2 = grid->GetNode(ElemToLocalNodeID(2))->location()
                      - grid->GetNode(ElemToLocalNodeID(0))->location();
  ZEROPTV cross;
  cross.crossProduct(edge1, edge2);
  if (cross.z() <= 0) {
    throw TALYException() << "Element " << elm_id() << " has clockwise "
       "connectivity. TalyFEM only supports counter-clockwise windings.";
  }*/
}

kBasisFunction ELEM2dTriangle::basis_function() const {
  switch (n_nodes()) {
    case 3: return BASIS_LINEAR;
    case 6: return BASIS_QUADRATIC;
    case 9: return BASIS_CUBIC;
    case 10: return BASIS_CUBIC;
    default: throw NotImplementedException();
  }
}

}  // namespace TALYFEMLIB
