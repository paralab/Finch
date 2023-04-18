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
#include <talyfem/grid/elem_types/elem1d.h>

#include <algorithm>  // for std::swap

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>

namespace TALYFEMLIB {

const int* ELEM1d::GetSurfaceCheckArray() const {
  static const int B2n1DCheckArray[] = {
  // {surface_id, node_id}
     -1, 0,
     +1, 1
  };
  return B2n1DCheckArray;
}

int ELEM1d::GetNodesPerSurface() const {
  return 1;
}

double ELEM1d::GetMeasure(const GRID* grid) const {
  return grid->GetNode(ElemToLocalNodeID(1))->location().x()
       - grid->GetNode(ElemToLocalNodeID(0))->location().x();
}

ZEROPTV ELEM1d::CalculateNormal(const GRID* grid, int surf_id) const {
  ZEROPTV normal = grid->GetNode(ElemToLocalNodeID(1))->location()
                 - grid->GetNode(ElemToLocalNodeID(0))->location();
  normal.Normalize();

  // flip normal depending on surface
  normal = normal * (surf_id == -1 ? -1.0 : 1.0);
  return normal;
}

kBasisFunction ELEM1d::basis_function() const {
  switch (n_nodes()) {
    case 2: return BASIS_LINEAR;
    case 3: return BASIS_QUADRATIC;
    case 4: return BASIS_CUBIC;
    default: throw NotImplementedException();
  }
}

}  // namespace TALYFEMLIB
