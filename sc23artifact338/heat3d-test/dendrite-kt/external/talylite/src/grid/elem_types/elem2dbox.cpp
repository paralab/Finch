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
#include <talyfem/grid/elem_types/elem2dbox.h>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>


namespace TALYFEMLIB {

const int* ELEM2dBox::GetSurfaceCheckArray() const {
  static const int B4n2DCheckArray[] = {
  // {surface_id, node_id1, node_id2}
     -1, 2, 0,
     +1, 1, 3,
     -2, 0, 1,
     +2, 3, 2
  };

  static const int B9n2DCheckArray[] = {
  // {surface_id, node_id1, node_id2}
     -1, 0, 3, 6,
     +1, 2, 5, 8,
     -2, 0, 1, 2,
     +2, 6, 7, 8
  };

  /* layout of local nodes in element:
   *
   *  3-- 9-- 8-- 2
   *  |   |   |   |
   * 10--15--14-- 7
   *  |   |   |   |
   * 11--12--13-- 6
   *  |   |   |   |
   *  0-- 4-- 5-- 1
   */
  static const int B16n2DCheckArray[] = {
  // {surface_id, node_id1, node_id2}
     -1, 3, 0, 10, 11,
     +1, 1, 2, 6, 7,
     -2, 0, 1, 4, 5,
     +2, 2, 3, 8, 9,
  };

  switch (n_nodes()) {
    case 4: return B4n2DCheckArray;
    case 9: return B9n2DCheckArray;
    case 16: return B16n2DCheckArray;
    default: throw NotImplementedException();
  }
}

int ELEM2dBox::GetNodesPerSurface() const {
  switch (n_nodes()) {
    case 4: return 2;
    case 9: return 3;
    case 16: return 4;
    default: throw NotImplementedException();
  }
}

/***
    * @author maksbh
    * The normals in the case with Dendro ordering is hard coded
    * based on the surf_id. The regular logic of cross product doesn't work
    * as we can't guarantee the anti-clockwise ordering for the nodes.
    *
    */
ZEROPTV ELEM2dBox::CalculateNormal(const GRID* grid, int surface_id) const {
  ZEROPTV normal;
  switch (surface_id){
    case (-1):
      normal(0) = -1;
      break;
    case (1):
      normal(0) = 1;
      break;
    case (-2):
      normal(1) = -1;
      break;
    case (2):
      normal(1) = 1;
      break;
    default:
      std::cout << "Something went wrong for Hexahedral elements. Exiting\n";
      assert(false);
      exit(-1);
  }
  return normal;

}

kBasisFunction ELEM2dBox::basis_function() const {
  switch (n_nodes()) {
    case 4: return BASIS_LINEAR;
    case 9: return BASIS_QUADRATIC;
    case 16: return BASIS_CUBIC;
    default: throw NotImplementedException();
  }
}

}  // namespace TALYFEMLIB
