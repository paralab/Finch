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
#include <talyfem/grid/elem_types/elem3dhexahedral.h>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>


namespace TALYFEMLIB {

const int* ELEM3dHexahedral::GetSurfaceCheckArray() const {
  static const int B8n3DCheckArray[] = {
  // {surface_id, node_id1, node_id2, node_id3, node_id4 }
       -1, 0, 2, 4, 6,
       +1, 1, 3, 5, 7,
       -2, 0, 1, 4, 5,
       +2, 2, 3, 6, 7,
       -3, 0, 1, 2, 3,
       +3, 4, 5, 6, 7
  };

  static const int B27n3DCheckArray[] = {
  // {surface_id, node_id1, node_id2, node_id3, node_id4 }
    -1, 0, 3, 6, 9, 12, 15, 18, 21, 24,
    +1, 2, 5, 8, 11, 14, 17, 20, 23, 26,
    -2, 0, 1, 2, 9, 10, 11, 18, 19, 20,
    +2, 6, 7, 8, 15, 16, 17, 24, 25, 26,
    -3, 0, 1, 2, 3, 4, 5, 6, 7, 8,
    +3, 18, 19, 20, 21, 22, 23, 24, 25, 26
  };

  static const int B64n3DCheckArray[] = {
  // {surface_id, node_id1, node_id2, node_id3, node_id4 }
    -1, 0, 4, 7, 3, 8, 12, 35, 34, 15, 11, 22, 23, 47, 59, 58, 46,
    +1, 5, 1, 2, 6, 13, 9, 18, 19, 10, 14, 31, 30, 54, 42, 43, 55,
    -2, 0, 1, 5, 4, 16, 17, 9, 13, 29, 28, 12, 8, 40, 41, 53, 52,
    +2, 7, 6, 2, 3, 33, 32, 14, 10, 20, 21, 11, 15, 57, 56, 44, 45,
    -3, 1, 0, 3, 2, 17, 16, 23, 22, 21, 20, 19, 18, 25, 24, 27, 26,
    +3, 4, 5, 6, 7, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
  };

  switch (n_nodes()) {
    case 8: return B8n3DCheckArray;
    case 27: return B27n3DCheckArray;
    case 64: return B64n3DCheckArray;
    default: throw NotImplementedException();
  }
}

int ELEM3dHexahedral::GetNodesPerSurface() const {
  switch (n_nodes()) {
    case 8: return 4;
    case 27: return 9;
    case 64: return 16;
    default: throw NotImplementedException();
  }
}

ZEROPTV ELEM3dHexahedral::CalculateNormal(const GRID* grid, int surf_id) const {
    /***
     * @author maksbh
     * The normals in the case with Dendro ordering is hard coded
     * based on the surf_id. The regular logic of cross product doesn't work
     * as we can't guarantee the anti-clockwise ordering for the nodes.
     *
     */
  ZEROPTV normal;
  switch (surf_id){
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
      case (-3):
          normal(2) = -1;
          break;
      case 3:
          normal(2) = 1;
          break;
      default:
          std::cout << "Something went wrong for Hexahedral elements. Exiting\n";
          assert(false);
          exit(-1);
  }
  return normal;
}

kBasisFunction ELEM3dHexahedral::basis_function() const {
  switch (n_nodes()) {
    case 8: return BASIS_LINEAR;
    case 27: return BASIS_QUADRATIC;
    case 64: return BASIS_CUBIC;
    default: throw NotImplementedException();
  }
}

}  // namespace TALYFEMLIB
