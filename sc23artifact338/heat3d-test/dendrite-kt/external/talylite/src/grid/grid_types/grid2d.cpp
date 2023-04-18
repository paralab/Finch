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
#include <talyfem/grid/grid_types/grid2d.h>

#include <talyfem/grid/elem.h>
#include <talyfem/grid/elem-types.h>  // make_elem_of_type
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

GRID2D::GRID2D(int basis_function_order)
    : GRID(basis_function_order) {
  set_grid_type(kGrid2dTriangle);
  set_nsd(2);
}

// overridden to double the number of elements because we use triangles
void GRID2D::SetNodeElemCounts(const int* n_elems, bool set_totals) {
  GRID::SetNodeElemCounts(n_elems, set_totals);
  if (set_totals) { set_n_elements(2 * n_elements()); }  // fix element count
}

// ******************************
//      Private functions
// ******************************

void GRID2D::CreateElementsBasis1() {
  int n_node_x = n_nodes_per_direction(0);

  ElemType elmType = kElem2dTriangle;

  const int kNodesPerElement = 3;
  LocalNodeID nodeIDArray1[kNodesPerElement];
  LocalNodeID nodeIDArray2[kNodesPerElement];

  for (int i = 0; i < n_elems_per_direction(0); i++) {
    for (int j = 0; j < n_elems_per_direction(1); j++) {
      /*
       * creates the two triangle elements with the local node ids as shown:
       *
       *  node 1  node 2
       *   2       2---1
       *   | \      \  |
       *   |  \      \ |
       *   0---1       0
       *
       */
      int elmID = 2 * (j * n_elems_per_direction(0) + i);

      nodeIDArray1[0] = n_node_x * j + i;
      nodeIDArray1[1] = n_node_x * j + i + 1;
      nodeIDArray1[2] = n_node_x * j + n_node_x + i;

      nodeIDArray2[0] = nodeIDArray1[1];
      nodeIDArray2[1] = n_node_x * j + n_node_x + i + 1;
      nodeIDArray2[2] = nodeIDArray1[2];

      ELEM* pElm1 = make_elem_of_type(elmType);
      elm_array_[elmID] = pElm1;
      pElm1->redim(kNodesPerElement, nodeIDArray1);
      pElm1->set_elm_id(elmID);
      pElm1->Validate(this);

      ELEM* pElm2 = make_elem_of_type(elmType);
      elm_array_[elmID + 1] = pElm2;
      pElm2->redim(kNodesPerElement, nodeIDArray2);
      pElm2->set_elm_id(elmID + 1);
      pElm2->Validate(this);
    }
  }
}

void GRID2D::CreateNodes(const double* dimensions) {
  double len_x = dimensions[0];
  double len_y = dimensions[1];

  double dx = len_x / n_elems_per_direction(0);  // x length per element
  double dy = len_y / n_elems_per_direction(1);  // y length per element
  for (int i = 0; i < n_nodes_per_direction(0); i++) {
    double x_coord = i * dx;
    for (int j = 0; j < n_nodes_per_direction(1); j++) {
      double y_coord = j * dy;
      int nodeID = n_nodes_per_direction(0) * j + i;

      NODE* pNode = new NODE();
      node_array_[nodeID] = pNode;
      pNode->setCoor(x_coord, y_coord);

      if (i == 0) { pNode->addIndicatorNum(1); }
      if (i == n_nodes_per_direction(0) - 1) { pNode->addIndicatorNum(2); }
      if (j == 0) { pNode->addIndicatorNum(3); }
      if (j == n_nodes_per_direction(1) - 1) { pNode->addIndicatorNum(4); }
    }
  }
}

}  // namespace TALYFEMLIB
