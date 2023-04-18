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
#include <talyfem/grid/grid_types/grid1d.h>

#include <talyfem/grid/elem_types/elem1d.h>  // for ELEM1d class
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

GRID1D::GRID1D(int basis_function_order)
    : GRID(basis_function_order) {
  set_nsd(1);
  set_grid_type(kGrid1d);
}

// ******************************
//      Private functions
// ******************************

void GRID1D::CreateNodes(const double* dimensions) {
  double len_x = dimensions[0];
  double dx = len_x / n_elems_per_direction(0);

  // Create the nodes
  for (LocalNodeID i_node = 0; i_node < n_nodes(); i_node++) {
    double x_coord = dx * i_node;
    NODE* pNewNode = new NODE();
    pNewNode->setCoor(x_coord, 0.0, 0.0);
    node_array_[i_node] = pNewNode;
  }
  // add boundary indicators to left and right end of mesh
  node_array_[0]->addIndicatorNum(1);  // left boundary
  node_array_[n_nodes() - 1]->addIndicatorNum(2);  // right boundary
}

void GRID1D::CreateElementsBasisGeneral(int basis) {
  assert(basis > 0);  // doesn't make sense to have a zero or negative basis
  int nodes_per_element = basis + 1;
  LocalNodeID *nodeIDArray = new LocalNodeID[nodes_per_element];
  for (int elm_id = 0; elm_id < n_elements(); elm_id++) {
    ELEM* pNewElm = new ELEM1d();

    // for basis one, elements will be [0, 1], [2, 3], [4, 5], ...
    // for basis two, elements will be [0, 1, 2], [3, 4, 5], [6, 7, 8], ...
    // etc...

    // connectivity for linear is [0, 1], for quadratic is [0, 2, 1],
    // for cubic is [0, 2, 3, 1], etc.

    nodeIDArray[0] = basis * elm_id;
    nodeIDArray[1] = nodeIDArray[0] + nodes_per_element - 1;
    for (ElemNodeID i_node = 2; i_node < nodes_per_element; i_node++) {
      nodeIDArray[i_node] = nodeIDArray[0] + i_node - 1;
    }
    pNewElm->redim(nodes_per_element, nodeIDArray);

    pNewElm->GenSurfaceIndicator(this, cared_surface_indicator_);
    pNewElm->set_elm_id(elm_id);
    elm_array_[elm_id] = pNewElm;
  }
  delete [] nodeIDArray;
}

}  // namespace TALYFEMLIB
