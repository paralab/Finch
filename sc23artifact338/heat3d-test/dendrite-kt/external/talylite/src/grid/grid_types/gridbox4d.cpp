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

#include <talyfem/grid/grid_types/gridbox4d.h>

#include <talyfem/domain_decomposition/mesh_partition.h>  // for pmesh class
#include <talyfem/grid/elem.h>
#include <talyfem/grid/elem-types.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/math/math.h>


namespace TALYFEMLIB {

GRIDBox4D::GRIDBox4D(int basis_function_order)
    : GRID(basis_function_order) {
  set_nsd(4);
}

void GRIDBox4D::redimDD(const double* dimensions, const int* n_elems) {
  double len_x = dimensions[0];
  double len_y = dimensions[1];
  double len_z = dimensions[2];
  double len_t = dimensions[3];

  SetNodeElemCounts(n_elems, false);

  CMeshPartition pmesh;
  pmesh.CreateBox4D(len_x, len_y, len_z, len_t,
                    n_elems_per_direction(0), n_elems_per_direction(1),
                    n_elems_per_direction(2), n_elems_per_direction(3), basis_order());
  pmesh.TransferToGrid(this);
  // ~ pmesh.GetISNodeCopies(&isCmpGlbNodes, &isShrNodes,
  // ~                       &bISGlbAndShrNodesCreated);
  pmesh.PartitionFree();

  // called again to reset n_elems_per_direction and n_nodes_per_direction
  // arrays which were wiped by set_nsd in pmesh.TransferToGrid(this)
  // in this case, we don't reset the node and element totals
  SetNodeElemCounts(n_elems, false);
  int n_node_x = n_nodes_per_direction(0);
  int n_node_y = n_nodes_per_direction(1);
  int n_node_z = n_nodes_per_direction(2);
  int n_node_t = n_nodes_per_direction(3);

  // set boundary indicators
  for (LocalNodeID nodeID = 0; nodeID < n_nodes(); nodeID++) {
    PhysicalNodeID glbID = physical_map(nodeID);
    PetscInt i = glbID % n_node_x;
    PetscInt j = (glbID / n_node_x) % n_node_y;
    PetscInt k = (glbID / (n_node_x * n_node_y)) % n_node_z;
    PetscInt l = (glbID / (n_node_x * n_node_y * n_node_z));

    NodeIndicator indicators = 0;
    if (i == 0) { indicators |= INDICATOR_NUM(1); }
    if (i == n_node_x - 1) { indicators |= INDICATOR_NUM(2); }
    if (j == 0) { indicators |= INDICATOR_NUM(3); }
    if (j == n_node_y - 1) { indicators |= INDICATOR_NUM(4); }
    if (k == 0) { indicators |= INDICATOR_NUM(5); }
    if (k == n_node_z - 1) { indicators |= INDICATOR_NUM(6); }
    if (l == 0) { indicators |= INDICATOR_NUM(7); }
    if (l == n_node_t - 1) { indicators |= INDICATOR_NUM(8); }

    this->GetNode(nodeID)->setIndicators(indicators);
  }

  // TODO(4D): Surface integration works, but it untested now.
  // Surface integration is enabled for testing only.
  SetCaredSurfaceIndicator();
  GenElmSurfaceIndicator();
}
}  // namespace TALYFEMLIB
