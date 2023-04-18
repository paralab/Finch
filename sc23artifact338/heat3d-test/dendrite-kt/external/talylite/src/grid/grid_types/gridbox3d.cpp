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
#include <talyfem/grid/grid_types/gridbox3d.h>

#include <talyfem/common/indexer3d.h>  // for Indexer3D class
#include <talyfem/domain_decomposition/mesh_partition.h>  // for pmesh class
#include <talyfem/grid/elem.h>
#include <talyfem/grid/elem-types.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/math/math.h>  // for qbf3DIDarr and cbf3DIDarr arrays


namespace TALYFEMLIB {

GRIDBox3D::GRIDBox3D(int basis_function_order)
    : GRID(basis_function_order) {
  set_nsd(3);
}

void GRIDBox3D::redimDD(const double* dimensions, const int* n_elems) {
  double len_x = dimensions[0];
  double len_y = dimensions[1];
  double len_z = dimensions[2];

  SetNodeElemCounts(n_elems, false);

  CMeshPartition pmesh;
  pmesh.CreateBox3D(len_x, len_y, len_z,
                    n_elems_per_direction(0), n_elems_per_direction(1),
                    n_elems_per_direction(2), basis_order());
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

  // set boundary indicators
  for (LocalNodeID nodeID = 0; nodeID < n_nodes(); nodeID++) {
    PhysicalNodeID glbID = physical_map(nodeID);
    PetscInt i = glbID / (n_node_x * n_node_y);
    PetscInt j = (glbID % (n_node_x * n_node_y)) / n_node_x;
    PetscInt k = (glbID % (n_node_x * n_node_y)) % n_node_x;

    NodeIndicator indicators = 0;
    if (k == 0) { indicators |= INDICATOR_NUM(1); }
    if (k == n_node_x - 1) { indicators |= INDICATOR_NUM(2); }
    if (j == 0) { indicators |= INDICATOR_NUM(3); }
    if (j == n_node_y - 1) { indicators |= INDICATOR_NUM(4); }
    if (i == 0) { indicators |= INDICATOR_NUM(5); }
    if (i == n_node_z - 1) { indicators |= INDICATOR_NUM(6); }

    this->GetNode(nodeID)->setIndicators(indicators);
  }

  SetCaredSurfaceIndicator();
  GenElmSurfaceIndicator();
}

// ******************************
//      Private functions
// ******************************

void GRIDBox3D::CreateElementsBasis1() {
  const ElemNodeID kNodesPerElem = 8;
  LocalNodeID nodeIDArray[kNodesPerElem];

  Indexer3D indexer(n_nodes_per_direction(0), n_nodes_per_direction(1),
                    n_nodes_per_direction(2));
  Indexer3D elem_indexer(n_elems_per_direction(0), n_elems_per_direction(1),
                         n_elems_per_direction(2));
  for (int i = 0; i < n_elems_per_direction(0); i++) {
    for (int j = 0; j < n_elems_per_direction(1); j++) {
      for (int k = 0; k < n_elems_per_direction(2); k++) {
        /* layout of local nodes in element:
         *
         * 7------6
         * |\     |\
         * | 3------2
         * | |    | |
         * 4-|----5 |
         *  \|     \|
         *   0------1
         */

        nodeIDArray[0] = indexer(i, j, k);
        nodeIDArray[1] = indexer(i + 1, j, k);
        nodeIDArray[2] = indexer(i + 1, j + 1, k);
        nodeIDArray[3] = indexer(i, j + 1, k);
        nodeIDArray[4] = indexer(i, j, k + 1);
        nodeIDArray[5] = indexer(i + 1, j, k + 1);
        nodeIDArray[6] = indexer(i + 1, j + 1, k + 1);
        nodeIDArray[7] = indexer(i, j + 1, k + 1);

        const int elem_id = elem_indexer(i, j, k);
        ELEM* pElm = make_elem_of_type(kElem3dHexahedral);
        elm_array_[elem_id] = pElm;
        pElm->redim(kNodesPerElem, nodeIDArray);
        pElm->set_elm_id(elem_id);
        pElm->Validate(this);
      }
    }
  }
}

void GRIDBox3D::CreateElementsBasis2() {
  const ElemNodeID kNodesPerElem = 27;
  LocalNodeID nodeIDArray[kNodesPerElem];

  Indexer3D indexer(n_nodes_per_direction(0), n_nodes_per_direction(1),
                    n_nodes_per_direction(2));
  Indexer3D elem_indexer(n_elems_per_direction(0) / 2,
                         n_elems_per_direction(1) / 2,
                         n_elems_per_direction(2) / 2);
  // actual number of elements in each direction is half of desired number
  for (int i = 0; i < n_elems_per_direction(0); i += 2) {
    for (int j = 0; j < n_elems_per_direction(1); j += 2) {
      for (int k = 0; k < n_elems_per_direction(2); k += 2) {
        /* layout of local nodes in element:
         *
         *   7-------19--------6
         *   |\       |\       |\
         *   |11-------24-------10
         *   | |\     | |\     | |\
         *   | | 3-------14--------2
         *  20-|-|---21-|-|---18 | |
         *   |\| |    |\| |    |\| |
         *   |25-|-----26-|-----23 |
         *   | |\|    | |\|    | |\|
         *   | |15-------16-------13
         *   4-|-|---17-|-|----5 | |
         *    \| |     \| |     \| |
         *     8-|-----22-|------9 |
         *      \|       \|       \|
         *       0-------12--------1
         */

        for (int k1 = 0; k1 < 3; k1++) {
          for (int j1 = 0; j1 < 3; j1++) {
            for (int i1 = 0; i1 < 3; i1++) {
              const int id = 9 * k1 + 3 * j1 + i1;
              nodeIDArray[qbf3DIDarr[id] - 1] = indexer(i + i1, j + j1, k + k1);
            }
          }
        }

        const int elem_id = elem_indexer(i / 2, j / 2, k / 2);
        ELEM* pElm = make_elem_of_type(kElem3dHexahedral);
        elm_array_[elem_id] = pElm;
        pElm->redim(kNodesPerElem, nodeIDArray);
        pElm->set_elm_id(elem_id);
        pElm->Validate(this);
      }
    }
  }
}

void GRIDBox3D::CreateElementsBasis3() {
  const ElemNodeID kNodesPerElem = 64;
  LocalNodeID nodeIDArray[kNodesPerElem];

  Indexer3D indexer(n_nodes_per_direction(0), n_nodes_per_direction(1),
                    n_nodes_per_direction(2));
  Indexer3D elem_indexer(n_elems_per_direction(0) / 3,
                         n_elems_per_direction(1) / 3,
                         n_elems_per_direction(2) / 3);
  // actual number of elements in each direction one third of desired number
  for (int i = 0; i < n_elems_per_direction(0); i += 3) {
    for (int j = 0; j < n_elems_per_direction(1); j += 3) {
      for (int k = 0; k < n_elems_per_direction(2); k += 3) {
        // construct the element
        for (int k1 = 0; k1 < 4; k1++) {
          for (int j1 = 0; j1 < 4; j1++) {
            for (int i1 = 0; i1 < 4; i1++) {
              const int id = 16 * k1 + 4 * j1 + i1;
              nodeIDArray[cbf3DIDarr[id] - 1] = indexer(i + i1, j + j1, k + k1);
            }
          }
        }

        const int elem_id = elem_indexer(i / 3, j / 3, k / 3);
        ELEM* pElm = make_elem_of_type(kElem3dHexahedral);
        elm_array_[elem_id] = pElm;
        pElm->redim(kNodesPerElem, nodeIDArray);
        pElm->set_elm_id(elem_id);
        pElm->Validate(this);
      }
    }
  }
}

void GRIDBox3D::CreateNodes(const double* dimensions) {
  double dx = dimensions[0] / n_elems_per_direction(0);
  double dy = dimensions[1] / n_elems_per_direction(1);
  double dz = dimensions[2] / n_elems_per_direction(2);

  for (int i = 0; i < n_nodes_per_direction(0); i++) {
    double x_coord = i * dx;
    for (int j = 0; j < n_nodes_per_direction(1); j++) {
      double y_coord = j * dy;
      for (int k = 0; k < n_nodes_per_direction(2); k++) {
        double z_coord = k * dz;
        int node_id = n_nodes_per_direction(0) * n_nodes_per_direction(1) * k +
                      n_nodes_per_direction(0) * j + i;

        NODE* pNode = new NODE();
        node_array_[node_id] = pNode;
        pNode->setCoor(x_coord, y_coord, z_coord);

        if (i == 0) { pNode->addIndicatorNum(1); }
        if (i == n_nodes_per_direction(0) - 1) { pNode->addIndicatorNum(2); }
        if (j == 0) { pNode->addIndicatorNum(3); }
        if (j == n_nodes_per_direction(1) - 1) { pNode->addIndicatorNum(4); }
        if (k == 0) { pNode->addIndicatorNum(5); }
        if (k == n_nodes_per_direction(2) - 1) { pNode->addIndicatorNum(6); }
      }
    }
  }
}

/* void GRIDBox3D::SetBndrIndicators(InputData* pIdata){

 for(int nodeID=1;nodeID<=n_nodes();nodeID++)
 {
 int inds[6];
 int nInd=0;
 double x = this->getCoor (nodeID, 1);

 double epsilon_x = 0.5 * (pIdata->L[0] / (pIdata->Nelem[0]));
 if(x< (0.0 + epsilon_x))
 {
 inds[nInd]=1;
 nInd++;
 }
 if(x> (pIdata->L[0] - epsilon_x))
 {
 inds[nInd]=2;
 nInd++;
 }

 double epsilon_y = 0.5 * (pIdata->L[1] / (pIdata->Nelem[1]));
 double y = this->getCoor (nodeID, 2);
 if(y< (0.0 + epsilon_y))
 {
 inds[nInd]=3;
 nInd++;
 }
 if(y> (pIdata->L[1] - epsilon_y))
 {
 inds[nInd]=4;
 nInd++;
 }

 double epsilon_z = 0.5 * (pIdata->L[2] / (pIdata->Nelem[2]));
 double z = this->getCoor (nodeID, 3);
 if(z<(0.0 + epsilon_z))
 {
 inds[nInd]=5;
 nInd++;
 }
 if((z>pIdata->L[2] - epsilon_z))
 {
 inds[nInd]=6;
 nInd++;
 }

 this->GetNode(nodeID)->setIndicator(nInd,inds);
 }
 }
*/

}  // namespace TALYFEMLIB
