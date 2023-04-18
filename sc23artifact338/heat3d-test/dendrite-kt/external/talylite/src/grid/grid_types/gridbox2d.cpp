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
#include <talyfem/grid/grid_types/gridbox2d.h>

#include <talyfem/common/indexer2d.h>  // for Indexer2D class
#include <talyfem/domain_decomposition/mesh_partition.h>  // for pmesh class
#include <talyfem/grid/elem.h>
#include <talyfem/grid/elem-types.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/math/math.h>  // for insertValue


namespace TALYFEMLIB {

GRIDBox2D::GRIDBox2D(int basis_function_order)
    : GRID(basis_function_order),
      bottomline_(),
      topline_(),
      len_x_(0.0) {
  set_nsd(2);
  set_grid_type(kGrid2dBox);
}

int GRIDBox2D::GetBottomNodeID(int nodeXID) const {
  return nodeXID;
}

int GRIDBox2D::GetElmID(int nXID, int nYID) const {
  return nYID * n_elems_per_direction(0) + nXID;
}

int GRIDBox2D::GetElmXID(int elmID) const {
  return elmID % n_elems_per_direction(0);
}

// note: this is per-node, not per element, so +1
int GRIDBox2D::GetNodeID(int i, int j) const {
  if (i < 0) { return -1; }
  if (i >= n_nodes_per_direction(0)) { return -2; }
  if (j < 0) { return -3; }
  if (j >= n_nodes_per_direction(1)) { return -4; }

  return j * (n_nodes_per_direction(0)) + i;
}

int GRIDBox2D::GetNodeXID(int node_id) const {
  return node_id % (n_nodes_per_direction(0));
}

int GRIDBox2D::GetTopNodeID(int nodeXID) const {
  return n_elems_per_direction(1) * (n_nodes_per_direction(0)) + nodeXID;
}

void GRIDBox2D::GetXYID(int node_id, int& I, int& J) const {
  I = GetNodeXID(node_id);
  J = (node_id - I) / (n_nodes_per_direction(0));
}

double GRIDBox2D::GetY(int i, int j) const {
  // insertValue is a templated function defined in Math/MathEx.h.
  return insertValue(1.0, static_cast<double>(n_nodes_per_direction(1)),
                     bottomline_(i), topline_(i), static_cast<double>(j + 1));
}

void GRIDBox2D::redim(const double* dimensions, const int* n_elems) {
  ValidateParams(dimensions, n_elems);

  double xlength = dimensions[0];
  double ylength = dimensions[1];
  int nGridX = n_elems[0];
  int nGridY = n_elems[1];

  redimCos(0, xlength, 1, ylength, nGridX, nGridY);
}

void GRIDBox2D::redimAbove(const GRIDBox2D& grid, double maxY, int nGridY) {
  int nGridX = grid.n_elems_per_direction(0);
  double* bottomline = new double[nGridX + 1];
  double* topline = new double[nGridX + 1];

  // generate the bottomline and topline
  double DX = grid.len_x_;
  for (int i = 0; i < nGridX + 1; i++) {
    topline[i] = maxY;
    bottomline[i] = grid.topline_(i);
  }

  redimGeneral(DX, nGridX, nGridY, bottomline, topline);
  delete[] bottomline;
  delete[] topline;
}

void GRIDBox2D::redimCos(double amplitude, double wavelength,
                         int noOfWavelength, double DY, int nGridX,
                         int nGridY) {
  double* bottomline = new double[nGridX + 1];
  double* topline = new double[nGridX + 1];

  // generate the bottomline and topline
  double DX = wavelength * noOfWavelength;
  double dx = DX / nGridX;
  for (int i = 0; i < nGridX + 1; i++) {
    bottomline[i] = amplitude * cos(i * dx / (wavelength) * 2 * 3.1415926535);
    topline[i] = DY;
  }

  redimGeneral(DX, nGridX, nGridY, bottomline, topline);
  delete[] bottomline;
  delete[] topline;
}

void GRIDBox2D::redimCosMold(double amplitude, double wavelength,
                             int noOfWavelength, double thick, int nGridX,
                             int nGridY) {
  double* bottomline = new double[nGridX + 1];
  double* topline = new double[nGridX + 1];

  // generate the bottomline and topline
  double DX = wavelength * noOfWavelength;
  double dx = DX / nGridX;
  for (int i = 0; i < nGridX + 1; i++) {
    topline[i] = amplitude * cos(i * dx / (wavelength) * 2 * 3.1415926535);
    bottomline[i] = -thick;
  }

  redimGeneral(DX, nGridX, nGridY, bottomline, topline);
  delete[] bottomline;
  delete[] topline;
}

void GRIDBox2D::redimCosShell(double amplitude, double wavelength,
                              int noOfWavelength, double thick, int nGridX,
                              int nGridY) {
  double* bottomline = new double[nGridX + 1];
  double* topline = new double[nGridX + 1];

  // generate the bottomline and topline
  double DX = wavelength * noOfWavelength;
  double dx = DX / nGridX;
  for (int i = 0; i < nGridX + 1; i++) {
    bottomline[i] = amplitude * cos(i * dx / (wavelength) * 2 * 3.1415926535);
    topline[i] = bottomline[i] + thick;
  }

  redimGeneral(DX, nGridX, nGridY, bottomline, topline);
  delete[] bottomline;
  delete[] topline;
}

void GRIDBox2D::redimDD(const double* dimensions, const int* n_elems) {
  double xlength = dimensions[0];
  double ylength = dimensions[1];

  SetNodeElemCounts(n_elems, false);

  CMeshPartition pmesh;
  pmesh.CreateBox2D(xlength, ylength,
                    n_elems_per_direction(0),
                    n_elems_per_direction(1), basis_order());
  pmesh.TransferToGrid(this);
  // ~ pmesh.GetISNodeCopies(&isCmpGlbNodes, &isShrNodes,
  // ~                       &bISGlbAndShrNodesCreated);
  // ~ pmesh.PrintToParallelFile("box.dat");
  pmesh.PartitionFree();

  // called again to reset n_elems_per_direction and n_nodes_per_direction
  // arrays which were wiped by set_nsd in pmesh.TransferToGrid(this)
  // in this case, we don't reset the node and element totals
  SetNodeElemCounts(n_elems, false);
  int n_node_x = n_nodes_per_direction(0);
  int n_node_y = n_nodes_per_direction(1);

  // set boundary indicators
  for (LocalNodeID nodeID = 0; nodeID < n_nodes(); nodeID++) {
    NodeIndicator indicators = 0;

    PhysicalNodeID glbID = physical_map(nodeID);
    PetscInt i = glbID / n_node_x;
    PetscInt j = glbID % n_node_x;

    if (j == 0) { indicators |= INDICATOR_NUM(1); }
    if (j == n_node_x - 1) { indicators |= INDICATOR_NUM(2); }
    if (i == 0) { indicators |= INDICATOR_NUM(3); }
    if (i == n_node_y - 1) { indicators |= INDICATOR_NUM(4); }

    this->GetNode(nodeID)->setIndicators(indicators);
  }
  SetCaredSurfaceIndicator();
  GenElmSurfaceIndicator();
}

void GRIDBox2D::redimGeneral(double DX, int nGridX, int nGridY,
                             const double bottomline[],
                             const double topline[]) {
  int n_elems[2];
  n_elems[0] = nGridX;
  n_elems[1] = nGridY;

  SetNodeElemCounts(n_elems);

  redimArrays(n_nodes(), n_elements());

  int n_node_x = n_nodes_per_direction(0);

  // copy the bottomline and topline
  bottomline_.fill_from_array(bottomline, n_node_x);
  topline_.fill_from_array(topline, n_node_x);

  len_x_ = DX;

  CreateNodes(NULL);

  switch (basis_order()) {
    case 1:
      CreateElementsBasis1();
      break;
    case 2:
      CreateElementsBasis2();
      break;
    case 3:
      CreateElementsBasis3();
      break;
    default:
      throw NotImplementedException() <<
          "Unsupported basis order for GridBox2D redim.";
  }
  SetCaredSurfaceIndicator();
  GenElmSurfaceIndicator();
  set_grid_type(kGrid2dBox);
}

// ******************************
//      Private functions
// ******************************

void GRIDBox2D::CreateElementsBasis1() {
  const ElemNodeID kNodesPerElem = 4;
  LocalNodeID nodeIDArray[kNodesPerElem];

  int n_node_x = n_nodes_per_direction(0);
  int n_node_y = n_nodes_per_direction(1);
  int n_true_elem_x = n_elems_per_direction(0);
  int n_true_elem_y = n_elems_per_direction(1);
  Indexer2D indexer(n_node_x, n_node_y);
  for (int i = 0; i < n_true_elem_x; i++) {
    for (int j = 0; j < n_true_elem_y; j++) {
      PetscInt elmID = n_true_elem_x * j + i;

      /* layout of local nodes in element:
       *
       * 3--2
       * |  |
       * 0--1
       *
       */
      nodeIDArray[0] = indexer(i, j);
      nodeIDArray[1] = indexer(i + 1, j);
      nodeIDArray[2] = indexer(i + 1, j + 1);
      nodeIDArray[3] = indexer(i, j + 1);

      ELEM* pElm = make_elem_of_type(kElem2dBox);
      elm_array_[elmID] = pElm;
      pElm->redim(kNodesPerElem, nodeIDArray);
      pElm->set_elm_id(elmID);
      pElm->Validate(this);
    }
  }
}

void GRIDBox2D::CreateElementsBasis2() {
  const ElemNodeID kNodesPerElem = 9;
  LocalNodeID nodeIDArray[kNodesPerElem];

  int n_node_x = n_nodes_per_direction(0);
  int n_node_y = n_nodes_per_direction(1);
  Indexer2D indexer(n_node_x, n_node_y);
  // actual number of elements in each direction is half of desired number
  for (int i = 0; i < n_elems_per_direction(0); i += 2) {
    for (int j = 0; j < n_elems_per_direction(1); j += 2) {
      // parentheses are necessary because this is integer division
      PetscInt elmID = (n_elems_per_direction(0) * j) / 4 + i / 2;
      /* layout of local nodes in element:
       *
       * 3--6--2
       * |  |  |
       * 7--8--5
       * |  |  |
       * 0--4--1
       *
       */
      nodeIDArray[0] = indexer(i, j);
      nodeIDArray[1] = indexer(i + 2, j);
      nodeIDArray[2] = indexer(i + 2, j + 2);
      nodeIDArray[3] = indexer(i, j + 2);
      nodeIDArray[4] = indexer(i + 1, j);
      nodeIDArray[5] = indexer(i + 2, j + 1);
      nodeIDArray[6] = indexer(i + 1, j + 2);
      nodeIDArray[7] = indexer(i, j + 1);
      nodeIDArray[8] = indexer(i + 1, j + 1);

      ELEM* pElm = make_elem_of_type(kElem2dBox);
      elm_array_[elmID] = pElm;
      pElm->redim(kNodesPerElem, nodeIDArray);
      pElm->set_elm_id(elmID);
      pElm->Validate(this);
    }
  }
}

void GRIDBox2D::CreateElementsBasis3() {
  const ElemNodeID kNodesPerElem = 16;
  LocalNodeID nodeIDArray[kNodesPerElem];

  int n_node_x = n_nodes_per_direction(0);
  int n_node_y = n_nodes_per_direction(1);
  Indexer2D indexer(n_node_x, n_node_y);
  // actual number of elements in each direction is one third of desired number
  for (int i = 0; i < n_elems_per_direction(0); i += 3) {
    for (int j = 0; j < n_elems_per_direction(1); j += 3) {
      // parentheses are necessary because this is integer division
      PetscInt elmID = (n_elems_per_direction(0) * j) / 9 + i / 3;

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
      nodeIDArray[0] = indexer(i, j);
      nodeIDArray[1] = indexer(i + 3, j);
      nodeIDArray[2] = indexer(i + 3, j + 3);
      nodeIDArray[3] = indexer(i, j + 3);
      nodeIDArray[4] = indexer(i + 1, j);
      nodeIDArray[5] = indexer(i + 2, j);
      nodeIDArray[6] = indexer(i + 3, j + 1);
      nodeIDArray[7] = indexer(i + 3, j + 2);
      nodeIDArray[8] = indexer(i + 2, j + 3);
      nodeIDArray[9] = indexer(i + 1, j + 3);
      nodeIDArray[10] = indexer(i, j + 2);
      nodeIDArray[11] = indexer(i, j + 1);
      nodeIDArray[12] = indexer(i + 1, j + 1);
      nodeIDArray[13] = indexer(i + 2, j + 1);
      nodeIDArray[14] = indexer(i + 2, j + 2);
      nodeIDArray[15] = indexer(i + 1, j + 2);

      ELEM* pElm = make_elem_of_type(kElem2dBox);
      elm_array_[elmID] = pElm;
      pElm->redim(kNodesPerElem, nodeIDArray);
      pElm->set_elm_id(elmID);
      pElm->Validate(this);
    }
  }
}

void GRIDBox2D::CreateNodes(const double* dimensions) {
  int n_node_x = n_nodes_per_direction(0);
  int n_node_y = n_nodes_per_direction(1);
  double dx = len_x_ / n_elems_per_direction(0);
  for (int i = 0; i < n_node_x; i++) {
    double x_coord = i * dx;
    for (int j = 0; j < n_node_y; j++) {
      double y_coord = GetY(i, j);

      LocalNodeID nodeID = n_node_x * j + i;

      NODE* pNode = new NODE();
      node_array_[nodeID] = pNode;
      pNode->setCoor(x_coord, y_coord);

      if (i == 0) { pNode->addIndicatorNum(1); }
      if (i == n_node_x - 1) { pNode->addIndicatorNum(2); }
      if (j == 0) { pNode->addIndicatorNum(3); }
      if (j == n_node_y - 1) { pNode->addIndicatorNum(4); }
    }
  }
}

/*
 void GRIDBox2D::SetBndrIndicators(InputData* pIdata){

 for(int nodeID=1;nodeID<=this->n_nodes();nodeID++)
 {
 int inds[4];
 int nInd=0;
 double x = this->GetCoord (nodeID, 1);

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
 double y = this->GetCoord (nodeID, 2);
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


 this->GetNode(nodeID)->setIndicator(nInd,inds);
 }
 }
 */

}  // namespace TALYFEMLIB
