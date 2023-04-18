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
#include <talyfem/fem/preallocator_original.h>

#include <time.h>

#include <algorithm>
#include <vector>  // for std::vector
#include <set>

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/elem.h>

namespace TALYFEMLIB {

PreallocatorOriginal::PreallocatorOriginal() {
  p_grid_ = NULL;
  n_dof_ = 0;
  prealloc_n_periodic_vars = 0;
  prealloc_n_total_vars = 0;
  prealloc_n_periodic_bounds = 0;
  dnz = NULL;
  onz = NULL;
}

PreallocatorOriginal::~PreallocatorOriginal() {
  Destroy();
}

void PreallocatorOriginal::Destroy() {
  if (dnz != NULL)
    PetscFree(dnz);
  if (onz != NULL)
    PetscFree(onz);
}

void PreallocatorOriginal::redim(GRID* pGrid, int ndof) {
  p_grid_ = pGrid;
  n_dof_ = ndof;
}

void PreallocatorOriginal::PresetPeriodicData(int n_per_bounds, int n_per_vars,
                                int n_total_vars) {
  prealloc_n_periodic_bounds = n_per_bounds;
  prealloc_n_periodic_vars = n_per_vars;
  prealloc_n_total_vars = n_total_vars;
  if (prealloc_n_periodic_bounds == 0 && prealloc_n_periodic_vars > 0) {
    throw TALYException()
        << "Setting periodic variables with no periodic boundaries!";
  }
  if (prealloc_n_periodic_vars == 0 && prealloc_n_periodic_bounds > 0) {
    PrintWarning("Setting periodic boundaries with no periodic values.");
  }
}

const PetscInt* PreallocatorOriginal::get_dnz() const {
  return dnz;
}

const PetscInt* PreallocatorOriginal::get_onz() const {
  return onz;
}

void PreallocatorOriginal::calc_preallocation(Mat& mat, PetscInt bs) {
  if (bs == -1) {
    bs = 1;
  }

  // These are used to track nodes that are in an element that on one
  // or more borders.
  // We keep track of these because they may be impacted by periodic
  // boundaries later. If they are, that will change the number of
  // diagonal and/or off diagonal entries which can invalidate the
  // preallocation and create errors. Any element on a periodic boundary
  // can have all of its nodes affected. To avoid this, we flag these
  // nodes and assign the maximum possible memory to these. The excess
  // memory ensures that there will always be enough space for entries.
  // the sets below are named based on the number of borders their nodes
  // are on.
  std::set < PetscInt > border1_elem_nodes;
  std::set < PetscInt > border2_elem_nodes;
  std::set < PetscInt > border3_elem_nodes;

  // share own nodes numbers among processors
  PetscInt nOwnNodes = p_grid_->GetNumOwnedNodes();
  PetscInt *nAllOwns;  // list of number of own nodes on each processor
  PetscInt ownNodesStart = 0, ownNodesEnd;
  PetscMalloc(GetMPISize() * sizeof(PetscInt), &nAllOwns);
  MPI_Allgather(&nOwnNodes, 1, MPI_TALYFEM_INT, nAllOwns, 1, MPI_TALYFEM_INT,
                PETSC_COMM_WORLD);

  for (int i = 0; i < GetMPIRank(); i++) {
    ownNodesStart += nAllOwns[i];
  }
  ownNodesEnd = ownNodesStart + nOwnNodes;

  const int num_rows = n_dof_ * p_grid_->n_nodes();
  std::vector < std::vector<PetscInt> > diagtermno(num_rows);
  std::vector < std::vector<PetscInt> > offdiagtermno(num_rows);
  PetscMalloc(n_dof_ * nOwnNodes * sizeof(PetscInt) / bs, &dnz);
  PetscMalloc(n_dof_ * nOwnNodes * sizeof(PetscInt) / bs, &onz);

  int AroundNeigborD = static_cast<int>(pow(3.0, p_grid_->nsd()));
  int AroundNeigborO = static_cast<int>(pow(3.0, p_grid_->nsd()));
  AroundNeigborD *= (3 * this->p_grid_->basis_order() - 2);
  AroundNeigborO *= (3 * this->p_grid_->basis_order() - 2);
  if (p_grid_->nsd() == 3 && this->p_grid_->basis_order() > 1) {
    AroundNeigborD *= this->p_grid_->basis_order();
    AroundNeigborO *= this->p_grid_->basis_order();
  }
  // periodic boundaries with some variables periodic and others not can
  // lead to exceeding the number of unknowns given above. Consider the
  // mesh below with the X direction periodic and the Y direction
  // non periodic. The original mesh is on the left, the mesh after
  // applying periodic boundaries to the nodes is on the right.
  // 1--2--3--4                           1--2--3--1
  // |  |  |  |        Y                  |  |  |  |
  // 5--6--7--8        ^                  5--6--7--5
  // |  |  |  |        |                  |  |  |  |
  // 9--10-11-12       |                  9--10-11-9
  // |  |  |  |        |                  |  |  |  |
  // 13-14-15-16       --------->  X      13-14-15-13
  //
  // The system has two unknowns, A and B.
  // First consider the case where both A and B are periodic variables
  // Look at node 5, it will have 18 unknowns (agrees with calculation
  // above).
  // Six come from the A variable on nodes 1,2,5,6,9,10.
  // Three come from the A variable from the periodic remapping of node
  // 5 to the other side, which connects it to the three nodes: 3,7,11.
  // In addition to those nine unknowns, there are nine more from
  // variable B.
  //
  // Now consider A to be periodic (uses mesh on right) and B to be
  // non-periodic (uses mesh on left). Again look at node 5.
  // As before, it has six unknowns from A on nodes 1,2,5,6,8,10.
  // It has six unknowns from B on nodes 1,2,5,6,8,10, again as before
  // As in the previous case, the periodic remapping of 5 gives three
  // unknowns for A from nodes 3,7,11 and three unknowns from B for the
  // same nodes.
  // This gives 18 unknowns, as above.
  // however, since B is not periodic, it uses the mesh on the left and
  // has different values on nodes 4,8,12 which are also added during the
  // assembly. This gives a total of **21** unknowns... three more than
  // the "maximum" value.
  //
  // This can only happen when there is more than one degree of freedom
  // and at least one is periodic. Every non-periodic degree of freedom
  // will add 3^(nsd-1) unknowns.
  //
  // When both X and Y are periodic, the corner nodes add 12 unknown for
  // each non-periodic degree of freedom. A similar process gives the
  // number of additional unknowns for different dimensionalities.
  //
  // TODO: these additional values may only hold for 1st order basis
  // functions with block size 1.

  // set the maximum unknown number for periodic systems.
  // **this includes the n_dof_ values**
  // The values are named by the number of boundaries a node is on.
  // ideally this would only include periodic boundaries
  int AroundNeigborO_periodic1 = AroundNeigborO * n_dof_;
  int AroundNeigborD_periodic1 = AroundNeigborD * n_dof_;
  int AroundNeigborO_periodic2 = AroundNeigborO * n_dof_;
  int AroundNeigborD_periodic2 = AroundNeigborD * n_dof_;
  int AroundNeigborO_periodic3 = AroundNeigborO * n_dof_;
  int AroundNeigborD_periodic3 = AroundNeigborD * n_dof_;
  int non_per_vars = prealloc_n_total_vars - prealloc_n_periodic_vars;

  if (non_per_vars > 0 && prealloc_n_periodic_bounds == 1) {
    // only 1 boundary is periodic, so the maximum memory to allocate
    // doesn't depend on how many boundaries the node is on.
    AroundNeigborO_periodic1 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborD_periodic1 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborO_periodic2 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborD_periodic2 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborO_periodic3 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborD_periodic3 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
  }
  if (non_per_vars > 0 && prealloc_n_periodic_bounds == 2) {
    // 2 boundaries are periodic. In 3D, the maximum memory to allocate
    // for a corner (on three boundaries) is the same as an edge
    // (on 2 boundaries).
    AroundNeigborO_periodic1 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborD_periodic1 += static_cast<int>(pow(3.0, (p_grid_->nsd() - 1)))
        * non_per_vars;
    AroundNeigborO_periodic2 +=
        4 * static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborD_periodic2 +=
        4 * static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborO_periodic3 +=
        4 * static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborD_periodic3 +=
        4 * static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
  }
  if (non_per_vars > 0 && prealloc_n_periodic_bounds == 3) {
    // 3 boundaries are periodic. The maximum memory will be different
    // for nodes on different numbers of boundaries.
    AroundNeigborO_periodic1 +=
        static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborD_periodic1 +=
        static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborO_periodic2 +=
        4 * static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborD_periodic2 +=
        4 * static_cast<int>(pow(3.0, (p_grid_->nsd() - 1))) * non_per_vars;
    AroundNeigborO_periodic3 += 56 * non_per_vars;
    AroundNeigborD_periodic3 += 56 * non_per_vars;
  }
  // make sure we don't exceed the number of unknowns for the
  // diagonal terms
  if (AroundNeigborD_periodic1 > p_grid_->GetNumOwnedNodes() * n_dof_)
    AroundNeigborD_periodic1 = p_grid_->GetNumOwnedNodes() * n_dof_;
  if (AroundNeigborD_periodic2 > p_grid_->GetNumOwnedNodes() * n_dof_)
    AroundNeigborD_periodic2 = p_grid_->GetNumOwnedNodes() * n_dof_;
  if (AroundNeigborD_periodic3 > p_grid_->GetNumOwnedNodes() * n_dof_)
    AroundNeigborD_periodic3 = p_grid_->GetNumOwnedNodes() * n_dof_;

  if (AroundNeigborD > p_grid_->GetNumOwnedNodes())
    AroundNeigborD = p_grid_->GetNumOwnedNodes();

  PetscInt *lclNodesTolclRows;
  PetscMalloc(p_grid_->n_nodes() * sizeof(PetscInt), &lclNodesTolclRows);
  for (int i = 0, j = 0; i < p_grid_->n_nodes(); i++) {
    if (p_grid_->solution_map(i) >= ownNodesStart
        && p_grid_->solution_map(i) < ownNodesEnd) {
      lclNodesTolclRows[i] = j++;
    } else {
      lclNodesTolclRows[i] = -1;
    }
  }

  // skeletal assembly of matrix
  for (int elmID = 0; elmID < p_grid_->n_elements(); elmID++) {
    const ELEM* elem = p_grid_->GetElm(elmID);
    //~ if (!IsMyElement(elmID)) continue;
    int vertexn = elem->n_nodes();
    ZEROARRAY < PetscInt > vertexArray;
    vertexArray.redim(vertexn);
    ZEROARRAY < PetscInt > vertexArrayLoc;
    vertexArrayLoc.redim(vertexn);
    PetscInt* vertex_ptr = vertexArray.data();
    PetscInt* vertex_loc_ptr = vertexArrayLoc.data();

    bool elem_on_1border = false;  // element on 1 border
    bool elem_on_2border = false;  // element on 2 borders
    bool elem_on_3border = false;  // element on 3 borders
    for (int i = 0; i < elem->n_nodes(); i++) {
      PetscInt lid = elem->node_id_array(i);
      PetscInt gid = p_grid_->solution_map(elem->node_id_array(i));

      vertex_ptr[i] = gid;
      vertex_loc_ptr[i] = lid;

      // check if any nodes are border nodes
      if (prealloc_n_periodic_bounds > 0) {
        int boundary_count = p_grid_->node_array_[lid]->getIndicatorNo();
        if (boundary_count >= 3) {
          elem_on_3border = true;
        } else if (boundary_count == 2) {
          elem_on_2border = true;
        } else if (boundary_count == 1) {
          elem_on_1border = true;
        }
      }
    }

    // if any node in the element is a border node, then all of the
    // nodes are near border nodes, so store them as such
    if (prealloc_n_periodic_bounds > 0) {
      if (elem_on_3border) {
        for (int i = 0; i < elem->n_nodes(); i++) {
          border3_elem_nodes.insert(elem->node_id_array(i));
        }
      } else if (elem_on_2border) {
        for (int i = 0; i < elem->n_nodes(); i++) {
          border2_elem_nodes.insert(elem->node_id_array(i));
        }
      } else if (elem_on_1border) {
        for (int i = 0; i < elem->n_nodes(); i++) {
          border1_elem_nodes.insert(elem->node_id_array(i));
        }
      }
    }

    for (int i = 0; i < elem->n_nodes(); i++) {
      if (lclNodesTolclRows[vertex_loc_ptr[i]] != -1) {
        for (int j = 0; j < elem->n_nodes(); j++) {
          for (int k = 0; k < n_dof_; k++) {
            if (lclNodesTolclRows[vertex_loc_ptr[j]] != -1) {
              std::vector<PetscInt>::iterator it = find(
                  diagtermno[vertex_loc_ptr[i] * n_dof_ + k].begin(),
                  diagtermno[vertex_loc_ptr[i] * n_dof_ + k].end(),
                  vertex_loc_ptr[j] * n_dof_ + k);
              if (it == diagtermno[vertex_loc_ptr[i] * n_dof_ + k].end()) {
                diagtermno[vertex_loc_ptr[i] * n_dof_ + k].push_back(
                    vertex_loc_ptr[j] * n_dof_ + k);
              }
            } else {
              std::vector<PetscInt>::iterator it = find(
                  offdiagtermno[vertex_loc_ptr[i] * n_dof_ + k].begin(),
                  offdiagtermno[vertex_loc_ptr[i] * n_dof_ + k].end(),
                  vertex_loc_ptr[j] * n_dof_ + k);
              if (it == offdiagtermno[vertex_loc_ptr[i] * n_dof_ + k].end()) {
                offdiagtermno[vertex_loc_ptr[i] * n_dof_ + k].push_back(
                    vertex_loc_ptr[j] * n_dof_ + k);
              }
            }
          }
        }
      }
    }
  }

  // write number of nonzeroes per row using data
  // collected during skeletal assembly of matrix
  for (int i = 0, j = 0; i < p_grid_->n_nodes(); i++) {
    if (lclNodesTolclRows[i] != -1) {
      // TODO: this needs to be properly fixed.
      // this is a temporary work around to the problem of preallocate
      // being unaware of the periodic boundaries. Without knowing
      // about the boundaries, it is unable to correctly allocate
      // memory. For now, we assume that any nodes near the boundary
      // of the domain will have the maximum number of neighbors.
      // Otherwise, we use the memory that was determined above.
      // This overestimates the memory needed, but ensures there is
      // always enough allocated.
      if (border3_elem_nodes.count(i) > 0) {
        for (int k = 0; k < n_dof_; k++) {
          // do *not* multiply by n_dof_ (was already done above)
          dnz[(j * n_dof_ + k) / bs] = AroundNeigborD_periodic3 / bs;
          onz[(j * n_dof_ + k) / bs] = AroundNeigborO_periodic3 / bs;
        }
      } else if (border2_elem_nodes.count(i) > 0) {
        for (int k = 0; k < n_dof_; k++) {
          // do *not* multiply by n_dof_ (was already done above)
          dnz[(j * n_dof_ + k) / bs] = AroundNeigborD_periodic2 / bs;
          onz[(j * n_dof_ + k) / bs] = AroundNeigborO_periodic2 / bs;
        }
      } else if (border1_elem_nodes.count(i) > 0) {
        for (int k = 0; k < n_dof_; k++) {
          // do *not* multiply by n_dof_ (was already done above)
          dnz[(j * n_dof_ + k) / bs] = AroundNeigborD_periodic1 / bs;
          onz[(j * n_dof_ + k) / bs] = AroundNeigborO_periodic1 / bs;
        }
      } else {
        for (int k = 0; k < n_dof_; k++) {
          dnz[(j * n_dof_ + k) / bs] = diagtermno[i * n_dof_ + k].size()
              * n_dof_ / bs;
          onz[(j * n_dof_ + k) / bs] = offdiagtermno[i * n_dof_ + k].size()
              * n_dof_ / bs;
        }
      }
      j++;
    }
  }

  // nodes on boundary have also entries from other processes (ghost elements)
  for (int i = 0; i < p_grid_->n_shared_nodes_(GetMPIRank()); i++) {
    if (lclNodesTolclRows[p_grid_->shared_nodes_(i)] != -1) {
      for (int k = 0; k < n_dof_; k++) {
        onz[(lclNodesTolclRows[p_grid_->shared_nodes_(i)] * n_dof_ + k) / bs] =
            AroundNeigborO * n_dof_ / bs;
        dnz[(lclNodesTolclRows[p_grid_->shared_nodes_(i)] * n_dof_ + k) / bs] =
            AroundNeigborD * n_dof_ / bs;
      }
    }
  }

  /*
   for (int i = 0; i < p_grid_->n_nodes(); i++) {
   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d] lclNodesTolclRows[%d]: %d\n",
   p_grid_->grid_id(), i, lclNodesTolclRows[i]);
   }
   PetscSynchronizedFlush(PETSC_COMM_WORLD);
   MPI_Barrier(PETSC_COMM_WORLD);

   for (int i = 0, j = 0; i < p_grid_->n_nodes(); i++) { //use j to skip not owned
   if (lclNodesTolclRows[i] != -1) {
   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d] %d(%d:%d) dnz: %d onz %d\n",
   p_grid_->grid_id(), j, p_grid_->solution_map(i), p_grid_->physical_map(i), dnz[j], onz[j]);
   j++;
   }
   }
   PetscSynchronizedFlush(PETSC_COMM_WORLD);
   MPI_Barrier(PETSC_COMM_WORLD);

   for (int i = 0; i < p_grid_->n_nodes(); i++) {
   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d] solution_map(%d):%d phyMap(%d):%d\n",
   p_grid_->grid_id(), i, p_grid_->solution_map(i), i, p_grid_->physical_map(i));
   }
   PetscSynchronizedFlush(PETSC_COMM_WORLD);
   MPI_Barrier(PETSC_COMM_WORLD);
   */

  PetscFree(nAllOwns);
  PetscFree(lclNodesTolclRows);
}

}  // namespace TALYFEMLIB

