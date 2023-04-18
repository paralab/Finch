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
#include <talyfem/fem/periodic_bounds.h>

#include <talyfem/common/exceptions.h>  // for throwing exception
#include <talyfem/grid/grid-types.h>
#include <talyfem/domain_decomposition/mesh_partition.h>  // for MPI_TALYFEM_INT

namespace TALYFEMLIB {

PeriodicBounds::PeriodicBounds(GRID *p_grid, int* ind, int indno) {
  p_grid_ = p_grid;

  is_periodic_ = false;

  int rank = GetMPIRank();
  int size = GetMPISize();

  // if there are no periodic bounds, do nothing else
  if (indno <= 0) { return; }

  is_periodic_ = true;
  // make note of which boundary indicators are periodic so we can easily
  // determine node parters later.
  this->is_boundary_periodic_.redim(p_grid_->nsd() * 2);
  this->is_boundary_periodic_.fill(false);
  for (int i = 0; i < indno; i++) {
    this->is_boundary_periodic_.set(ind[i]-2, true);
    this->is_boundary_periodic_.set(ind[i]-1, true);
  }

  this->n_periodic_bounds_ = indno;  // number of periodic boundaries

  // set the partners for the boundaries
  is_node_periodic_.redim(this->p_grid_->n_nodes());
  is_node_periodic_.fill(false);
  for (LocalNodeID oldID = 0; oldID < this->p_grid_->n_nodes(); oldID++) {
    PhysicalNodeID newID = NewNodeID(oldID);
    if (p_grid_->parallel_type_ == kWithDomainDecomp) {
      // newID is a *global* id, oldID is a *local* one. we need to compare
      // the *global* version of oldID with newID. Otherwise, we may miss a
      // periodic node if newID == oldID.
      if (newID != p_grid_->physical_map(oldID)) {
        is_node_periodic_.set(oldID, true);
        pbc_partners_[oldID] = newID;
      }
    } else {
      // newID is a *global* id, oldID is both the local and global id, no
      // need to try to get a global value in this case.
      // without domain decomposition, pbc_partners_ and pbc_sol_partners_
      // should be the same.
      if (newID != oldID) {
        is_node_periodic_.set(oldID, true);
        pbc_partners_[oldID] = newID;
        pbc_sol_partners_[oldID] = newID;
      }
    }
  }

  // for domain decomposition, the partner ID must be the global index
  // of the variable within the petsc solution vector. To find this, we
  // need to find the global location of the node within the solution array.
  // This is given by p_grid_->solution_map(localID).
  // However, this will only work if we use the localID that corresponds to
  // the node that is a periodic partner to the desired node. In general
  // this will not be on the same processor as the node we are replacing. So
  // we need to have the processor that owns the partner node call
  // p_grid_->solution_map on that node and tell us the result. Since we don't
  // know which processor owns which node, we need to have every processor check
  // for the node. To do this:
  // 1) we make a list, local_nodes, of the global nodeIDs we want the
  //    solution indices for.
  // 2) this list is sent to every process and stored in the search_nodes
  //    array. Each process checks whether it has each node in the array
  //    by checking if the globalID is in p_grid_->physical_map. If the desired
  //    node is found on a process, that process runs p_grid_->solution_map on
  //    it and stores the result in search_nodes. All nodes that were not found
  //    by a process are set to -1 in search_nodes.
  // 3) The orginal node gathers the found_node arrays from all other
  //    processes and uses the solution indices to set the idx_from_ array
  //    for domain decomposition and the pbc_sol_partners_ map for use during
  //    assembly.
  // 4) This sequence is repeated for each running process.
  //

  if (p_grid_->parallel_type_ == kWithDomainDecomp) {
    // find number of items needed by each process
    // number of desired partners for this process
    int n_desired_nodes = pbc_partners_.size();
    int *search_counts = new int[size];  // # of ids each process needs.
    MPI_Allgather(&n_desired_nodes, 1, MPI_INT, search_counts, 1,
                  MPI_INT, MPI_COMM_WORLD);

    // make arrays for transfering data. The length is equal to the
    // number of periodic target nodes we have.
    LocalNodeID *local_indices = new LocalNodeID[n_desired_nodes];
    // global ids of nodes we want to find
    PhysicalNodeID *desired_nodes = new PhysicalNodeID[n_desired_nodes];
    // gathered results of partner search
    SolutionNodeID *result_nodes = new SolutionNodeID[n_desired_nodes];

    int i_idx = 0;
    PeriodicPhysicalMap::iterator iter;
    for (iter = pbc_partners_.begin(); iter != pbc_partners_.end(); ++iter) {
      local_indices[i_idx] = iter->first;
      desired_nodes[i_idx] = iter->second;
      i_idx++;
    }

    // loop over processors and find map values for periodic boundaries
    for (int i_rank = 0; i_rank < size; i_rank++) {
      int array_size = search_counts[i_rank];  // size of array to recieve
      if (array_size == 0) { continue; }  // nothing needed, skip

      // nodes to look for
      PhysicalNodeID* search_nodes = new PhysicalNodeID[array_size];
      if (rank == i_rank) {
        // process we're looking for needs to copy list of nodes so it
        // can be broadcast to the other nodes.
        size_t mem_size = sizeof(desired_nodes[0]) * array_size;
        memcpy(search_nodes, desired_nodes, mem_size);
      }

      // send list of nodes we want
      MPI_Bcast(search_nodes, array_size, MPI_TALYFEM_INT,
                i_rank, MPI_COMM_WORLD);

      // Each processor loops over the received values. If it has the
      // value in its physical_map, it will find the global solution id and
      // set it in the return array. Otherwise, it sets -1 to indicate
      // it was not found.
      // Each time a node is found, it is added to the list of nodes that
      // need to be sent.
      for (int i = 0; i < array_size; i++) {
        LocalNodeID local_id = p_grid_->GetPhysicalMap()->find(search_nodes[i]);
        if (local_id != -1) {
          SolutionNodeID global_id = p_grid_->solution_map(local_id);
          search_nodes[i] = global_id;  // index of solution array
        } else {  // not in our list of nodes
          search_nodes[i] = -1;
        }
      }

      // Gather all the values to the main process. Every node that was not
      // found will be set to -1 in the array. Nodes that were found will
      // be set to the appropriate index value. Multiple processes may have
      // the desired node, but they will all have the same index value for
      // it. Thus there are only two possible values for each array index:
      // -1 and our desired index. Since the desired index must be >= 0,
      // MPI_MAX yields the correct result.
      MPI_Reduce(search_nodes, result_nodes, array_size, MPI_TALYFEM_INT,
                 MPI_MAX, i_rank, MPI_COMM_WORLD);

      if (i_rank == rank) {
        for (int k = 0; k < array_size; k++) {
          // set the correct solution partners and values for the domain
          // decomposition 'from_' vector
          LocalNodeID src_node = local_indices[k];
          SolutionNodeID dst_node = result_nodes[k];
          pbc_sol_partners_[src_node] = dst_node;
        }
      }
      delete [] search_nodes;
    }
    delete [] local_indices;
    delete [] desired_nodes;
    delete [] result_nodes;
    delete [] search_counts;
  }
}

PhysicalNodeID PeriodicBounds::NewNodeID(LocalNodeID lclnodeID) const {
  PhysicalNodeID nodeID = lclnodeID;

  if (p_grid_->parallel_type_ == kWithDomainDecomp) {  // find global ID
    nodeID = p_grid_->physical_map(nodeID);
  }

  PhysicalNodeID newID = nodeID;
  switch (this->p_grid_->grid_type()) {
    case kGrid3dBox: {  // 3D Box Grid
      GRIDBox3D* p_grid = static_cast<GRIDBox3D*>(this->p_grid_);
      // Check and apply boundary conditions for each dimension
      // individually. This correctly handles cases with mixed boundary
      // conditions without explicitly treating every possible combination.
      if (this->is_boundary_periodic_.get(1) &&  // TODO should this match
                                                 // with BoNode? (it used to)
          this->p_grid_->BoNode(lclnodeID, 2)) {
        newID -= p_grid->n_elems_per_direction(0);
      }
      if (this->is_boundary_periodic_.get(3) &&
          this->p_grid_->BoNode(lclnodeID, 4)) {
        newID -= p_grid->n_elems_per_direction(1) *
                 (p_grid->n_elems_per_direction(0) + 1);
      }
      if (this->is_boundary_periodic_.get(5) &&
          this->p_grid_->BoNode(lclnodeID, 6)) {
        newID -= p_grid->n_elems_per_direction(2) *
                 (p_grid->n_elems_per_direction(1) + 1) *
                 (p_grid->n_elems_per_direction(0) + 1);
      }
      return newID;
    }
    case kGrid2dBox: {  // 2D Box Grid
      GRIDBox2D* p_grid = static_cast<GRIDBox2D*>(this->p_grid_);
      if (this->is_boundary_periodic_.get(1) &&
          this->p_grid_->BoNode(lclnodeID, 2)) {
        newID -= p_grid->n_elems_per_direction(0);
      }
      if (this->is_boundary_periodic_.get(3) &&
          this->p_grid_->BoNode(lclnodeID, 4)) {
        newID -= p_grid->n_elems_per_direction(1) *
                 (p_grid->n_elems_per_direction(0) + 1);
      }
      return newID;
    }
    case kGrid1d: {  // 1D Grid
      GRID1D* p_grid = static_cast<GRID1D*>(this->p_grid_);
      if (this->is_boundary_periodic_.get(1) &&
          this->p_grid_->BoNode(lclnodeID, 2)) {
        newID -= p_grid->n_elements();
      }
      return newID;
    }
    default: {
      throw NotImplementedException() << "Unsupported gridType in NewNodeID.";
    }
  }
}

}  // namespace TALYFEMLIB
