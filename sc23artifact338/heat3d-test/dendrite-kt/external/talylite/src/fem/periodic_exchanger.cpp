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
#include <talyfem/fem/periodic_exchanger.h>

#include <vector>  // for std::vector

#include <talyfem/fem/periodic_bounds.h>  // class of member variable
#include <talyfem/grid/grid_types/grid.h>  // for items from GRID class
#include <talyfem/grid/node.h>  // for elem_id_ in Node class


namespace TALYFEMLIB {

PeriodicExchanger::PeriodicExchanger(PeriodicBounds *periodic_bounds,
                                     Vec *transfer_array)
    : send_buffer_(NULL),
      recv_buffer_(NULL),
      transfer_array_(transfer_array),
      send_to_array_(),
      recv_from_array_(),
      send_vec_(),
      recv_vec_(),
      recv_buffer_to_local_indices_(),
      local_indices_to_send_(),
      periodic_bounds_(periodic_bounds),
      is_initialized_(false),
      n_recv_values_(0),
      n_send_values_(0) { }

PeriodicExchanger::~PeriodicExchanger() {
  if (is_initialized_) {
    // clear receive data structures
    // recv_vec_ needs to be destroyed *before* deleting recv_buffer_
    VecDestroy(&recv_vec_);
    VecScatterDestroy(&recv_from_array_);
    delete [] recv_buffer_;

    // clear send data structures
    // send_vec_ needs to be destroyed *before* deleting send_buffer_
    VecDestroy(&send_vec_);
    VecScatterDestroy(&send_to_array_);
    delete [] send_buffer_;
  }
}

void PeriodicExchanger::Initialize() {
  // all of the initialization code assumes the system has not been
  // previously initialized.
  if (is_initialized_) {
    throw TALYException() << "Illegal attempt to re-initialize "
                          << "periodic exchanger.";
  }

  InitializeRecv();
  InitializeSend();
  is_initialized_ = true;
}

void PeriodicExchanger::DoPBCExchange() {
  // TODO: unsure whether this will work when there is nothing to scatter
  // send data to transfer array
  VecScatterBegin(send_to_array_, send_vec_, *transfer_array_,
                  INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(send_to_array_, send_vec_, *transfer_array_,
                INSERT_VALUES, SCATTER_FORWARD);
  // grab data from transfer array
  VecScatterBegin(recv_from_array_, *transfer_array_, recv_vec_,
                  INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(recv_from_array_, *transfer_array_, recv_vec_,
                INSERT_VALUES, SCATTER_FORWARD);
}

/**
 * To initialize the receive structures, we need to identify which of the values
 * this process owns are needed by another process. When doing an exchange, this
 * process will need to send these items to other processes.
 *
 * In order to find these values, we use the following strategy:
 * 1) Set all values in the transfer array to zero
 * 2) Every process sets the value of the transfer array to 1 at each
 *    index it needs data from.
 * 3) Every process then pulls the values from the transfer array for all
 *    nodes that it owns.
 * 4) The process stores the indices of all entries with a value of 1. These
 *    are the values that other processes need.
 */
void PeriodicExchanger::InitializeRecv() {
  std::vector<PhysicalNodeID> desired_indices;  // global indices we want
  PeriodicPhysicalMap::const_iterator iter;
  const PeriodicPhysicalMap *partners = periodic_bounds_->pbc_partners();
  for (iter = partners->begin(); iter!= partners->end(); ++iter) {
    // local index to put value in
    recv_buffer_to_local_indices_.push_back(iter->first);
    desired_indices.push_back(iter->second);
  }
  n_recv_values_ = desired_indices.size();

  // create index set of values we want to receive from other.
  IS desired_index_set;
  ISCreateGeneral(PETSC_COMM_SELF, n_recv_values_, &desired_indices[0],
                  PETSC_COPY_VALUES, &desired_index_set);
  // create buffer to hold values received.
  recv_buffer_ = new PetscScalar[n_recv_values_];
  // create Vec to hold receive values for scatter
  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n_recv_values_, recv_buffer_,
                        &recv_vec_);
  // create receive scatter
  VecScatterCreate(*transfer_array_, desired_index_set, recv_vec_, NULL,
                   &recv_from_array_);
  ISDestroy(&desired_index_set);  // no longer needed
}

void PeriodicExchanger::InitializeSend() {
  LocalNodeID n_nodes = periodic_bounds_->Grid()->n_nodes();

  // Step 1: set all values in array to zero
  VecSet(*transfer_array_, 0.0);

  // Step 2: set values we want from elsewhere to '1'.
  // we do this by doing a reverse scatter using the receive scatter
  // first set the values we want to store (1 for all values we want)
  VecSet(recv_vec_, 1.0);
  // send the values
  VecScatterBegin(recv_from_array_, recv_vec_, *transfer_array_,
                  INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(recv_from_array_, recv_vec_, *transfer_array_,
                INSERT_VALUES, SCATTER_REVERSE);

  // Step 3: pull the values for the nodes that we own.
  // first, make list of nodes we own (this is the physical ID of the node,
  // which corresponds to the output of pbc_partners()
  IS global_indices_in_process;
  if (periodic_bounds_->Grid()->parallel_type_ == kWithDomainDecomp) {
    ISCreateGeneral(PETSC_COMM_SELF, n_nodes,
                    periodic_bounds_->Grid()->GetPhysicalMap()->data(),
                    PETSC_USE_POINTER, &global_indices_in_process);
  } else {
    // construct list of nodes owned by this process
    std::vector<PetscInt> owned_node_list;
    for (LocalNodeID i = 0; i < periodic_bounds_->Grid()->n_nodes(); i++) {
      int e = periodic_bounds_->Grid()->node_array_[i]->elem_id_;
      if (periodic_bounds_->Grid()->IsMyElement(e)) {
        owned_node_list.push_back(i);
      }
    }
    // correct number of nodes to count only nodes this process uses
    n_nodes = owned_node_list.size();
    ISCreateGeneral(PETSC_COMM_SELF, owned_node_list.size(),
                    &owned_node_list[0],
                    PETSC_COPY_VALUES, &global_indices_in_process);
  }
  // Now create vector to store resulting values pulled from array
  Vec values_from_array;
  VecCreateSeq(PETSC_COMM_SELF, n_nodes, &values_from_array);

  // create scatter context to pull values from transfer_array
  VecScatter scatter;
  VecScatterCreate(*transfer_array_, global_indices_in_process,
                   values_from_array, NULL, &scatter);
  ISDestroy(&global_indices_in_process);  // no longer needed

  // take values from the transfer array and put them in the local Vec
  VecScatterBegin(scatter, *transfer_array_, values_from_array, INSERT_VALUES,
                  SCATTER_FORWARD);
  VecScatterEnd(scatter, *transfer_array_, values_from_array, INSERT_VALUES,
                SCATTER_FORWARD);
  VecScatterDestroy(&scatter);  // no longer needed

  // Step 4: identify any non-zero values. These are indices that some other
  // process is interested in obtaining during exchange.
  PetscScalar *values_arr;
  VecGetArray(values_from_array, &values_arr);
  for (LocalNodeID i = 0; i < n_nodes; i++) {
    if (values_arr[i] != 0.0) {
      local_indices_to_send_.push_back(i);
    }
  }
  VecRestoreArray(values_from_array, &values_arr);
  VecDestroy(&values_from_array);  // no longer needed

  // set up data structures for use during exchange
  n_send_values_ = local_indices_to_send_.size();
  // create buffer to hold values sent.
  send_buffer_ = new PetscScalar[n_send_values_];
  // create Vec to hold send values for scatter
  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n_send_values_, send_buffer_,
                        &send_vec_);

  // now we need a scatter to push the items we want and a
  // scatter to pull the ones other processes want.
  std::vector<PhysicalNodeID> target_indices(n_send_values_);
  if (periodic_bounds_->Grid()->parallel_type_ == kWithDomainDecomp) {
    for (int i = 0; i < n_send_values_; i++) {
      LocalNodeID local_index = local_indices_to_send_[i];
      target_indices[i] = periodic_bounds_->Grid()->physical_map(local_index);
    }
  } else {
    for (int i = 0; i < n_send_values_; i++) {
      LocalNodeID local_index = local_indices_to_send_[i];
      target_indices[i] = local_index;
    }
  }
  IS targets;
  ISCreateGeneral(PETSC_COMM_SELF, n_send_values_, &target_indices[0],
                  PETSC_COPY_VALUES, &targets);
  VecScatterCreate(send_vec_, NULL, *transfer_array_, targets,
                   &send_to_array_);
  ISDestroy(&targets);  // no longer needed
}

}  // namespace TALYFEMLIB
