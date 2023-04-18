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
#include <talyfem/domain_decomposition/cross_domain_communicator.h>

#include <talyfem/grid/shareinfo.h>
#include <talyfem/grid/grid_types/grid.h>


namespace TALYFEMLIB {

void CrossDomainCommunicator::Initialize(GRID* p_grid) {
  CleanCommData();  // remove any existing data
  ZEROARRAY<PetscInt> send_indices_source;  // indices in local data array
  ZEROARRAY<PetscInt> send_indices;  // indices of where to put data to send
  ZEROARRAY<PetscInt> recv_indices;  // where to get received data from
  ZEROARRAY<PetscInt> recv_indices_destination;  // where to put recevied data
  IS is_recv_destination;  // index set formed from recv_indices_destination
  IS is_recv;  // index set formed from recv_indices
  IS is_send_source;  // indes set formed from send_indices_source
  IS is_send;  // indes set formed from send_indices

  const int grid_id = p_grid->grid_id();  // id of this process
  const int n_grids = p_grid->n_subgrids();
  const int n_shared_nodes = p_grid->n_shared_nodes_(grid_id);

  send_indices.redim(n_shared_nodes);
  send_indices_source.fill_sequence(n_shared_nodes);  // 0, 1, 2, 3, ...

  // index where data for grid "i" starts
  ZEROARRAY<PetscInt> index_of_grid_data;
  index_of_grid_data.redim(n_grids);
  index_of_grid_data(0) = 0;  // first process starts at the beginning
  for (int i_grid = 0; i_grid < n_grids - 1; i_grid++) {
    // location of the next grid's data starts at the end of the current
    // grid's data
    index_of_grid_data(i_grid + 1) = index_of_grid_data(i_grid) +
                                     p_grid->n_shared_nodes_(i_grid);
  }
  // where to put data items in the transfer matrix
  for (int offset = 0; offset < n_shared_nodes; offset++) {
    // first term is the starting location of this processor's data in the
    // global transfer array
    send_indices(offset) = index_of_grid_data(grid_id) + offset;
  }
  send_node_indices_ = &(p_grid->shared_nodes_);

  // number of nodes we'll get from other processes
  int n_nodes_from_others = p_grid->CalcNumNodesOthers();

  recv_indices.redim(n_nodes_from_others);
  recv_indices_destination.redim(n_nodes_from_others);

  recv_node_indices_.redim(n_nodes_from_others);
  int index = 0;
  for (LocalNodeID A = 0; A < p_grid->n_nodes(); A++) {
    if (!(p_grid->node_belong_(A).is_owned)) {  // grab data for unowned nodes
      int unowned_grid_id = p_grid->node_belong_(A).share_data(0).grid_id();
      int cmuID = index_of_grid_data(unowned_grid_id)
          + p_grid->node_belong_(A).share_data(0).comm_id();
      recv_indices(index) = cmuID;
      recv_indices_destination(index) = index;
      recv_node_indices_(index) = A;
      index++;
    }
  }

  // Create Vector for communication
  // this is the total size of the transfer vector. The first term is the
  // starting offset of the data for the last process and the second term
  // is the number of data items for that process.
  int total_shared_nodes = index_of_grid_data(n_grids - 1) +
                           p_grid->n_shared_nodes_(n_grids - 1);

  VecCreateMPI(PETSC_COMM_WORLD, n_shared_nodes, total_shared_nodes,
                        &transfer_vec_);
  ISCreateGeneral(PETSC_COMM_SELF, n_shared_nodes,
                           send_indices_source.data(), PETSC_COPY_VALUES,
                           &is_send_source);
  ISCreateGeneral(PETSC_COMM_SELF, n_shared_nodes,
                           send_indices.data(), PETSC_COPY_VALUES, &is_send);
  ISCreateGeneral(PETSC_COMM_SELF, n_nodes_from_others,
                           recv_indices.data(), PETSC_COPY_VALUES,
                           &is_recv_destination);
  ISCreateGeneral(PETSC_COMM_SELF, n_nodes_from_others,
                           recv_indices_destination.data(), PETSC_COPY_VALUES,
                           &is_recv);

  send_data_.redim(n_shared_nodes);  // storage space for outgoing data
  recv_data_.redim(n_nodes_from_others);  // storage space for incoming data

  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, send_data_.size(),
                                 send_data_.data(), &send_data_vec_);
  VecScatterCreate(send_data_vec_, is_send_source, transfer_vec_,
                            is_send, &scatter_send_);
  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, recv_data_.size(),
                                 recv_data_.data(), &recv_data_vec_);
  VecScatterCreate(transfer_vec_, is_recv_destination,
                            recv_data_vec_, is_recv, &scatter_recv_);
  ISDestroy(&is_send_source);
  ISDestroy(&is_send);
  ISDestroy(&is_recv_destination);
  ISDestroy(&is_recv);

  initialized_ = true;
}

void CrossDomainCommunicator::SendRecvData() {
  // send values
  VecScatterBegin(scatter_send_, send_data_vec_, transfer_vec_,
                           INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scatter_send_, send_data_vec_, transfer_vec_,
                         INSERT_VALUES, SCATTER_FORWARD);
  // receive values
  VecScatterBegin(scatter_recv_, transfer_vec_, recv_data_vec_,
                           INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scatter_recv_, transfer_vec_, recv_data_vec_,
                         INSERT_VALUES, SCATTER_FORWARD);
}

void CrossDomainCommunicator::CleanCommData() {
  if (!initialized_) { return; }
  VecDestroy(&transfer_vec_);
  VecDestroy(&send_data_vec_);
  VecDestroy(&recv_data_vec_);
  VecScatterDestroy(&scatter_send_);
  VecScatterDestroy(&scatter_recv_);
}

}  // namespace TALYFEMLIB
