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
#ifndef FEM_PERIODIC_EXCHANGER_H_
#define FEM_PERIODIC_EXCHANGER_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif

#include <vector>  // for std::vector

#include <talyfem/grid/grid_types/grid.h>  // for items from GRID class
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

class PeriodicBounds;

/**
 * Class to exchange data values across periodic boundaries.
 *
 * This is used to ensure that periodic node partners have the same data.
 * For example, if initial conditions for a solver are determined randomly
 * then a node and its periodic partner may have a different value. This
 * will lead to inconsistencies in the data. This function exchanges the
 * data values across periodic boundaries. It handles domain decomposition
 * where the periodic nodes might be on different processes.
 * This does *not* do inter-process exchange across process boundaries
 * within the interior of the domain.
 *
 * Note that this may not be reproducible. When multiple processes share the
 * same node, there is no mechanism to determine which process's value will
 * be used as the actual value. However, it is guaranteed that after the
 * exchange each process will have the same value on the shared nodes.
 *
 * To reduce memory usage, this class is designed to use an external PETSc Vec
 * to store values during exchange. This has two important implications:
 * 1) Prior to initializing this class, the supplied Vec must already be created
 * and have a length equal to or greater than the number of physical nodes in
 * the system.
 * 2) The initialization and exchange process will overwrite the contents of
 * this Vec.
 *
 * Typical usage of this class is as follows:
 * // Create object using previously created PeriodicBounds object and
 * // previously created and initialized array to transfer data
 * PeriodicExchanger exchanger(periodic_bounds, array);
 * exchanger.Initialize();  // initialize object
 *
 * // to exchange data, we take date that we want to sync, in this case data
 * // from the 'data_array' array, store it in the send buffer, then call the
 * // exchange buffer data to send needs to be stored in the send buffer
 * PetscScalar *send_buffer = exchanger.send_buffer();
 * const int *send_indices = exchanger.send_indices();
 * for (int i = 0; i < exchanger.n_send_values(); i++) {
 *   int node_id = send_indices[i];  // index of ith source node
 *   // store the data to exchange in the send buffer
 *   send_buffer[i] = data_array[node_id];
 * }
 *
 * // do the actual exchange
 * exchanger.DoPBCExchange();
 *
 * // pull the final exchanged values and store them in the proper places
 * PetscScalar *recv_buffer = exchanger.recv_buffer();
 * const int *recv_indices = exchanger.recv_indices();
 * for (int i = 0; i < exchanger.n_recv_values(); i++) {
 *   int node_id = recv_indices[i];  // index of ith target node
 *   // store the data from exchange in the proper location
 *   data_array[node_id] = recv_buffer[i];
 * }
 *
 */
class PeriodicExchanger {
 public:
  /**
   * Construct the object using given arguments
   *
   * @param periodic_bounds PeriodicBounds object describing periodic system
   * @param transfer_array PETSc Vec to use in transferring data
   */
  PeriodicExchanger(PeriodicBounds *periodic_bounds, Vec *transfer_array);

  ~PeriodicExchanger();

  /**
   * Initialize the structures needed to allow exchange.
   *
   * Note: this will wipe out the contents of the array used as a transfer
   * array for this exchanger.
   */
  void Initialize();

  /**
   * Sends values across periodic boundaries
   *
   * See class comment for example usage
   */
  void DoPBCExchange();

  /**
   * Returns the number of values this process will receive.
   */
  inline int n_recv_values() const { return n_recv_values_; }

  /**
   * Returns the number of values this process will send.
   */
  inline int n_send_values() const { return n_send_values_; }

  /**
   * Returns a pointer to the receive buffer.
   */
  inline PetscScalar* recv_buffer() { return recv_buffer_; }

  /**
   * Returns a pointer to the send buffer.
   */
  inline PetscScalar* send_buffer() { return send_buffer_; }

  /**
   * Returns a pointer to the array of receiving node indices.
   */
  inline const LocalNodeID* recv_indices() const {
    return &recv_buffer_to_local_indices_[0];
  }

  /**
   * Returns a pointer to the array of sending node indices.
   */
  inline const LocalNodeID* send_indices() const {
    return &local_indices_to_send_[0];
  }

 private:
  /**
   * Create the data structures for receiving exchange data
   *
   * This entails determining which periodic values this process needs to
   * obtain from other processes.
   */
  void InitializeRecv();

  /**
   * Create the data structures for receiving exchange data
   */
  void InitializeSend();

  PetscScalar *send_buffer_;  ///< buffer for data to send
  PetscScalar *recv_buffer_;  ///< buffer for data to receive
  Vec *transfer_array_;  ///< Vec to use when transferring data
  VecScatter send_to_array_;  ///< scatter for sending data
  VecScatter recv_from_array_;  ///< scatter for receiving data
  Vec send_vec_;  ///< Vec to use when sending data
  Vec recv_vec_;  ///< Vec to use when receiving data
  ///< where to put data
  std::vector<LocalNodeID> recv_buffer_to_local_indices_;
  std::vector<LocalNodeID> local_indices_to_send_;  ///< where to get data from
  PeriodicBounds *periodic_bounds_;  ///< periodic boundary details
  bool is_initialized_;  ///< whether exchanger is initialized
  int n_recv_values_;  ///< number of values to receive
  int n_send_values_;  ///< number of values to send
};

}  // namespace TALYFEMLIB

#endif  // FEM_PERIODIC_EXCHANGER_H_
