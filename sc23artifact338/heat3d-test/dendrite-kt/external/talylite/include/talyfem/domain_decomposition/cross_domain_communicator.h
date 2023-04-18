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
#ifndef DOMAINDECOMPOSITION_CROSS_DOMAIN_COMMUNICATOR_H_
#define DOMAINDECOMPOSITION_CROSS_DOMAIN_COMMUNICATOR_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>  // for Petsc types
#endif

#include <talyfem/data_structures/zeroarray.h>


namespace TALYFEMLIB {

class GRID;

/**
 * Class to sync data across domain boundaries
 *
 * On the border between domains, there will be nodes that are shared between
 * two or more processes. The purpose of this class is to sync the data in
 * these nodes. This is done by using PETSc functions to scatter the data
 * to a global array where other processes can grab the relevant data. Only
 * the data that needs to be shared is sent.
 *
 * Prior to sending data, the communicators need to be set up. These
 * communicators can be reused. Once created, the data is shared by having
 * each process store its data in the class's send buffer, do the send, then
 * pull the data from the class's receive buffer. Example usage is:
 *
 * CrossDomainCommunicator comm;
 * comm.Initialize(p_grid);  // set up the data structures
 *
 * // store the data to send
 * const int n_send_nodes = comm.n_send_nodes();
 * const int *send_nodes = comm.send_node_indices();  // nodes to get data from
 * double *send_data = comm.send_data();  // the send buffer to put data in
 * for (int i = 0; i < n_send_nodes; i++) {
 *   int node_id = send_nodes[i];  // get the node index
 *   send_data[i] = node_data_array[node_id].data_to_sync;  // store data
 * }
 *
 * comm.SendRecvData();  // exchange the data
 *
 * // copy the received data
 * const int n_recv_nodes = comm.n_recv_nodes();
 * const int *recv_nodes = comm.recv_node_indices();  // nodes to put data in
 * const double *recv_data = comm.recv_data();  // receive buffer with data
 * for (int i = 0; i < n_recv_nodes; i++) {
 *   int node_id = recv_nodes[i];  // get the node index
 *   node_data_array[node_id].data_to_sync = recv_data[i];  // update data
 * }
 *
 */
class CrossDomainCommunicator {
 public:
  CrossDomainCommunicator() : initialized_(false) { }

  ~CrossDomainCommunicator() {
    CleanCommData();
  }

  /**
   * Initialize the object
   *
   * This function creates the data structures used to send and receive data
   * during the communication and also the PETSc communication structures.
   *
   * @param p_grid pointer to grid structure
   */
  void Initialize(GRID* p_grid);

  /**
   * Exchange data across processes
   *
   * This calls the PETSc communication routines to exchange the data.
   * Prior to calling this function, the relevant data must have been placed in
   * the appropriate location in the send_data_ array.
   * After calling this function, the desired data will be available in the
   * recv_data_ array. See the class comment for an example.
   */
  void SendRecvData();

  /**
   * Returns pointer to array of the indices of nodes whose data will be sent
   *
   * @return node indices of send nodes
   */
  inline const int* send_node_indices() const {
    return (*send_node_indices_).data();
  }

  /**
   * Returns pointer to array of the indices of nodes with data to recieve
   *
   * @return node indices of receive nodes
   */
  inline const PetscInt* recv_node_indices() const {
    return recv_node_indices_.data();
  }

  /**
   * Returns number of nodes whose data will be sent
   *
   * @return number of nodes with data to send
   */
  inline int n_send_nodes() const {
    return (*send_node_indices_).size();
  }

  /**
   * Returns number of nodes with data to receive
   *
   * @return number of nodes with data to receive
   */
  inline int n_recv_nodes() const {
    return recv_node_indices_.size();
  }

  /**
   * Returns pointer to send buffer
   *
   * @return pointer to send buffer
   */
  inline double* send_data() const {
    return send_data_.data();
  }

  /**
   * Returns pointer to receive buffer
   *
   * @return pointer to receive buffer
   */
  inline const double* recv_data() const {
    return recv_data_.data();
  }

 private:
  /**
   * Frees memory used for sending and receiving data
   */
  void CleanCommData();

  // for sharing data
  ZEROARRAY<int> *send_node_indices_;  ///< indices of nodes to send from
  ZEROARRAY<double> send_data_;  ///< array to store data in send_data_vec_
  Vec send_data_vec_;  ///< vec to store data to send
  VecScatter scatter_send_;  ///< scatter context to send data
  // for getting data
  ZEROARRAY<PetscInt> recv_node_indices_;  ///< indices of nodes to receive to
  ZEROARRAY<double> recv_data_;  ///< array to store data in recv_data_vec_
  Vec recv_data_vec_;  ///< vec to store data to send
  VecScatter scatter_recv_;  ///< scatter context to receive data

  Vec transfer_vec_;  ///< Vec used for transferring data during Communication

  bool initialized_;  ///< whether communication data is initialized
};

}  // namespace TALYFEMLIB

#endif  // DOMAINDECOMPOSITION_CROSS_DOMAIN_COMMUNICATOR_H_
