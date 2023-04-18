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
#ifndef TESTS_PERIODICEXCHANGE_INCLUDE_PBC_EQUATION_H_
#define TESTS_PERIODICEXCHANGE_INCLUDE_PBC_EQUATION_H_

#include "include/pbc_input_data.h"
#include "include/pbc_node_data.h"


class PBCEquation : public CEquation<PBCNodeData> {
 public:
  PBCInputData* input_data_;

  explicit PBCEquation(PBCInputData* input_data)
      : input_data_(input_data) {
  }

  virtual void fillEssBC() {
    // Set node values based on boundary conditions.
    // the values are set to the following:
    // * periodic source nodes are set to the node number * the degree
    //   of freedom
    // * periodic target nodes are set to the negative of the node number * the
    //   degree of freedom
    // * nodes with direchlet boundaries are set to 999000 + the node number
    //
    // All node numbers are the global values.
    // Where more than one boundary condition applies, the value is set based
    // on the follworing priority (highest to lowest):
    // periodic source, periodic target, direchlet condition.
    for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
      bool is_direchlet = false;
      bool is_periodic_target = false;
      bool is_periodic_source = false;
      for (int k = 0; k < input_data_->nsd; k++) {
        if (input_data_->boundaryType[k] == "direchlet") {
          if (p_grid_->BoNode(nodeID, 2 * k + 1)
              || p_grid_->BoNode(nodeID, 2 * k + 2)) {
            is_direchlet = true;
          }
        } else if (input_data_->boundaryType[k] == "periodic") {
          if (p_grid_->BoNode(nodeID, 2 * k + 1)) {
            is_periodic_source = true;
          }
          if (p_grid_->BoNode(nodeID, 2 * k + 2)) {
            is_periodic_target = true;
          }
        }
        if (is_periodic_source) {
          // PrintInfo("periodic source: ", p_grid_->physical_map(nodeID));

          double base_value;
          if (p_grid_->parallel_type_ == kWithDomainDecomp) {
            base_value = 1.0 * (this->p_grid_->physical_map(nodeID) + 1);
          } else {
            base_value = 1.0 * (nodeID + 1);
          }
          for (int i_dof = 0; i_dof < n_dof_; i_dof++) {
            const double value = (i_dof + 1) * base_value;
            this->p_data_->GetNodeData(nodeID).u[i_dof] = value;
          }
        } else if (is_periodic_target) {
          // PrintInfo("periodic target: ", p_grid_->physical_map(nodeID));

          double base_value;
          if (p_grid_->parallel_type_ == kWithDomainDecomp) {
            base_value = -1.0 * (this->p_grid_->physical_map(nodeID) + 1);
          } else {
            base_value = -1.0 * (nodeID + 1);
          }
          for (int i_dof = 0; i_dof < n_dof_; i_dof++) {
            const double value = (i_dof + 1) * base_value;
            this->p_data_->GetNodeData(nodeID).u[i_dof] = value;
          }
        } else if (is_direchlet) {
          // PrintInfo("direchlet: ", p_grid_->physical_map(nodeID));

          double node_value;
          if (p_grid_->parallel_type_ == kWithDomainDecomp) {
            node_value = 999000 + (this->p_grid_->physical_map(nodeID) + 1);
          } else {
            node_value = 999000 + (nodeID + 1);
          }
          for (int i_dof = 0; i_dof < n_dof_; i_dof++) {
            this->p_data_->GetNodeData(nodeID).u[i_dof] = node_value;
          }
        }
      }
    }
  }

  virtual void Solve(double dt = 0.0, double t = 0.0) {
    SetInitialConditions();  // set all variables to zero
    initEssBC();
    fillEssBC();

    ApplyEssBCToSolution();
    Assemble();
    ApplyAllBC();

    // the exchange of values is done by the following lines. Variables one
    // and three in the node data are exchanged across periodic boundaries.
    // The result is that the values are the same for a node and its periodic
    // partner.
    // Value two in the node data is not changed in this example.
    //
    // There are two methods of doing the exchange.
    // The first method loops over the periodic source nodes and stores the data
    // we want to send in pbc_source_data_. In this case, we are storing the
    // first value in the node data.
    // next we call DoPBCExchange() to do the exchange. This will put the data
    // for periodic nodes in the pbc_target_data_.
    // finally we loop over the periodic target nodes and put the values we
    // received in the desired location. In this case, we put them in the
    // first values in the node data.
    // The result of this is that the first data value of the node data
    // has been synced across the periodic boundaries. The value we sync across
    // the boundary does not need to be in the node data structure.
    if (exchanger_ != NULL) {
      PetscScalar *send_buffer = exchanger_->send_buffer();
      const int *send_indices = exchanger_->send_indices();
      for (int i = 0; i < exchanger_->n_send_values(); i++) {
        int node_id = send_indices[i];  // index of ith source node
        // set pbc_source_data_[i] to the data value we want to exchange
        send_buffer[i] = p_data_->GetNodeData(node_id).u[0];
      }

      DoPBCExchange();

      PetscScalar *recv_buffer = exchanger_->recv_buffer();
      const int *recv_indices = exchanger_->recv_indices();
      for (int i = 0; i < exchanger_->n_recv_values(); i++) {
        int node_id = recv_indices[i];  // index of ith target node
        // put recv_buffer[i] in the location we want to store it
        p_data_->GetNodeData(node_id).u[0] = recv_buffer[i];
      }

      // The second method is a short hand for when the value we want to
      // exchange is in the node data structure. DoPBCExchangeByID(0) will
      // exchange the ith data item in the node data. The value of i starts at
      // zero. So in this case, we are exchanging the third item in the node
      // data structure.
      DoPBCExchangeByID(2);
    }
  }

  // set values to zero everywhere
  void SetInitialConditions() {
    for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
      for (int i_dof = 0; i_dof < n_dof_; i_dof++) {
        this->p_data_->GetNodeData(nodeID).u[i_dof] = 0.0;
      }
    }
  }

  void Integrands(const FEMElm& fe, ZeroMatrix<double>& Ae,
                  ZEROARRAY<double>& be) {}
};

#endif  // TESTS_PERIODICEXCHANGE_INCLUDE_PBC_EQUATION_H_
