/*
  Copyright 2014-2017 Baskar Ganapathysubramanian

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
#ifndef BT_GRIDFIELD_HPP
#define BT_GRIDFIELD_HPP

#include <iostream>
#include <fstream>
#include "BTNodeData.h"
#include "BTInputData.h"

class BTGridField : public GridField<BTNodeData> {
 public:
  BTGridField() { }
  virtual ~BTGridField() { }

  void SetIC(int nsd) {
    for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
      BTNodeData* pData = &(GetNodeData(nodeID));

      // initialize u and du (for hermite case) to 0
      for (int i = 0; i < 8; i++) {
        pData->u[i] = 0.0;
      }
    }
    PrintStatusStream(std::cerr, "", "IC set ");
  }

  void SetAnalyticalSolution() {
    for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
      const NODE* node = this->p_grid_->GetNode(nodeID);
      this->GetNodeData(nodeID).u_a = GetAnalyticalSolution(node->location());
    }
  }

  double GetAnalyticalSolution(const ZEROPTV& pt) const {
    return -log(1 + cos(M_PI * (0.5 + pt.x())));
  }

  /**
   * Calculate and return the maximum node error for the system.
   *
   * The error in this case is defined as the nodal value of:
   *  |u - u_exact|
   *
   * @param if_domain_decomp whether the system uses domain decomposition
   * @return the maximum nodal error in the system
   */
  double CalcMaxError(bool if_domain_decomp) {
    double local_error = 0.0;  // local to this MPI process
    double global_error = 0.0;  // global across all MPI processes

    // loop over all nodes on this process and find the maximum error
    for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      const BTNodeData* pData = &(GetNodeData(node_id));
      const double error = fabs(pData->u[0] - pData->u_a);
      if (error > local_error) { local_error = error; }
    }

    // for meshes with domain decomposition, find maximum global error
    // across all processes using MPI.
    if (if_domain_decomp) {
      MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      local_error = global_error;
    } else {
      global_error = local_error;
    }
    return global_error;
  }

  double CalcL2Error(const BTInputData* input_data) {
    FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION | BASIS_DIMENSION_REDUCTION);

    double l2_error = 0.0;
    const double n_elements = p_grid_->n_elements();
    for (int elm_id = 0; elm_id < n_elements; elm_id++) {
      if (!p_grid_->IsMyElement(elm_id))
        continue;

      fe.refill(elm_id, 0);
      while (fe.next_itg_pt()) {
        const double detJxW = fe.detJxW();
        double val_c = valueFEM(fe, 0);
        double val_a = GetAnalyticalSolution(fe.position());
        l2_error += (val_c - val_a) * (val_c - val_a) * detJxW;
      }
    }

    double global_error = 0.0;
    MPI_Allreduce(&l2_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    return sqrt(global_error);
  }
};

#endif
