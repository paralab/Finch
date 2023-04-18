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
#ifndef TESTS_PERIODICNODE_INCLUDE_PBC_EQUATION_H_
#define TESTS_PERIODICNODE_INCLUDE_PBC_EQUATION_H_

#include "include/pbc_input_data.h"
#include "include/pbc_node_data.h"
#include "FEM/preallocator_perfect.h"

class PBCEquation : public CEquation<PBCNodeData> {
 public:
  PBCInputData* input_data_;

  explicit PBCEquation(PBCInputData* input_data)
      : input_data_(input_data) {
    SetPreallocator(new TALYFEMLIB::PreallocatorPerfect<PBCNodeData>(this));
  }

  virtual void fillEssBC() {
    for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
      for (int k = 0; k < input_data_->nsd; k++) {
        if (input_data_->boundaryType[k] == "direchlet") {
          // index of boundary nodes (BoNode values) starts at 1, while the
          // index of boundary values starts at zero.
          if (p_grid_->BoNode(nodeID, 2 * k + 1)) {
            specifyValue(nodeID, 0, 0.0);
          }
          if (p_grid_->BoNode(nodeID, 2 * k + 2)) {
            specifyValue(nodeID, 0, 0.0);
          }
        }
      }
    }
  }

  virtual void Solve(double dt = 0.0, double t = 0.0) {
    initEssBC();
    fillEssBC();

    ApplyEssBCToSolution();
    Assemble();
    ApplyAllBC();

    FillNodeData();
  }

  // put the index of the node, after applying periodic boundary conditions,
  // into the node data.
  void FillNodeData() {
    for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
      if (has_per_bc_.get(nodeID)) {
        int partner = periodic_bounds_->GetPeriodicPartner(nodeID) + 1;
        this->p_data_->GetNodeData(nodeID).u[0] = partner;
      } else {
        if (p_grid_->parallel_type_ == kWithDomainDecomp) {
          this->p_data_->GetNodeData(nodeID).u[0] =
              p_grid_->physical_map(nodeID) + 1;
        } else {
          this->p_data_->GetNodeData(nodeID).u[0] = nodeID + 1;
        }
      }
    }
  }

  void Integrands(const FEMElm& fe, ZeroMatrix<double>& Ae,
                  ZEROARRAY<double>& be) {}

  void IntegrandsPreallocator(const FEMElm& fe, ZeroMatrix<bool>& Ae) {
    for (int i = 0; i < Ae.nx(); i++) {
      for (int j = 0; j < Ae.ny(); j++) {
        Ae(i, j) = true;
      }
    }
  }
};

#endif  // TESTS_PERIODICNODE_INCLUDE_PBC_EQUATION_H_
