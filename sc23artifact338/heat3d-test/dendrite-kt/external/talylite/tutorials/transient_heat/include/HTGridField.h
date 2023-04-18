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
#ifndef INCLUDE_HTGRIDFIELD_H_
#define INCLUDE_HTGRIDFIELD_H_

#include <iostream>
#include <fstream>

#include "HTAnalyticSolution.h"
#include "HTNodeData.h"


/**
 * Class stores the grid field data for the transient heat system.
 */
class HTGridField : public GridField<HTNodeData> {
 public:
  HTGridField(HTAnalyticSolution *analytic_solution)
      : analytic_solution_(analytic_solution) { }

  virtual ~HTGridField() { }

  /**
   * Sets boundary indicators for the system
   *
   * This is used when we load in a grid from e.g. a PLT file,
   * which doesn't come with boundary indicator information.
   * This is only needed is inputData.ifBoxGrid is set to 0 and the mesh is read
   * from a file. This assumes the system is a box grid with sides parallel to
   * the xyz axis.
   *
   * NOTE: This implementation will not work with domain decomposition
   */
  void SetBoundaries() {
    // The edges are determined by the maximum and minimum cooridinate values
    // which we're storing as points. As a starting point, we choose the first
    // node point in the system.
    ZEROPTV min_pos, max_pos;
    min_pos(0) = p_grid_->GetNode(0)->getCoor(0);
    min_pos(1) = p_grid_->GetNode(0)->getCoor(1);
    min_pos(2) = p_grid_->GetNode(0)->getCoor(2);
    max_pos(0) = p_grid_->GetNode(0)->getCoor(0);
    max_pos(1) = p_grid_->GetNode(0)->getCoor(1);
    max_pos(2) = p_grid_->GetNode(0)->getCoor(2);

    // Next we loop over all nodes in the system and see if any of the position
    // values are greater the current maximum or less than the current minimum
    // values. If so we store those new values. When this loop is finished, we
    // will have stored the minimum and maximum values for all dimensions.
    // This loop starts at one because we've already stored the values for the
    // zero-th node in the point data above.
    for (int node_id = 1; node_id < p_grid_->n_nodes(); node_id++) {
      for (int n = 0; n < p_grid_->nsd(); n++) {
        double len = p_grid_->GetCoord(node_id, n);
        if (len < min_pos(n)) { min_pos(n) = len; }
        if (len > max_pos(n)) { max_pos(n) = len; }
      }
    }

    // Now that we know the minimum and maximum values for each dimension, we
    // loop over the nodes and test each coordinate to see if it lies on a
    // edge. If so, we add the appropriate boundary indicator to the node.
    // The indicators are numbers which depend on the number of dimensions of
    // the system.
    for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      for (int n = 0; n < p_grid_->nsd(); n++) {
        // Check if this is on the lower edges (i.e. the coordinate is a equal
        // to the minimum value). If so, set the indicator based on which
        // direction this is: x=1, y=2, z=3
        if (p_grid_->GetCoord(node_id, n) <= min_pos(n)) {
          p_grid_->GetNode(node_id)->addIndicatorNum(n + 1);
        }
        // Check if this is on the lower edges (i.e. the coordinate is a equal
        // to the minimum value). If so, set the indicator based on which
        // direction this is:
        // x=1+n_dimensions, y=2+n_dimensions, z=3+n_dimensions
        if (p_grid_->GetCoord(node_id, n) >= max_pos(n)) {
          p_grid_->GetNode(node_id)->addIndicatorNum(n + 1 + p_grid_->nsd());
        }
      }
    }

    // generate surfaces
    p_grid_->SetCaredSurfaceIndicator();
    for (int elem_id = 0; elem_id < p_grid_->n_elements(); elem_id++) {
      p_grid_->GetElm(elem_id)->GenSurfaceIndicator(p_grid_, p_grid_->cared_surface_indicator_);
    }
  }

  /**
   * Sets the initial condition for all node points
   *
   * This calculates the initial condition for every node point and stores
   * those values in the node data objects.
   */
  void SetInitialConditions() {
    for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      HTNodeData* pData = &(GetNodeData(node_id));  // node data for this node
      double t_zero_value = analytic_solution_->ValueAt(node_id);
      pData->u = pData->u_pre = pData->u_analytical = t_zero_value;

      // initialize hermite derivatives to 0
      for (int i = 0; i < 7; i++) {
        pData->du[i] = pData->du_pre[i] = 0;
      }
    }
    PrintStatus("Initial conditions set ");
  }

  /**
   * Sets the analytic solution for all node points
   *
   * This calculates the analytic solution at each node point and stores the
   * value in the node data objects.
   *
   * @param t the current time
   * @param n_dimensions number of dimensions
   */
  void SetAnalyticalSolution(double t, int n_dimensions) {
    for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      HTNodeData* pData = &(GetNodeData(node_id));  // node data for this node
      pData->u_analytical = analytic_solution_->ValueAt(node_id, t);
    }
  }

  /**
   * Calculate and return the maximum node error for the system.
   *
   * The error in this case is defined as the nodal value of:
   *  |u - u_exact|
   *  -------------
   *     u_exact
   *
   * @param if_domain_decomp whether the system uses domain decomposition
   * @return the maximum nodal error in the system
   */
  double CalcMaxError(bool if_domain_decomp) {
    double local_error = 0.0;  // local to this MPI process
    double global_error = 0.0;  // global across all MPI processes

    // loop over all nodes on this process and find the maximum error
    for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      HTNodeData* pData = &(GetNodeData(node_id));
      double error = fabs((pData->u - pData->u_analytical) / pData->u_analytical);
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

  double CalcL2Error(const InputData* input_data, double t) const {
    FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION | BASIS_DIMENSION_REDUCTION);

    double l2_error = 0.0;
    const double n_elements = p_grid_->n_elements();
    for (int elm_id = 0; elm_id < n_elements; elm_id++) {
      if (!p_grid_->IsMyElement(elm_id))
        continue;

      fe.refill(elm_id, 0);
      while (fe.next_itg_pt()) {
        const double detJxW = fe.detJxW();
        double val_c = valueFEM(fe, U);
        double val_a = analytic_solution_->ValueAt(fe.position(), t);
        l2_error += (val_c - val_a) * (val_c - val_a) * detJxW;
      }
    }

    double global_error;
    MPI_Allreduce(&l2_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    return sqrt(global_error);
  }

 private:
  HTAnalyticSolution *analytic_solution_;
};

#endif  // INCLUDE_HTGRIDFIELD_H_
