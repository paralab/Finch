/*
  Copyright 2017 Baskar Ganapathysubramanian

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
#pragma once

#include <talyfem/integrator/function_integrator.h>
#include <talyfem/grid/gridfield.h>


namespace TALYFEMLIB {

/**
 * Class to calculate the integrate of a NodeData item over the system
 */
template<class NodeData>
class ValueFunction : public FunctionIntegrator {
 public:
  /**
   * Construct the value integrator object.
   *
   * @param[in] basis_function basis function to use
   * @param[in] rel_order relative order to use
   * @param[in] p_data pointer to grid field data
   * @param[in] index index of NodeData item to integrate
   * @param[in] do_accelerate whether to use acceleration
   * @param[in] method method of integration (see IntegrationMethod for options)
   */
  ValueFunction(kBasisFunction basis_function, int rel_order,
                GridField<NodeData> *p_data, int index,
                bool do_accelerate = false,
                AssemblyMethod method = kAssembleGaussPoints)
    : FunctionIntegrator(rel_order, do_accelerate, method),
      index_(index),
      p_data_(p_data) { }

  /**
   * Calculate the NodeData value at the Gauss point.
   *
   * @param[in] fe the finite element with the appropriate Gauss point selected
   * @return the value of the NodeData item at the gauss point
   */
  virtual double CalcGaussPointIntegral(FEMElm& fe) override {
    double value = 0.0;
    for (int a = 0, max_a = fe.nbf(); a < max_a; a++) {
      const int node_id = fe.elem()->ElemToLocalNodeID(a);
      const double local_value = p_data_->GetNodeData(node_id).value(index_);
      value += local_value * fe.N(a);
    }
    return value * fe.detJxW();
  }

  int index_;  ///< index of variable to integrate over
  GridField<NodeData> *p_data_;  ///< pointer to grid field data
};

}  // namespace TALYFEMLIB
