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

#include <integrator_test_node_data.h>


/**
 * Class stores the grid field data for the test system.
 */
class IntegratorTestGridField : public GridField<IntegratorTestNodeData> {
 public:
  /**
   * Sets values for all nodes.
   *
   * value 0 of each node is product[pi * sin(pi * x_i)] for i in [1 to nsd]
   * value 1 of each node is product[pi * cos(pi * x_i)] for i in [1 to nsd]
   */
  void SetValues() {
    const int nsd = p_grid_->nsd();
    for (int node_id = 0, n_nodes = p_grid_->n_nodes();
         node_id < n_nodes; node_id++) {
      const NODE* p_node = p_grid_->GetNode(node_id);  // node we care about
      IntegratorTestNodeData* p_data = &(GetNodeData(node_id));

      double x = p_node->x();
      double sin_field = M_PI * sin(M_PI*x);
      double cos_field = M_PI * cos(M_PI*x);
      if (nsd > 1) {
        double y = p_node->y();
        sin_field *= M_PI * sin(M_PI*y);
        cos_field *= M_PI * cos(M_PI*y);
      }
      if (nsd > 2) {
        double z = p_node->z();
        sin_field *= M_PI * sin(M_PI*z);
        cos_field *= M_PI * cos(M_PI*z);
      }

      p_data->value(0) = sin_field;
      p_data->value(1) = cos_field;
    }
  }
};
