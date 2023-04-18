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
#pragma once

#include <vector>

#include <talyfem/grid/elem.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>

namespace TALYFEMLIB {

/**
 * Accessor for accessing the nodes of a particular element.
 * Simplifies passing around an (ELEM*, GRID*) pair all the time.
 */
struct ElemNodes {
 public:
  inline ElemNodes() {}

  /**
   * @param elem element
   * @param grid grid containing node data for elem
   */
  inline ElemNodes(const ELEM* elem, const GRID* grid) {
    const int number_nodes = elem->n_nodes();
    nodes_.resize(number_nodes);
    for (int i = 0; i < number_nodes; i++) {
      nodes_[i] = grid->node_array_[elem->node_id_array(i)]->location();
    }
  }

  /**
   * @param nodes node locations
   * @param number_nodes number of node locations
   */
  inline ElemNodes(const ZEROPTV* nodes, int number_nodes) {
    nodes_.resize(number_nodes);
    for (int i = 0; i < number_nodes; i++) {
      nodes_[i] = nodes[i];
    }
  }

  /**
   * Get the location of node i.
   * @param i element-local node index
   * @returns node i
   */
  inline const ZEROPTV& node_pos(int i) const {
    return nodes_.at(i);
  }

  /**
   * @returns number of nodes
   */
  inline unsigned int n_nodes() const {
    return nodes_.size();
  }

 private:
  std::vector<ZEROPTV> nodes_;
};

}  // namespace TALYFEMLIB
