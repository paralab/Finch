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
#ifndef GRID_GRID_TYPES_GRID2D_H_
#define GRID_GRID_TYPES_GRID2D_H_

#include <string>

#include <talyfem/grid/grid_types/grid.h>  // parent class


namespace TALYFEMLIB {

class InputData;

/**
 * Simple box grid in 2d with triangle elements.
 *
 * The redim function generates a normal 2D box grid with each box split
 * into two triangles:
 *   +--+       +--+
 *   |  | ==>   |\ |
 *   |  |       | \|
 *   +--+       +--+
 * Only basis order 1 is currently supported
 */
class GRID2D : public GRID {
 public:
  /**
   * Create a GRID2D object
   *
   * @param basis_function_order order of basis functions to use
   */
  explicit GRID2D(int basis_function_order = 1);

  ~GRID2D() { }

 private:
  void CreateElementsBasis1();
  void CreateNodes(const double *dimensions);

  std::string GridTypeName() const { return std::string("GRID2D"); }

  void SetNodeElemCounts(const int* n_elems, bool set_totals = true);
};

}  // namespace TALYFEMLIB

#endif  // GRID_GRID_TYPES_GRID2D_H_
