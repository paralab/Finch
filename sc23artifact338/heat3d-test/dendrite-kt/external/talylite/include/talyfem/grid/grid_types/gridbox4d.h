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

#ifndef TALYFEM_GRIDBOX4D_H
#define TALYFEM_GRIDBOX4D_H
#include <string>

#include <talyfem/grid/grid_types/grid.h>  // parent class


namespace TALYFEMLIB {

/**
 * Box grid in 4d
 */
class GRIDBox4D : public GRID {
 public:
  /**
   * Create a GRIDBox4D object
   *
   * @param basis_function_order: order of basis functions to use
   */
  explicit GRIDBox4D(int basis_function_order = 1);

  ~GRIDBox4D() { }

 private:

  std::string GridTypeName() const { return std::string("GRIDBOX4D"); }

  virtual void redimDD(const double* dimensions, const int* n_elems);
};

}  // namespace TALYFEMLIB
#endif //TALYFEM_GRIDBOX4D_H
