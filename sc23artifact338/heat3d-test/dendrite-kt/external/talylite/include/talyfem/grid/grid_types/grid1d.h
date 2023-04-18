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
#ifndef GRID_GRID_TYPES_GRID1D_H_
#define GRID_GRID_TYPES_GRID1D_H_

#include <string>

#include <talyfem/grid/grid_types/grid.h>  // parent class


namespace TALYFEMLIB {

/**
 * Simple grid in 1D
 *
 * redim will create a 1D linear grid with the given number of elements. Note
 * that actual number of elements created is the requested number divided
 * by the order of the basis function. The number of elements requested in each
 * direction must be a multiple of the order of the basis function.
 */
class GRID1D : public GRID {
 public:
  /**
   * Create a GRID1D object
   *
   * @param basis_function_order order of basis functions to use
   */
  explicit GRID1D(int basis_function_order = 1);

  ~GRID1D() { }

 private:
  /**
   * Generates elements for an arbitrary basis function order
   *
   * @param basis order of the basis function
   */
  void CreateElementsBasisGeneral(int basis);

  void CreateElementsBasis1() { CreateElementsBasisGeneral(1); }
  void CreateElementsBasis2() { CreateElementsBasisGeneral(2); }
  void CreateElementsBasis3() { CreateElementsBasisGeneral(3); }

  void CreateNodes(const double *dimensions);

  std::string GridTypeName() const { return std::string("GRID1D"); }
};

}  // namespace TALYFEMLIB

#endif  // GRID_GRID_TYPES_GRID1D_H_
