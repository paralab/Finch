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
#ifndef FEM_BASE_BOUNDARY_CONDITION_H_
#define FEM_BASE_BOUNDARY_CONDITION_H_

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

/**
 * Base boundary condition
 *
 * Base class for boundary conditions. This provides some common funtions.
 * It is not meant to be used directly.
 */
class BaseBoundaryCondition {
 public:
  /**
   * Initializes the boundary condition.
   *
   * @param index the index giving the location of where to apply the boundary
   *              condition
   */
  explicit BaseBoundaryCondition(LocalVarIdx index) : index_(index) { }

  /**
   * Returns the index of the boundary condition as specified during creation
   *
   * @return index of boundary condition
   */
  inline LocalVarIdx get_index() const { return index_; }

 protected:
  /**
   * Returns the index of the boundary condition in the Petsc structures.
   *
   * This uses the passed mapping if needed.
   *
   * @param pXMap the global numbering of each variable (member variable
   *              index is using a local numbering, which may different
   *              from global numbering)
   * @return index location of the boundary condition in the Petsc structures
   */
  GlobalVarIdx PetscIndex(const ZEROARRAY<GlobalVarIdx>* pXMap) const;

  LocalVarIdx index_;  ///< row index of the boundary condition
};

}  // namespace TALYFEMLIB

#endif  // FEM_BASE_BOUNDARY_CONDITION_H_
