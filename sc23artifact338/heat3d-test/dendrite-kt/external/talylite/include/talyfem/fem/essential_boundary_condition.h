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
#ifndef FEM_ESSENTIAL_BOUNDARY_CONDITION_H_
#define FEM_ESSENTIAL_BOUNDARY_CONDITION_H_


#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/fem/base_boundary_condition.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

/**
 * Essential boundary condition
 *
 * This class is responsible for storing and applying an essential
 * boundary condition to the system.
 *
 * The condition is applied to the matrix solution by setting the value at
 * the desired point to a value of 1 (this can be changed if needed for
 * stability). This will keep them from being solved by the solver.
 * The condition is applied to the solution vector by replacing the entry
 * with the actual value at that point. This value will not be changed during
 * the solution process.
 */
class EssentialBoundaryCondition : public BaseBoundaryCondition {
 public:
  /**
   * Initializes the object to the given essential boundary condition
   *
   * The condition is specified by giving a value and a location to put it.
   * The location is the index of an entry in the matrix. The value is the
   * value of the solution for that entry. The matrix coefficient can also be
   * specified. This is the value to place in the Ag Matrix. Typically this is
   * one; However, this can apparently lead to the matrix being ill-conditioned.
   * To avoid this, a different matrix coefficient value can be given. In this
   * case the value of the essential boundary condition that is specified should
   * be matrix_coefficent * essential_boundary_value.
   * With the normal case where the matrix value is 1, the boundary value that
   * should be specified is simply equal to the essential boundary value.
   *
   * @param index the index of involved variables start (index starts from 1,
   *              But this array is a C array, which starts from 0 )
   * @param matrix_value coefficient to place in the matrix entry
   * @param value value of essential_condition * matrix_value (for cases where
   *              matrix_value = 1, this will be equal to the value of essential
   *              boundary condition)
   */
  EssentialBoundaryCondition(LocalVarIdx index, double matrix_value = 1.0,
                             double value = 0.0);

  ~EssentialBoundaryCondition() { }

  /**
   * Apply this essential boundary condition to Global matrix Ag and vector bg
   *
   * @param Ag Global stiffness matrix (Petsc format)
   * @param bg Global load vector (Petsc format)
   * @param pXMap the global numbering of each variable (member variable
   *              index is using a local numbering, which may different
   *              from global numbering)
   * @param recalc_matrix whether the matrix needs to be recalculated
   */
  void Apply(Mat& Ag, Vec& bg, const ZEROARRAY<GlobalVarIdx>* pXMap = NULL,
             bool recalc_matrix = true) const;

  /**
   * Apply this essential boundary condition to the solution vector
   *
   * @param solution the solution vector to apply the condition to
   */
  void ApplyEssBCToSolution(ZEROARRAY<double>& solution) const;

  /**
   * Pass object to stream
   *
   * @param out output stream to send object to
   * @param obj object to output to stream
   * @return modified output stream
   */
  friend std::ostream& operator<<(std::ostream& out,
                                  const EssentialBoundaryCondition& obj) {
    out << obj.matrix_value_ << "*x[" << obj.index_ << "]";
    out << "=" << obj.essential_value_ << " loc";
    return out;
  }

 private:
  double essential_value_;  ///< the value to put in the bg vector
  double matrix_value_;  ///< coefficient to put in the Ag matrix
};

}  // namespace TALYFEMLIB

#endif  // FEM_ESSENTIAL_BOUNDARY_CONDITION_H_
