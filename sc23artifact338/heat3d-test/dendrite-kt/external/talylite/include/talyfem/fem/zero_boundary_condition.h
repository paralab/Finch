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
#ifndef FEM_ZERO_BOUNDARY_CONDITION_H_
#define FEM_ZERO_BOUNDARY_CONDITION_H_


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
 * Zero boundary condition
 *
 * This is an essential boundary condition enforcing a solution value of zero.
 *
 * The purpose of this class is disable solving of a matrix row in support of
 * a different boundary condition. It is intended for internal library use.
 * To get the equivilent effect in user code, use ESSCondition instead.
 *
 * The condition is applied to the Ax=b equation by setting the matrix value at
 * the desired point to a value of 1.0. The corresponding entry in the b
 * vector is set to 0.0, yielding a trivial solution of 0.0 for the x value.
 */
class ZeroBoundaryCondition : public BaseBoundaryCondition {
 public:
  /**
   * Initializes the zero boundary condition.
   *
   * The condition is specified by giving the location of an entry in the
   * matrix.
   *
   * @param index the index giving the location of where to apply the bounday
   *              condition
   */
  explicit ZeroBoundaryCondition(LocalVarIdx index)
      : BaseBoundaryCondition(index) { }

  /**
   * Apply this boundary condition to Global matrix Ag and vector bg
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
   * Apply this boundary condition to the solution vector
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
                                  const ZeroBoundaryCondition& obj) {
    out << "x[" << obj.index_ << "]" << "= 0" << " loc";
    return out;
  }
};

}  // namespace TALYFEMLIB

#endif  // FEM_ZERO_BOUNDARY_CONDITION_H_
