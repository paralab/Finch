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
#ifndef FEM_BOUNDARY_CONDITIONS_H_
#define FEM_BOUNDARY_CONDITIONS_H_

#include <vector>

#include <talyfem/fem/essential_boundary_condition.h>
#include <talyfem/fem/zero_boundary_condition.h>
#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/utils/utils.h>  // for Print* functions


namespace TALYFEMLIB {

/**
 * Boundary Condition Container
 *
 * This class is responsible for storing and applying all boundary condtions.
 */
class BoundaryConditions {
 public:
  BoundaryConditions() { }

  ~BoundaryConditions() {
    DeleteAllConditions();
  }

  /**
   * Deletes all boundary condition objects stored in this object
   */
  void DeleteAllConditions() {
    DeleteDirichletConditions();
    DeletePeriodicConditions();
  }

  /**
   * Deletes dirichlet boundary condition objects stored in this object
   */
  void DeleteDirichletConditions() {
    for (dirichlet_vect_type::iterator it = dirichlet_boundaries_.begin();
         it != dirichlet_boundaries_.end(); it++)
      delete *it;
    dirichlet_boundaries_.clear();
  }

  /**
   * Deletes periodic boundary condition objects stored in this object
   */
  void DeletePeriodicConditions() {
    for (periodic_vect_type::iterator it = periodic_boundaries_.begin();
         it != periodic_boundaries_.end(); it++)
      delete *it;
    periodic_boundaries_.clear();
  }

  /**
   * Prints all boundary conditions stored in this object
   */
  void PrintAllConditions() const {
    PrintDirichletConditions();
    PrintPeriodicConditions();
  }

  /**
   * Prints dirichlet boundary conditions stored in this object
   */
  void PrintDirichletConditions() const {
    for (dirichlet_vect_type::const_iterator it = dirichlet_boundaries_.begin();
         it != dirichlet_boundaries_.end(); it++) {
      PrintInfo(*(*it));
    }
  }

  /**
   * Prints periodic boundary conditions stored in this object
   */
  void PrintPeriodicConditions() const {
    for (periodic_vect_type::const_iterator it = periodic_boundaries_.begin();
         it != periodic_boundaries_.end(); it++) {
      PrintInfo(*(*it));
    }
  }

  /**
   * Apply dirichlet boundary conditions to global Ag matrix and bg vector
   *
   * @param Ag Global stiffness matrix (Petsc format)
   * @param bg Global load vector (Petsc format)
   * @param pXMap the global numbering of each variable (member variable
   *              index is using a local numbering, which may different
   *              from global numbering)
   * @param recalc_matrix whether the matrix needs to be recalculated
   */
  void ApplyDirichletConditions(Mat& Ag, Vec& bg,
                                const ZEROARRAY<GlobalVarIdx>* pXMap = NULL,
                                bool recalc_matrix = true) {
    for (dirichlet_vect_type::iterator it = dirichlet_boundaries_.begin();
         it != dirichlet_boundaries_.end(); it++) {
      (*it)->Apply(Ag, bg,  pXMap, recalc_matrix);
    }
  }

  /**
   * Apply periodic boundary conditions to global Ag matrix and bg vector
   *
   * Note that this does NOT copy values across periodic boundaries, it only
   * sets the matrix entries to prevent solving for periodic partners whose
   * values will be overwritten later.
   *
   * @param Ag Global stiffness matrix (Petsc format)
   * @param bg Global load vector (Petsc format)
   * @param pXMap the global numbering of each variable (member variable
   *              index is using a local numbering, which may different
   *              from global numbering)
   * @param recalc_matrix whether the matrix needs to be recalculated
   */
  void ApplyPeriodicConditions(Mat& Ag, Vec& bg,
                                const ZEROARRAY<GlobalVarIdx>* pXMap = NULL,
                                bool recalc_matrix = true) {
    for (periodic_vect_type::iterator it = periodic_boundaries_.begin();
         it != periodic_boundaries_.end(); it++) {
      (*it)->Apply(Ag, bg,  pXMap, recalc_matrix);
    }
  }

  /**
   * Apply this essential boundary condition to the solution vector
   *
   * @param solution the solution vector to apply the condition to
   */
  void ApplyDirichletConditionsToSolution(ZEROARRAY<double>& solution) {
    for (dirichlet_vect_type::iterator it = dirichlet_boundaries_.begin();
         it != dirichlet_boundaries_.end(); it++) {
      (*it)->ApplyEssBCToSolution(solution);
    }
  }

  /**
   * Apply periodic boundary conditions to the solution vector
   *
   * @param solution the solution vector to apply the condition to
   */
  void ApplyPeriodicConditionsToSolution(ZEROARRAY<double>& solution) {
    for (periodic_vect_type::iterator it = periodic_boundaries_.begin();
         it != periodic_boundaries_.end(); it++) {
      (*it)->ApplyEssBCToSolution(solution);
    }
  }

  /**
   * Adds a Dirichlet boundary condition
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
  void AddDirichletCondition(LocalVarIdx index, double matrix_value = 1.0,
                              double value = 0.0) {
    EssentialBoundaryCondition* pBC = new EssentialBoundaryCondition(
        index, matrix_value, value);
    dirichlet_boundaries_.push_back(pBC);
  }

  /**
   * Add an essential boundary condition for a periodic variable
   *
   * This condition will set the diagonal matrix entry to 1.0 and the vector
   * entry to 0.0, leading to a solution of 0.0. This is basically a fake
   * equation to prevent PETSc from solving for an non-existant point.
   *
   * @param index index in the periodic variable
   */
  void AddPeriodicCondition(LocalVarIdx index) {
    ZeroBoundaryCondition* pBC = new ZeroBoundaryCondition(index);
    periodic_boundaries_.push_back(pBC);
  }

 private:
  typedef std::vector<EssentialBoundaryCondition*> dirichlet_vect_type;
  ///< type of dirichlet  boundary list

  typedef std::vector<ZeroBoundaryCondition*> periodic_vect_type;
  ///< type of periodic boundary list

  dirichlet_vect_type dirichlet_boundaries_;  ///< Vector of all the dirichlet
                                              ///< boundary conditions
  periodic_vect_type periodic_boundaries_;  ///< Vector of all the periodic
                                            ///< boundary conditions
};

}  // namespace TALYFEMLIB

#endif  // FEM_BOUNDARY_CONDITIONS_H_
